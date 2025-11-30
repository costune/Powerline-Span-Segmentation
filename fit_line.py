#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
电力线拟合脚本
从LAS点云数据中提取并拟合电力线，输出Shapefile格式结果
"""

import argparse
import math
import os
from pathlib import Path

import laspy
import numpy as np
from numpy import sqrt
from scipy.optimize import curve_fit
from sklearn.cluster import DBSCAN
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from tqdm import tqdm
import shapefile
from pyproj import CRS, Transformer


# ==================== 拟合函数 ====================

def fit_catenary_xz(x, z):
    """
    使用悬链线模型拟合 x-z 平面投影
    
    悬链线模型: z(x) = a*cosh((x-x0)/a) + c
    
    参数:
        x: ndarray, x坐标数组
        z: ndarray, z坐标数组
    
    返回:
        tuple: (a, x0, c) 悬链线模型参数
    """
    def catenary(x, a, x0, c):
        return a * np.cosh((x - x0) / a) + c

    # 初始参数猜测
    a0 = (np.max(z) - np.min(z)) / 2 + 1e-6  # 防止除零
    x0_0 = np.mean(x)
    c0 = np.mean(z) - a0
    p0 = [a0, x0_0, c0]

    # 非线性最小二乘拟合
    params, _ = curve_fit(catenary, x, z, p0=p0, maxfev=20000)
    a, x0, c = params
    
    return a, x0, c


def ransac_fit_line(x, y):
    """
    使用线性回归拟合直线
    
    参数:
        x: ndarray, shape (n, 1), x坐标
        y: ndarray, shape (n, 1), y坐标
    
    返回:
        tuple: (a, b) 直线参数，y = a*x + b
    """
    model = LinearRegression()
    model.fit(x, y)
    
    a = model.coef_
    b = model.intercept_
    
    return a, b


def filter_points_by_fit(x, y, z, a, b, catenary_params, threshold, n=1):
    """
    根据拟合结果过滤离群点
    
    参数:
        x, y, z: 坐标数组
        a, b: xy平面直线拟合参数
        catenary_params: xz平面悬链线拟合参数 (a, x0, c)
        threshold: 过滤阈值
        n: 阈值倍数
    
    返回:
        tuple: 过滤后的 (x, y, z, Y, Z)
    """
    # y方向约束
    Y = (a * x + b).flatten()
    idx_y = np.where(np.abs(y - Y) < n * threshold)
    x = x[idx_y]
    y = y[idx_y]
    z = z[idx_y]
    
    # z方向约束
    a_cat, x0_cat, c_cat = catenary_params
    Z = a_cat * np.cosh((x - x0_cat) / a_cat) + c_cat
    idx_z = np.where(np.abs(z - Z) < n * threshold)
    x = x[idx_z]
    y = y[idx_z]
    z = z[idx_z]
    
    # 重新计算拟合值
    Y = (a * x + b).flatten()
    Z = a_cat * np.cosh((x - x0_cat) / a_cat) + c_cat
    
    return x, y, z, Y, Z


# ==================== 点云处理函数 ====================

def get_tower_centers(tower_path, normal_tower_height):
    """
    从电力塔点云中提取电力塔中心坐标
    
    参数:
        tower_path: str, 电力塔点云文件路径
        normal_tower_height: float, 正常电力塔高度阈值
    
    返回:
        list: 电力塔中心坐标列表 [(x1, y1), (x2, y2), ...]
    """
    las = laspy.read(tower_path)
    las_x = np.array(las.x)
    las_y = np.array(las.y)
    las_z = np.array(las.z)
    data = np.stack([las_x, las_y, las_z], axis=1)

    # DBSCAN聚类
    dbscan = DBSCAN(eps=3, min_samples=15)
    labels = dbscan.fit_predict(data)

    # 筛选电力塔簇
    tower_clusters = []
    for label in np.unique(labels):
        if label == -1:  # 跳过噪声点
            continue
            
        indices = np.where(labels == label)
        cluster = data[indices]
        
        # 计算簇的尺寸
        height = np.max(cluster[:, 2]) - np.min(cluster[:, 2])
        width = np.max(cluster[:, 0]) - np.min(cluster[:, 0])
        length = np.max(cluster[:, 1]) - np.min(cluster[:, 1])

        # 根据高度特征识别电力塔
        # 110kV: 15-30m, 220kV: 30-40m, 500kV: 40-60m, 1000kV: 60-80+m
        if height > normal_tower_height and height > width and height > length:
            tower_clusters.append(cluster)

    # 计算每个电力塔的中心坐标
    centers = []
    for cluster in tower_clusters:
        x_mean = np.mean(cluster[:, 0])
        y_mean = np.mean(cluster[:, 1])
        centers.append((x_mean, y_mean))

    print(f'检测到电力塔个数: {len(centers)}')
    return centers


def cluster_and_fit_powerlines(seg_lines):
    """
    对分段电力线进行聚类和拟合
    
    参数:
        seg_lines: list, 分段电力线点云列表
    
    返回:
        tuple: (endpoints, all_a, all_b, all_p) 所有电力线的端点和拟合参数
    """
    all_a = []
    all_b = []
    all_p = []
    endpoints = []

    for test_line in tqdm(seg_lines, desc='拟合电力线', total=len(seg_lines)):
        # 对每段电力线再次聚类，分离单根导线
        dbscan = DBSCAN(eps=1.5, min_samples=6)
        labels = dbscan.fit_predict(test_line)
        
        max_point_x = np.max(test_line[:, 0])
        min_point_x = np.min(test_line[:, 0])
        
        # 按标签分离单根导线
        clusters = []
        for label in np.unique(labels):
            indices = np.where(labels == label)
            cluster = test_line[indices]
            clusters.append(cluster)
        
        # 对每根导线进行三维拟合
        for cluster in clusters:
            if cluster.shape[0] < 1000:  # 点数太少，跳过
                continue
            
            x = cluster[:, 0]
            y = cluster[:, 1]
            z = cluster[:, 2]
            
            # 第一次拟合
            a, b = ransac_fit_line(x.reshape(-1, 1), y.reshape(-1, 1))
            a_param, x0_param, c_param = fit_catenary_xz(x, z)
            
            # 计算拟合误差
            Y = (a * x + b).flatten()
            Z = a_param * np.cosh((x - x0_param) / a_param) + c_param
            
            delta_y = Y - y
            delta_z = Z - z
            delta = sqrt(delta_y ** 2 + delta_z ** 2)
            rmse = math.sqrt(np.sum(delta ** 2) / delta.size)
            
            # 过滤离群点
            x, y, z, Y, Z = filter_points_by_fit(
                x, y, z, a, b, (a_param, x0_param, c_param), rmse, n=2
            )
            
            # 第二次拟合（使用过滤后的点）
            a, b = ransac_fit_line(x.reshape(-1, 1), y.reshape(-1, 1))
            a_param, x0_param, c_param = fit_catenary_xz(x, z)
            
            # 保存拟合参数
            endpoints.append((min_point_x, max_point_x))
            all_a.append(a)
            all_b.append(b)
            all_p.append([a_param, x0_param, c_param])
    
    return endpoints, all_a, all_b, all_p


def save_fitting_parameters(endpoints, all_a, all_b, all_p, save_path):
    """
    保存拟合参数到文件
    
    参数:
        endpoints: list, 端点坐标列表
        all_a, all_b, all_p: 拟合参数
        save_path: str, 保存路径
    
    返回:
        ndarray: 拼接后的参数矩阵
    """
    all_endpoints_x = np.array(endpoints).reshape(-1, 2)
    all_a = np.array(all_a).reshape(-1, 1)
    all_b = np.array(all_b).reshape(-1, 1)
    all_p = np.array(all_p).reshape(-1, 3)
    
    # 拼接所有参数: [x_start, x_end, a, b, a_cat, x0_cat, c_cat]
    params_matrix = np.concatenate([all_endpoints_x, all_a, all_b, all_p], axis=1)
    
    param_file = Path(save_path) / 'all_argument.txt'
    np.savetxt(param_file, params_matrix)
    
    return all_endpoints_x, all_a, all_b, all_p


def export_to_shapefile(all_endpoints_x, all_a, all_b, all_p, save_path, 
                        output_geographic, utm_zone, is_northern):
    """
    将拟合的电力线导出为Shapefile格式
    
    参数:
        all_endpoints_x: ndarray, 端点x坐标
        all_a, all_b, all_p: 拟合参数
        save_path: str, 保存路径
        output_geographic: bool, 是否输出经纬度坐标
        utm_zone: int, UTM投影带号
        is_northern: bool, 是否为北半球
    """
    # 配置坐标转换器
    if output_geographic:
        # UTM -> WGS84 (经纬度)
        utm_epsg = f"EPSG:326{utm_zone}" if is_northern else f"EPSG:327{utm_zone}"
        transformer = Transformer.from_crs(
            utm_epsg,
            "EPSG:4326",  # WGS84
            always_xy=True
        )
        output_crs = CRS.from_epsg(4326)
    else:
        # 保持UTM坐标
        transformer = None
        utm_epsg_code = 32600 + utm_zone if is_northern else 32700 + utm_zone
        output_crs = CRS.from_epsg(utm_epsg_code)
    
    # 创建Shapefile
    shp_path = Path(save_path) / 'output'
    with shapefile.Writer(shp_path) as w:
        w.field("name", "C")
        
        for i in range(len(all_endpoints_x)):
            # 获取拟合参数
            endpoint = all_endpoints_x[i]
            a = all_a[i]
            b = all_b[i]
            p = all_p[i]
            endpoint_begin, endpoint_end = endpoint
            
            # 生成电力线上的采样点
            x_range = np.linspace(endpoint_begin, endpoint_end, 1000)
            y_range = a * x_range + b
            z_range = p[0] * np.cosh((x_range - p[1]) / p[0]) + p[2]
            
            # 坐标转换
            if output_geographic:
                longitude, latitude = transformer.transform(x_range, y_range)
                point_coords = np.stack([longitude, latitude, z_range], axis=-1).tolist()
            else:
                point_coords = np.stack([x_range, y_range, z_range], axis=-1).tolist()
            
            # 写入Shapefile
            w.linez([point_coords])
            w.record(f"line{i}")
    
    # 创建投影文件
    proj_path = Path(save_path) / "output.prj"
    wkt = output_crs.to_wkt()
    with open(proj_path, 'w') as prj_file:
        prj_file.write(wkt)
    
    coord_type = "经纬度坐标" if output_geographic else f"UTM {utm_zone}{'N' if is_northern else 'S'} 投影坐标"
    print(f'拟合电力线并输出Shapefile文件完成 (坐标系: {coord_type})')
    print(f'输出位置: {save_path}')


# ==================== 主函数 ====================

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='电力线拟合工具 - 从LAS点云提取并拟合电力线',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--line_path",
        type=str,
        default="./line1.las",
        help="电力线点云文件路径"
    )
    
    parser.add_argument(
        "--tower_path",
        type=str,
        default="./tower1.las",
        help="电力塔点云文件路径"
    )
    
    parser.add_argument(
        "--normal_tower_height",
        type=float,
        default=10.0,
        help="正常电力塔高度阈值(米), 用于识别电力塔"
    )
    
    parser.add_argument(
        "--save_path",
        type=str,
        default="./outputs",
        help="结果保存路径"
    )
    
    parser.add_argument(
        "--output_geographic",
        action="store_true",
        help="输出Shapefile使用经纬度坐标(WGS84), 默认使用UTM投影坐标"
    )
    
    parser.add_argument(
        "--utm_zone",
        type=int,
        default=47,
        help="UTM投影带号 (1-60), 默认为47N"
    )
    
    parser.add_argument(
        "--utm_southern",
        action="store_true",
        help="指定UTM投影为南半球, 默认为北半球"
    )
    
    return parser.parse_args()


def main():
    """主函数"""
    # 解析参数
    args = parse_arguments()
    
    # 创建输出目录
    os.makedirs(args.save_path, exist_ok=True)
    
    # 判断UTM半球
    is_northern = not args.utm_southern
    
    print("=" * 60)
    print("电力线拟合工具")
    print("=" * 60)
    print(f"电力线点云: {args.line_path}")
    print(f"电力塔点云: {args.tower_path}")
    print(f"电力塔高度阈值: {args.normal_tower_height}米")
    print(f"输出路径: {args.save_path}")
    print(f"坐标系配置: UTM {args.utm_zone}{'N' if is_northern else 'S'}")
    print(f"输出格式: {'经纬度坐标' if args.output_geographic else 'UTM投影坐标'}")
    print("=" * 60)
    
    # 加载电力线分段模块
    try:
        from pcgrid import SpanSegment
    except ImportError as e:
        print(f"错误: {e}")
        return
    
    # 分段提取电力线
    print("\n正在分段提取电力线...")
    span_segment = SpanSegment()
    seg_lines = span_segment.segment(
        args.line_path,
        args.tower_path,
    )
    print(f"成功提取 {len(seg_lines)} 段电力线")
    
    # 聚类和拟合电力线
    print("\n正在拟合电力线...")
    endpoints, all_a, all_b, all_p = cluster_and_fit_powerlines(seg_lines)
    print(f"成功拟合 {len(endpoints)} 根导线")
    
    # 保存拟合参数
    print("\n正在保存拟合参数...")
    all_endpoints_x, all_a, all_b, all_p = save_fitting_parameters(
        endpoints, all_a, all_b, all_p, args.save_path
    )
    
    # 导出Shapefile
    print("\n正在导出Shapefile...")
    export_to_shapefile(
        all_endpoints_x, all_a, all_b, all_p,
        args.save_path,
        args.output_geographic,
        args.utm_zone,
        is_northern
    )
    
    print("\n" + "=" * 60)
    print("处理完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
