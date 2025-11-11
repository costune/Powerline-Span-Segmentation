#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
电力线隐患检测系统
Power Line Hazard Detection System

功能：检测点云中距离电力线过近的地面点（树障检测）
"""

import os
import argparse
from pathlib import Path
from typing import List, Union, Dict, Tuple

import numpy as np
from numpy.typing import ArrayLike
import laspy
import shapefile
from tqdm import tqdm
from pyproj import CRS, Transformer
from shapely.geometry import Point, Polygon
from multiprocessing import Pool
import psutil


# ============================================================================
# 常量定义
# ============================================================================
DEFAULT_UTM_ZONE = 47  # 默认UTM分区47N
DEFAULT_HEMISPHERE = 'N'  # 默认北半球
DEFAULT_THRESHOLD = 8.0  # 默认距离阈值(米)
DEFAULT_LINE_PARAM_PATH = "./outputs/all_argument.txt"
DEFAULT_OUTPUT_FOLDER = "./alarm"
CHUNK_SIZE = 10000  # 并行计算块大小
MAX_WORKERS = 8  # 最大并行进程数


# ============================================================================
# 工具类
# ============================================================================
class ProjTransformer:
    """坐标系转换工具类"""
    
    def __init__(self, crs_from: CRS, crs_to: CRS):
        """
        初始化坐标转换器
        
        Args:
            crs_from: 源坐标系
            crs_to: 目标坐标系
        """
        self.crs_from = crs_from
        self.crs_to = crs_to
        self.transformer = Transformer.from_crs(
            crs_from, crs_to, always_xy=True, accuracy=1e-10
        )

    def transform(self, points: np.ndarray) -> np.ndarray:
        """
        批量转换二维点坐标
        
        Args:
            points: 二维点集 [n, 2]
            
        Returns:
            转换后的点集 [n, 2]
        """
        result_points = []
        for point in points:
            x, y = self.transformer.transform(point[0], point[1])
            result_points.append([x, y])
        return np.array(result_points)


# ============================================================================
# 并行计算函数
# ============================================================================
def calculate_distances(chunk: np.ndarray, line_points: np.ndarray, 
                       radius: float) -> np.ndarray:
    """
    计算点云到电力线的距离（并行化）
    
    Args:
        chunk: 点云数据块
        line_points: 电力线采样点
        radius: 距离阈值
        
    Returns:
        是否在阈值内的布尔掩码
    """
    distances = np.linalg.norm(
        chunk[:, None, :] - line_points[None, :, :], axis=-1
    )
    return np.any(distances <= radius, axis=1)


def calculate_distances_with_values(chunk: np.ndarray, line_points: np.ndarray,
                                    radius: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算点云到电力线的距离并返回最小距离值（并行化）
    
    Args:
        chunk: 点云数据块
        line_points: 电力线采样点
        radius: 距离阈值
        
    Returns:
        (是否在阈值内的布尔掩码, 最小距离值)
    """
    distances = np.linalg.norm(
        chunk[:, None, :] - line_points[None, :, :], axis=-1
    )
    within_threshold = np.any(distances <= radius, axis=1)
    min_distances = np.min(distances, axis=1)
    return within_threshold, min_distances


def parallel_distance_computation(lidar_xyz: np.ndarray, line_points: np.ndarray,
                                  radius: float, n_jobs: int = 4,
                                  chunk_size: int = CHUNK_SIZE,
                                  need_distances: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    并行计算点云与电力线的距离关系
    
    Args:
        lidar_xyz: 点云坐标
        line_points: 电力线采样点
        radius: 距离阈值
        n_jobs: 并行进程数
        chunk_size: 数据块大小
        need_distances: 是否返回距离值
        
    Returns:
        布尔掩码或(布尔掩码, 距离值)
    """
    # 分块处理
    chunks = [
        lidar_xyz[i:i + chunk_size] 
        for i in range(0, lidar_xyz.shape[0], chunk_size)
    ]
    inputs = [(chunk, line_points, radius) for chunk in chunks]
    
    if need_distances:
        with Pool(n_jobs) as pool:
            results = pool.starmap(calculate_distances_with_values, inputs)
        
        within_threshold = [result[0] for result in results]
        distances = [result[1] for result in results]
        
        idx_within = np.unique(np.concatenate(within_threshold).nonzero()[0])
        dist_values = np.concatenate(distances)[idx_within]
        
        mask = np.zeros(lidar_xyz.shape[0], dtype=bool)
        mask[idx_within] = True
        
        return mask, dist_values
    else:
        with Pool(n_jobs) as pool:
            results = pool.starmap(calculate_distances, inputs)
        
        idx_within = np.unique(np.concatenate(results).nonzero()[0])
        mask = np.zeros(lidar_xyz.shape[0], dtype=bool)
        mask[idx_within] = True
        
        return mask


# ============================================================================
# 核心功能函数
# ============================================================================
def generate_line_sample_points(line_param: np.ndarray, 
                               sigma: float) -> np.ndarray:
    """
    生成电力线采样点
    
    Args:
        line_param: 电力线参数 [n_lines, 7] (x_begin, x_end, a, b, p0, p1, p2)
        sigma: 采样间距
        
    Returns:
        电力线采样点 [n_lines, n_samples, 3]
    """
    end_points = line_param[:, :2]
    a = line_param[:, 2:3]
    b = line_param[:, 3:4]
    p = line_param[:, 4:]
    
    # 计算采样点数
    min_num = int(
        np.ceil(np.abs(end_points[:, 0] - end_points[:, 1]) / sigma).max() + 1
    )
    
    # 生成采样点
    x = np.linspace(end_points[:, 0], end_points[:, 1], min_num).T
    y = x * a + b
    z = p[:, 0:1] * np.cosh((x - p[:, 1:2]) / p[:, 0:1]) + p[:, 2:3]
    
    line_points = np.stack([x, y, z], axis=-1)
    return line_points


def remove_power_line_points(lidar_xyz: np.ndarray, line_param: np.ndarray,
                             radius: float = 1.0) -> np.ndarray:
    """
    从点云中移除电力线点
    
    Args:
        lidar_xyz: 原始点云坐标
        line_param: 电力线参数
        radius: 移除半径阈值
        
    Returns:
        移除电力线后的点云
    """
    line_points = generate_line_sample_points(line_param, radius)
    
    # 获取系统CPU数量
    cpu_count = psutil.cpu_count(logical=False)
    n_jobs = min(cpu_count, MAX_WORKERS)
    
    for one_line in tqdm(line_points, desc="移除电力线点云", total=len(line_points)):
        x_begin, x_end = one_line[0, 0], one_line[-1, 0]
        
        # 筛选在当前电力线X范围内的点
        lidar_mask = np.logical_and(
            lidar_xyz[:, 0] >= x_begin,
            lidar_xyz[:, 0] <= x_end
        )
        mask_indices = np.where(lidar_mask)[0]
        lidar_one_line = lidar_xyz[lidar_mask]
        
        # 并行计算需要移除的点
        mask_to_remove = parallel_distance_computation(
            lidar_one_line, one_line, radius, n_jobs=n_jobs, chunk_size=CHUNK_SIZE
        )
        
        # 更新点云
        true_indices_to_remove = mask_indices[mask_to_remove]
        keep_mask = np.ones(lidar_xyz.shape[0], dtype=bool)
        keep_mask[true_indices_to_remove] = False
        lidar_xyz = lidar_xyz[keep_mask]
    
    return lidar_xyz


def detect_hazard_points(lidar_xyz: np.ndarray, line_param: np.ndarray,
                        radius: float = 20.0, 
                        sigma: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    检测隐患点
    
    Args:
        lidar_xyz: 地面点云坐标
        line_param: 电力线参数
        radius: 检测半径阈值
        sigma: 电力线采样间距
        
    Returns:
        (隐患点坐标, 距离值, 所属电力线编号)
    """
    line_points = generate_line_sample_points(line_param, sigma)
    
    cpu_count = psutil.cpu_count(logical=False)
    n_jobs = min(cpu_count, MAX_WORKERS)
    
    alarm_masks = []
    distance_arrays = []
    line_num_mapping = {}
    
    for line_num, one_line in tqdm(enumerate(line_points), 
                                   desc="检测隐患点", 
                                   total=len(line_points)):
        x_begin, x_end = one_line[0, 0], one_line[-1, 0]
        
        # 筛选在当前电力线X范围内的点
        lidar_mask = np.logical_and(
            lidar_xyz[:, 0] >= x_begin,
            lidar_xyz[:, 0] <= x_end
        )
        mask_indices = np.where(lidar_mask)[0]
        lidar_one_line = lidar_xyz[lidar_mask]
        
        # 并行计算距离
        mask_to_alarm, distances = parallel_distance_computation(
            lidar_one_line, one_line, radius, 
            n_jobs=n_jobs, chunk_size=CHUNK_SIZE, need_distances=True
        )
        
        # 构建全局掩码和距离数组
        true_indices_to_alarm = mask_indices[mask_to_alarm]
        global_mask = np.zeros(lidar_xyz.shape[0], dtype=bool)
        global_mask[true_indices_to_alarm] = True
        
        global_distances = np.full(lidar_xyz.shape[0], np.inf)
        global_distances[true_indices_to_alarm] = distances
        
        if np.any(global_mask):
            alarm_masks.append(global_mask)
            distance_arrays.append(global_distances)
            line_num_mapping[len(alarm_masks) - 1] = line_num
    
    # 合并所有电力线的结果
    if not alarm_masks:
        return np.array([]), np.array([]), np.array([])
    
    alarm_mask = np.any(alarm_masks, axis=0)
    alarm_points = lidar_xyz[alarm_mask]
    
    # 计算每个点的最小距离和对应电力线
    distance_matrix = np.array(distance_arrays).T
    min_distances = np.min(distance_matrix, axis=1)
    min_line_indices = np.argmin(distance_matrix, axis=1)
    
    valid_mask = np.isfinite(min_distances)
    min_distances = min_distances[valid_mask]
    min_line_indices = min_line_indices[valid_mask]
    
    line_numbers = np.fromiter(
        (line_num_mapping[idx] for idx in min_line_indices),
        dtype=int
    )
    
    return alarm_points, min_distances, line_numbers


def get_distance_level(distance: float) -> int:
    """
    根据距离计算危险等级
    
    Args:
        distance: 距离值
        
    Returns:
        危险等级 (0-3)
    """
    if distance <= 6.0:
        return 0
    elif distance <= 8.0:
        return 1
    elif distance <= 10.0:
        return 2
    else:
        return 3


def output_alarm_shapefile(alarm_points: np.ndarray, distances: np.ndarray,
                          line_numbers: np.ndarray, output_file: str,
                          crs_output: CRS, crs_source: CRS = None):
    """
    输出隐患点Shapefile
    
    Args:
        alarm_points: 隐患点坐标
        distances: 距离值
        line_numbers: 电力线编号
        output_file: 输出文件路径
        crs_output: 输出坐标系
        crs_source: 源坐标系（如需转换）
    """
    # 坐标转换
    if crs_source and crs_source != crs_output:
        transformer = Transformer.from_crs(
            crs_source.srs, crs_output.srs, always_xy=True
        )
        x, y = transformer.transform(alarm_points[:, 0], alarm_points[:, 1])
        output_points = np.stack([x, y, alarm_points[:, 2]], axis=-1)
    else:
        output_points = alarm_points
    
    # 输出Shapefile
    if output_file.endswith('.shp') or output_file.endswith('.dbf'):
        output_file = output_file[:-4]
    
    output_path = Path(output_file)
    os.makedirs(output_path.parent, exist_ok=True)
    
    with shapefile.Writer(str(output_path), shapeType=shapefile.POINTZ) as shp:
        shp.field("ID", "N")
        shp.field("Distance", "N", decimal=2)
        shp.field("DistLevel", "N")
        shp.field("LineNum", "C")
        
        for idx, (x, y, z) in enumerate(output_points):
            shp.pointz(x, y, z)
            shp.record(
                ID=idx,
                Distance=distances[idx],
                DistLevel=get_distance_level(distances[idx]),
                LineNum=f'line{line_numbers[idx]}'
            )
    
    # 输出投影文件
    prj_path = str(output_path) + '.prj'
    with open(prj_path, "w") as prj:
        prj.write(crs_output.to_wkt())
    
    print(f"隐患点已成功保存到 {output_file}.shp")


# ============================================================================
# 主处理函数
# ============================================================================
def detect_hazards(lidar_path: str, line_parameters: str, threshold: float,
                  output_folder: str, output_crs: str = 'utm',
                  utm_zone: int = DEFAULT_UTM_ZONE,
                  utm_hemisphere: str = DEFAULT_HEMISPHERE):
    """
    主检测函数
    
    Args:
        lidar_path: 点云文件路径
        line_parameters: 电力线参数文件路径
        threshold: 距离阈值
        output_folder: 输出文件夹
        output_crs: 输出坐标系 ('utm' 或 'latlon')
        utm_zone: UTM分区号
        utm_hemisphere: UTM半球 ('N' 或 'S')
    """
    print("=" * 60)
    print("电力线隐患检测系统")
    print("=" * 60)
    
    # 读取点云数据
    print("\n[1/5] 读取点云文件...")
    lidar = laspy.read(lidar_path)
    
    # 确定点云坐标系
    try:
        lidar_proj = lidar.header.vlrs[0].string
        crs_lidar = CRS.from_string(lidar_proj)
        print(f"点云坐标系: {crs_lidar.name}")
    except Exception:
        epsg_code = 32600 + utm_zone if utm_hemisphere == 'N' else 32700 + utm_zone
        crs_lidar = CRS.from_epsg(epsg_code)
        print(f"[警告] 文件中不存在投影信息，使用默认: WGS84 / UTM {utm_zone}{utm_hemisphere}")
    
    # 读取电力线参数
    print("\n[2/5] 读取电力线参数...")
    line_param = np.loadtxt(line_parameters)
    print(f"电力线数量: {len(line_param)}")
    
    # 移除电力线点云
    print("\n[3/5] 移除电力线点云...")
    removed_dir = Path(lidar_path).parent / 'removed'
    removed_path = removed_dir / (Path(lidar_path).stem + '.npy')
    
    if removed_path.exists():
        print(f"加载已处理文件: {removed_path}")
        lidar_xyz = np.load(removed_path)
    else:
        removed_dir.mkdir(parents=True, exist_ok=True)
        # 每隔100个点采样一次
        lidar_xyz = remove_power_line_points(
            np.array(lidar.xyz)[::100], line_param
        )
        np.save(removed_path, lidar_xyz)
        print(f"处理结果已保存: {removed_path}")
    
    print(f"处理后点云数量: {len(lidar_xyz)}")
    
    # 检测隐患点
    print(f"\n[4/5] 检测隐患点 (阈值: {threshold}m)...")
    alarm_points, distances, line_numbers = detect_hazard_points(
        lidar_xyz, line_param, radius=threshold
    )
    print(f"检测到隐患点数量: {len(alarm_points)}")
    
    if len(alarm_points) == 0:
        print("未检测到隐患点")
        return
    
    # 确定输出坐标系
    if output_crs.lower() == 'latlon':
        crs_output = CRS.from_epsg(4326)  # WGS84经纬度
        print("\n输出坐标系: WGS84 经纬度")
        crs_source = crs_lidar
    else:  # utm
        epsg_code = 32600 + utm_zone if utm_hemisphere == 'N' else 32700 + utm_zone
        crs_output = CRS.from_epsg(epsg_code)
        print(f"\n输出坐标系: WGS84 / UTM {utm_zone}{utm_hemisphere}")
        crs_source = crs_lidar if crs_lidar != crs_output else None
    
    # 输出结果
    print(f"\n[5/5] 输出隐患点Shapefile...")
    output_file = Path(output_folder) / f"alarm_{int(threshold)}m"
    output_alarm_shapefile(
        alarm_points, distances, line_numbers,
        str(output_file), crs_output, crs_source
    )
    
    print("\n" + "=" * 60)
    print("树障检测完成")
    print("=" * 60)


# ============================================================================
# 命令行接口
# ============================================================================
def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='电力线隐患检测系统 - 树障检测',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 使用默认参数（输出UTM 47N坐标系）
  python cal_distance.py --lidar_path data/point_cloud.las
  
  # 输出经纬度坐标系
  python cal_distance.py --lidar_path data/point_cloud.las --output_crs latlon
  
  # 指定UTM分区
  python cal_distance.py --lidar_path data/point_cloud.las --utm_zone 48
        """
    )
    
    parser.add_argument(
        "--lidar_path",
        type=str,
        required=True,
        help="点云文件路径 (.las/.laz格式)"
    )
    
    parser.add_argument(
        "--line_parameters", "-l",
        type=str,
        default=DEFAULT_LINE_PARAM_PATH,
        help=f"电力线参数文件路径 (默认: {DEFAULT_LINE_PARAM_PATH})"
    )
    
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"距离报警阈值(米) (默认: {DEFAULT_THRESHOLD})"
    )
    
    parser.add_argument(
        "--output_folder", "-o",
        type=str,
        default=DEFAULT_OUTPUT_FOLDER,
        help=f"输出文件夹路径 (默认: {DEFAULT_OUTPUT_FOLDER})"
    )
    
    parser.add_argument(
        "--output_crs",
        type=str,
        choices=['utm', 'latlon'],
        default='utm',
        help="输出坐标系类型: utm=UTM投影坐标, latlon=WGS84经纬度 (默认: utm)"
    )
    
    parser.add_argument(
        "--utm_zone",
        type=int,
        default=DEFAULT_UTM_ZONE,
        help=f"UTM分区号 (1-60) (默认: {DEFAULT_UTM_ZONE})"
    )
    
    parser.add_argument(
        "--utm_hemisphere",
        type=str,
        choices=['N', 'S'],
        default=DEFAULT_HEMISPHERE,
        help=f"UTM半球: N=北半球, S=南半球 (默认: {DEFAULT_HEMISPHERE})"
    )
    
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()
    
    # 参数验证
    if not os.path.exists(args.lidar_path):
        print(f"错误: 点云文件不存在: {args.lidar_path}")
        return 1
    
    if not os.path.exists(args.line_parameters):
        print(f"错误: 电力线参数文件不存在: {args.line_parameters}")
        return 1
    
    if not (1 <= args.utm_zone <= 60):
        print(f"错误: UTM分区号必须在1-60之间")
        return 1
    
    # 执行检测
    try:
        detect_hazards(
            lidar_path=args.lidar_path,
            line_parameters=args.line_parameters,
            threshold=args.threshold,
            output_folder=args.output_folder,
            output_crs=args.output_crs,
            utm_zone=args.utm_zone,
            utm_hemisphere=args.utm_hemisphere
        )
        return 0
    except Exception as e:
        print(f"\n错误: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
