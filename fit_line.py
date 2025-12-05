#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Power line fitting script
Extracts and fits power lines from LAS point cloud data and outputs results in Shapefile format
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


# ==================== Fitting Functions ====================

def fit_catenary_xz(x, z):
    """
    Fit the x-z plane projection using a catenary model.

    Catenary model: z(x) = a*cosh((x-x0)/a) + c

    Parameters:
        x: ndarray, x coordinates
        z: ndarray, z coordinates

    Returns:
        tuple: (a, x0, c) parameters of the catenary model
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
    Fit a line using linear regression.

    Parameters:
        x: ndarray, shape (n, 1), x coordinates
        y: ndarray, shape (n, 1), y coordinates

    Returns:
        tuple: (a, b) line parameters, y = a*x + b
    """
    model = LinearRegression()
    model.fit(x, y)
    
    a = model.coef_
    b = model.intercept_
    
    return a, b


def filter_points_by_fit(x, y, z, a, b, catenary_params, threshold, n=1):
    """
    Filter outliers based on the fitted models.

    Parameters:
        x, y, z: coordinate arrays
        a, b: line fit parameters in xy plane
        catenary_params: catenary parameters for xz plane (a, x0, c)
        threshold: filtering threshold
        n: threshold multiplier

    Returns:
        tuple: filtered (x, y, z, Y, Z)
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


# ==================== Point Cloud Processing Functions ====================

def get_tower_centers(tower_path, normal_tower_height):
    """
    Extract tower center coordinates from tower point cloud.

    Parameters:
        tower_path: str, path to the tower point cloud file
        normal_tower_height: float, normal tower height threshold

    Returns:
        list: tower center coordinates [(x1, y1), (x2, y2), ...]
    """
    las = laspy.read(tower_path)
    las_x = np.array(las.x)
    las_y = np.array(las.y)
    las_z = np.array(las.z)
    data = np.stack([las_x, las_y, las_z], axis=1)

    # DBSCAN clustering
    dbscan = DBSCAN(eps=3, min_samples=15)
    labels = dbscan.fit_predict(data)

    # Filter tower clusters
    tower_clusters = []
    for label in np.unique(labels):
        if label == -1:  # skip noise
            continue
            
        indices = np.where(labels == label)
        cluster = data[indices]
        
        # compute cluster dimensions
        height = np.max(cluster[:, 2]) - np.min(cluster[:, 2])
        width = np.max(cluster[:, 0]) - np.min(cluster[:, 0])
        length = np.max(cluster[:, 1]) - np.min(cluster[:, 1])

        # identify towers by height characteristics
        # example rough ranges: 110kV: 15-30m, 220kV: 30-40m, 500kV: 40-60m, 1000kV: 60-80+m
        if height > normal_tower_height and height > width and height > length:
            tower_clusters.append(cluster)

    # compute the center coordinate of each tower cluster
    centers = []
    for cluster in tower_clusters:
        x_mean = np.mean(cluster[:, 0])
        y_mean = np.mean(cluster[:, 1])
        centers.append((x_mean, y_mean))

    print(f'Detected towers: {len(centers)}')
    return centers


def cluster_and_fit_powerlines(seg_lines):
    """
    Cluster segmented lines and fit models to each conductor.

    Parameters:
        seg_lines: list of segmented line point clouds

    Returns:
        tuple: (endpoints, all_a, all_b, all_p) endpoints and fitting parameters
    """
    all_a = []
    all_b = []
    all_p = []
    endpoints = []

    for test_line in tqdm(seg_lines, desc='Fitting power lines', total=len(seg_lines)):
        # cluster each segment again to separate individual conductors
        dbscan = DBSCAN(eps=1.5, min_samples=6)
        labels = dbscan.fit_predict(test_line)
        
        max_point_x = np.max(test_line[:, 0])
        min_point_x = np.min(test_line[:, 0])
        
        # separate individual conductors by label
        clusters = []
        for label in np.unique(labels):
            indices = np.where(labels == label)
            cluster = test_line[indices]
            clusters.append(cluster)
        
        # perform 3D fitting for each conductor
        for cluster in clusters:
            if cluster.shape[0] < 1000:  # too few points, skip
                continue
            
            x = cluster[:, 0]
            y = cluster[:, 1]
            z = cluster[:, 2]
            
            # first fit
            a, b = ransac_fit_line(x.reshape(-1, 1), y.reshape(-1, 1))
            a_param, x0_param, c_param = fit_catenary_xz(x, z)
            
            # compute fitting error
            Y = (a * x + b).flatten()
            Z = a_param * np.cosh((x - x0_param) / a_param) + c_param
            
            delta_y = Y - y
            delta_z = Z - z
            delta = sqrt(delta_y ** 2 + delta_z ** 2)
            rmse = math.sqrt(np.sum(delta ** 2) / delta.size)
            
            # filter outliers
            x, y, z, Y, Z = filter_points_by_fit(
                x, y, z, a, b, (a_param, x0_param, c_param), rmse, n=2
            )
            
            # second fit (using filtered points)
            a, b = ransac_fit_line(x.reshape(-1, 1), y.reshape(-1, 1))
            a_param, x0_param, c_param = fit_catenary_xz(x, z)
            
            # store fitted parameters
            endpoints.append((min_point_x, max_point_x))
            all_a.append(a)
            all_b.append(b)
            all_p.append([a_param, x0_param, c_param])
    
    return endpoints, all_a, all_b, all_p


def save_fitting_parameters(endpoints, all_a, all_b, all_p, save_path):
    """
    Save fitting parameters to disk.

    Parameters:
        endpoints: list of endpoint coordinate pairs
        all_a, all_b, all_p: fitted parameters
        save_path: str, directory to save results

    Returns:
        ndarray: concatenated parameter matrix
    """
    all_endpoints_x = np.array(endpoints).reshape(-1, 2)
    all_a = np.array(all_a).reshape(-1, 1)
    all_b = np.array(all_b).reshape(-1, 1)
    all_p = np.array(all_p).reshape(-1, 3)
    
    # concatenate all parameters: [x_start, x_end, a, b, a_cat, x0_cat, c_cat]
    params_matrix = np.concatenate([all_endpoints_x, all_a, all_b, all_p], axis=1)
    
    param_file = Path(save_path) / 'all_argument.txt'
    np.savetxt(param_file, params_matrix)
    
    return all_endpoints_x, all_a, all_b, all_p


def export_to_shapefile(all_endpoints_x, all_a, all_b, all_p, save_path, 
                        output_geographic, utm_zone, is_northern):
    """
    Export the fitted power lines as a Shapefile.

    Parameters:
        all_endpoints_x: ndarray, endpoints x coordinates
        all_a, all_b, all_p: fitting parameters
        save_path: str, directory to save results
        output_geographic: bool, whether to output geographic coordinates (WGS84)
        utm_zone: int, UTM zone number
        is_northern: bool, True if northern hemisphere
    """
    # 配置坐标转换器
    if output_geographic:
        # UTM -> WGS84
        utm_epsg = f"EPSG:326{utm_zone}" if is_northern else f"EPSG:327{utm_zone}"
        transformer = Transformer.from_crs(
            utm_epsg,
            "EPSG:4326",  # WGS84
            always_xy=True
        )
        output_crs = CRS.from_epsg(4326)
    else:
        # keep UTM coordinates
        transformer = None
        utm_epsg_code = 32600 + utm_zone if is_northern else 32700 + utm_zone
        output_crs = CRS.from_epsg(utm_epsg_code)
    
    # create Shapefile
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
            
            # generate sampled points along the power line
            x_range = np.linspace(endpoint_begin, endpoint_end, 1000)
            y_range = a * x_range + b
            z_range = p[0] * np.cosh((x_range - p[1]) / p[0]) + p[2]
            
            # coordinate transformation
            if output_geographic:
                longitude, latitude = transformer.transform(x_range, y_range)
                point_coords = np.stack([longitude, latitude, z_range], axis=-1).tolist()
            else:
                point_coords = np.stack([x_range, y_range, z_range], axis=-1).tolist()
            
            # 写入Shapefile
            w.linez([point_coords])
            w.record(f"line{i}")
    
    # create projection file (.prj)
    proj_path = Path(save_path) / "output.prj"
    wkt = output_crs.to_wkt()
    with open(proj_path, 'w') as prj_file:
        prj_file.write(wkt)

    coord_type = "Geographic coordinates (WGS84)" if output_geographic else f"UTM {utm_zone}{'N' if is_northern else 'S'}"
    print(f'Power line fitting and Shapefile export completed (CRS: {coord_type})')
    print(f'Output path: {save_path}')


# ==================== 主函数 ====================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Power Line Fitting Tool - extract and fit power lines from LAS point clouds',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--line_path",
        type=str,
        default="./line1.las",
        help="Path to power line point cloud file"
    )
    
    parser.add_argument(
        "--tower_path",
        type=str,
        default="./tower1.las",
        help="Path to tower point cloud file"
    )
    
    parser.add_argument(
        "--normal_tower_height",
        type=float,
        default=10.0,
        help="Normal tower height threshold (meters), used to identify towers"
    )
    
    parser.add_argument(
        "--save_path",
        type=str,
        default="./outputs",
        help="Directory to save results"
    )
    
    parser.add_argument(
        "--output_geographic",
        action="store_true",
        help="Output Shapefile in geographic coordinates (WGS84). Default is UTM projection coordinates"
    )
    
    parser.add_argument(
        "--utm_zone",
        type=int,
        default=47,
        help="UTM zone number (1-60), default is 47N"
    )
    
    parser.add_argument(
        "--utm_southern",
        action="store_true",
        help="Specify UTM projection as southern hemisphere (default is northern hemisphere)"
    )
    
    return parser.parse_args()


def main():
    """Main entry point"""
    # 解析参数
    args = parse_arguments()
    
    # 创建输出目录
    os.makedirs(args.save_path, exist_ok=True)
    
    # 判断UTM半球
    is_northern = not args.utm_southern
    
    print("=" * 60)
    print("Power Line Fitting Tool")
    print("=" * 60)
    print(f"Power line point cloud: {args.line_path}")
    print(f"Tower point cloud: {args.tower_path}")
    print(f"Tower height threshold: {args.normal_tower_height} m")
    print(f"Output directory: {args.save_path}")
    print(f"CRS configuration: UTM {args.utm_zone}{'N' if is_northern else 'S'}")
    print(f"Output format: {'Geographic coordinates (WGS84)' if args.output_geographic else 'UTM projection coordinates'}")
    print("=" * 60)
    
    # 加载电力线分段模块
    try:
        from pcgrid import SpanSegment
    except ImportError as e:
        print(f"错误: {e}")
        return
    
    # segment the power lines
    print("\nSegmenting power lines...")
    span_segment = SpanSegment()
    seg_lines = span_segment.segment(
        args.line_path,
        args.tower_path,
    )
    print(f"Successfully extracted {len(seg_lines)} power line segments")
    
    # cluster and fit power lines
    print("\nFitting power lines...")
    endpoints, all_a, all_b, all_p = cluster_and_fit_powerlines(seg_lines)
    print(f"Successfully fitted {len(endpoints)} conductors")
    
    # save fitting parameters
    print("\nSaving fitting parameters...")
    all_endpoints_x, all_a, all_b, all_p = save_fitting_parameters(
        endpoints, all_a, all_b, all_p, args.save_path
    )
    
    # export to Shapefile
    print("\nExporting Shapefile...")
    export_to_shapefile(
        all_endpoints_x, all_a, all_b, all_p,
        args.save_path,
        args.output_geographic,
        args.utm_zone,
        is_northern
    )
    
    print("\n" + "=" * 60)
    print("Processing completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
