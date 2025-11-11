# cal_distance.py 使用说明

## 功能概述
电力线隐患检测系统 - 从点云中检测距离电力线过近的地面点（树障检测）

## 主要改进

### 1. 新增参数
- `--output_crs`: 输出坐标系类型
  - `utm`: UTM投影坐标（默认）
  - `latlon`: WGS84经纬度
- `--utm_zone`: UTM分区号 (1-60)，默认47
- `--utm_hemisphere`: UTM半球 (N/S)，默认N

### 2. 代码优化
- 清晰的模块化结构
- 完整的函数文档
- 改进的错误处理
- 移除冗余调试代码
- 统一的命名规范

## 使用示例

### 基本用法（使用默认UTM 47N）
```bash
python cal_distance.py --lidar_path data/lijiang_line1/110kV-Merge.las
```

### 输出经纬度坐标
```bash
python cal_distance.py --lidar_path data/lijiang_line1/110kV-Merge.las --output_crs latlon
```

### 指定不同UTM分区
```bash
python cal_distance.py --lidar_path data/lijiang_line1/110kV-Merge.las --utm_zone 48
```

### 完整参数示例
```bash
python cal_distance.py \
    --lidar_path data/lijiang_line1/110kV-Merge.las \
    --line_parameters ./test_result_lijiang1/all_argument.txt \
    --threshold 10 \
    --output_folder ./alarm_results \
    --output_crs utm \
    --utm_zone 47 \
    --utm_hemisphere N
```

### 查看帮助
```bash
python cal_distance.py --help
```

## 参数说明

| 参数 | 简写 | 类型 | 默认值 | 说明 |
|------|------|------|--------|------|
| `--lidar_path` | - | str | 必填 | 点云文件路径 (.las/.laz) |
| `--line_parameters` | `-l` | str | ./outputs/all_argument.txt | 电力线参数文件路径 |
| `--threshold` | `-t` | float | 8.0 | 距离报警阈值(米) |
| `--output_folder` | `-o` | str | ./alarm | 输出文件夹路径 |
| `--output_crs` | - | str | utm | 输出坐标系(utm/latlon) |
| `--utm_zone` | - | int | 47 | UTM分区号(1-60) |
| `--utm_hemisphere` | - | str | N | UTM半球(N/S) |

## 输出说明

### 输出文件
- `alarm_{阈值}m.shp`: 隐患点Shapefile（点要素）
- `alarm_{阈值}m.dbf`: 属性表
- `alarm_{阈值}m.shx`: 索引文件
- `alarm_{阈值}m.prj`: 投影信息

### 属性字段
- `ID`: 隐患点编号
- `Distance`: 到电力线的最小距离(米)
- `DistLevel`: 危险等级 (0-3)
  - 0: ≤6m (高危)
  - 1: 6-8m (较高)
  - 2: 8-10m (中等)
  - 3: >10m (较低)
- `LineNum`: 所属电力线编号

## 坐标系说明

### UTM坐标系 (默认)
- 适合大范围测量和计算
- 单位为米，方便距离计算
- 默认使用 WGS84 / UTM 47N
- 可通过 `--utm_zone` 指定分区

### 经纬度坐标系
- WGS84地理坐标系
- 适合在地图上直接显示
- 使用 `--output_crs latlon` 指定

## 注意事项

1. 如果点云文件包含投影信息，将自动使用文件中的坐标系
2. 如果文件无投影信息，将使用指定的UTM分区
3. 处理后的点云会缓存在 `removed/` 文件夹，加快后续处理
4. 建议根据实际情况调整 `--threshold` 参数
5. 输出的shapefile坐标系由 `--output_crs` 参数控制

## 代码结构

```
cal_distance.py
├── 常量定义
├── 工具类 (ProjTransformer)
├── 并行计算函数
│   ├── calculate_distances
│   ├── calculate_distances_with_values
│   └── parallel_distance_computation
├── 核心功能函数
│   ├── generate_line_sample_points
│   ├── remove_power_line_points
│   ├── detect_hazard_points
│   ├── get_distance_level
│   └── output_alarm_shapefile
├── 主处理函数 (detect_hazards)
└── 命令行接口 (main)
