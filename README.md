## fit_line: 电力线提取及树障检测

### 环境设置

```bash
$ conda env create -f enviroment.yml
```

### 命令行运行

**1. 电力线提取**

```bash
$ python fit_line.py --line_path path/to/line --tower_path path/to/las 
```

必需参数：`--line_path` 电力线文件位置，`--tower_path` 电力塔文件位置；可选参数：
`--normal_tower_height` 线路电力塔的最低高度，默认为10；`--save_path` 成果保存位置，默认为 `./outputs` 文件夹。

注意：算法仅能够处理在UTM坐标系东西向坐标值 x 不发生折回的点云，即对于每个 x 坐标，仅对应于一个电力走廊切面。

**2. 树障检测**

```bash
$ python cal_distance.py --lidar_path path/to/original/lidar --line_parameters outputs/all_argument.txt
```

必需参数：`--lidar_path` 原始点云文件位置，`--line_parameters` 上一部拟合得到的电力线参数文件位置，在上一步得到的成果文件夹中；可选参数：`--threshold` 检测距离阈值，默认为8；`--output_folder` 成果保存位置，默认为 `./alarm` 文件夹。



boost 1.65.0