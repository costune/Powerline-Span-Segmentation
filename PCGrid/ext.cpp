
#include <cstdio>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <queue>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "src/base.h"
#include "src/PCGrid.h"
#include "src/Las2PowerLine.h"

namespace py = pybind11;

int EraseRepeatedPoints(std::vector<cLasPOINT>& lasPoints)
{
	if (lasPoints.empty())
		return 0;

	std::sort(lasPoints.begin(), lasPoints.end(), [](const cLasPOINT& a, const cLasPOINT& b) {
		return (a.pt.x() < b.pt.x()) || (a.pt.x() == b.pt.x() && a.pt.y() < b.pt.y()) || (a.pt.x() == b.pt.x() && a.pt.y() == b.pt.y() && a.pt.z() < b.pt.z());
		});

	auto last = std::unique(lasPoints.begin(), lasPoints.end(), [](const cLasPOINT& a, const cLasPOINT& b) {
		return (a.pt.x() == b.pt.x() && a.pt.y() == b.pt.y() && a.pt.z() == b.pt.z());
		});

	lasPoints.erase(last, lasPoints.end());

	return static_cast<int>(lasPoints.size());
}

py::array_t<float> convert_vector_of_eigen_to_numpy(std::vector<Eigen::Vector3f>& vec) {
    
    // 1. 获取向量的总大小
    size_t vec_size = vec.size();

    // 2. 计算 NumPy 数组的尺寸
    // 每个 Eigen::Vector3f 包含 3 个浮点数
    std::vector<py::ssize_t> shape = { static_cast<py::ssize_t>(vec_size), 3 };

    // 3. 计算总元素数
    size_t total_elements = vec_size * 3;
    
    // 4. 创建一个新的 py::array_t，它会分配新的内存
    auto result_array = py::array_t<float>(shape);

    // 5. 访问 NumPy 数组的底层数据
    auto result_buffer = result_array.request();
    float* result_ptr = static_cast<float*>(result_buffer.ptr);

    // 6. 逐个 Eigen::Vector3f 复制数据
    for (size_t i = 0; i < vec_size; ++i) {
        memcpy(result_ptr + i * 3, vec[i].data(), 3 * sizeof(float));
    }
    
    return result_array;
}


py::list process_points(py::array_t<float> lineLasPoints_numpy, py::array_t<float> towerLasPoints_numpy)
{
    // process line points
	std::vector<cLasPOINT> lineLasPoints;
    int nLineLasPoints0;
    {
        if (lineLasPoints_numpy.ndim() != 2 || lineLasPoints_numpy.shape(1) != 3) {
            throw std::runtime_error("Input array must have shape (N, 3)");
        }

        auto n = lineLasPoints_numpy.shape(0);
        nLineLasPoints0 = static_cast<int>(n);
        auto buf = lineLasPoints_numpy.unchecked<2>();  // 只读，无边界检查
        lineLasPoints.reserve(n);

        for (ssize_t i = 0; i < n; i++) {
            cLasPOINT point(Eigen::Vector3f{buf(i, 0), buf(i, 1), buf(i, 2)});
            lineLasPoints.push_back(point);
        }
    }
	int nLineLasPoints1 = EraseRepeatedPoints(lineLasPoints);
	std::cout << "  "
		<< nLineLasPoints1
		<< " valid line points are read, "
		<< std::fixed << std::setprecision(1)
		<< 100.0 * (nLineLasPoints0 - nLineLasPoints1) / nLineLasPoints0
		<< "% points are erased"
		<< std::endl;

     // process tower points
     std::vector<cLasPOINT> towerLasPoints;
	int nTowerLasPoints0;
	{
        if (towerLasPoints_numpy.ndim() != 2 || towerLasPoints_numpy.shape(1) != 3) {
            throw std::runtime_error("Input array must have shape (N, 3)");
        }

        auto n = towerLasPoints_numpy.shape(0);
        nTowerLasPoints0 = static_cast<int>(n);
        auto buf = towerLasPoints_numpy.unchecked<2>();  // 只读，无边界检查
        towerLasPoints.reserve(n);

        for (ssize_t i = 0; i < n; i++) {
            cLasPOINT point(Eigen::Vector3f{buf(i, 0), buf(i, 1), buf(i, 2)});
            towerLasPoints.push_back(point);
        }
     }
	int nTowerLasPoints1 = EraseRepeatedPoints(towerLasPoints);
	std::cout << "  "
		<< nTowerLasPoints1
		<< " valid tower points are read, "
		<< std::fixed << std::setprecision(1)
		<< 100.0 * (nTowerLasPoints0 - nTowerLasPoints1) / nTowerLasPoints0
		<< "% points are erased"
		<< std::endl;

	// Main procedure
    if (nLineLasPoints1 < 100 || nTowerLasPoints1 < 100)
        throw std::runtime_error("Line points or tower points are less than 100");
    
    // 1.1 Tower center extraction using PCGrid 2D
    std::vector<clusterCENTER> towerCenters;
    int nTowers = GetTowerCenter(towerLasPoints, 3.0, towerCenters, false);  // set min height of tower to 3.0

    std::cout << "Found " << (int)towerCenters.size() << " tower centers with 2D grid based DBScan." << std::endl;
    PrintTowerList(towerCenters);

    // 1.2 Mainline powerline extraction using PCGrid 3D
    std::vector<cLasPOINT> mainLineLasPoints;
    {
        std::vector<clusterCENTER> mainLineCenters;
        FilterLinePointsWithPCGrid(lineLasPoints, 3.0, 15, mainLineCenters, mainLineLasPoints);

    }

    // 1.3 The towers are paired and verified using power lines to check if there are sufficient power line points between them.
    CheckTowersWithLinePoints(mainLineLasPoints, 1.0, towerCenters, "");
    std::cout << towerCenters.size() << " towers passed checking" << std::endl;
    PrintTowerList(towerCenters);


    // 2. Power line segmentation between towers
    // 2.1 Init span according to tower centers
    std::vector<lineSEG> lineSegs;
    InitLineSegs(towerCenters, lineSegs);

    std::vector<std::vector<int>> vvSegLinePtIndices;

    // 2.2 根据点杆塔连线的垂直距离，对电力线进行分档
    int nGroups = GroupPowerLinePoints(mainLineLasPoints, towerCenters, lineSegs, vvSegLinePtIndices);

    py::list outputs;
    {
        for (auto i = 0; i < vvSegLinePtIndices.size(); ++i)
        {
            std::vector<Eigen::Vector3f> pts;
            pts.reserve(vvSegLinePtIndices[i].size());

            for (int idx : vvSegLinePtIndices[i])
            {
                if (idx >= 0 && idx < static_cast<int>(mainLineLasPoints.size()))
                    pts.push_back(mainLineLasPoints[idx].pt);
            }

            if (pts.empty()) continue;

            outputs.append(convert_vector_of_eigen_to_numpy(pts));
        }
    }
    return outputs;
}




PYBIND11_MODULE(_C, m) {
     m.def("process_points", &process_points);
}