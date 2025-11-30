#include <cstdio>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <queue>
#include <cmath>
#include <chrono>  // Include time header
#include <string>
#include <iomanip>

#include "base.h"
#include "PCGrid.h"
#include "Las2PowerLine.h"




void PrintTowerList(const std::vector<clusterCENTER>& towerCenters)
{
	// Print the tower centers
	for (size_t i = 0; i < towerCenters.size(); ++i)
	{
		if (0 == i)
			std::cout << "\tTower Center " << std::setw(5) << i
				<< ": X=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters[i].pt.x()
				<< ", Y=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters[i].pt.y()
				<< ", maxZ=" << std::setw(8) << std::fixed << std::setprecision(2) << towerCenters[i].maxZ
				<< ", radius=" << std::setw(5) << std::fixed << std::setprecision(1) << towerCenters[i].get2DRadius()
				<< ", hei=" << std::setw(5) << std::fixed << std::setprecision(1) << towerCenters[i].getHeight()
				<< ", nPts=" << std::setw(6) << towerCenters[i].numOfPts
				<< std::endl;

		else {
			float dX = towerCenters[i].pt.x() - towerCenters[i - 1].pt.x();
			float dY = towerCenters[i].pt.y() - towerCenters[i - 1].pt.y();
			float dZ = towerCenters[i].pt.z() - towerCenters[i - 1].pt.z();
			float dist = sqrt((double)dX * dX + (double)dY * dY);
			std::cout << "\tTower Center " << std::setw(5) << i
				<< ": X=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters[i].pt.x()
				<< ", Y=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters[i].pt.y()
				<< ", maxZ=" << std::setw(8) << std::fixed << std::setprecision(2) << towerCenters[i].maxZ
				<< ", radius=" << std::setw(5) << std::fixed << std::setprecision(1) << towerCenters[i].get2DRadius()
				<< ", hei=" << std::setw(5) << std::fixed << std::setprecision(1) << towerCenters[i].getHeight()
				<< ", nPts=" << std::setw(6) << towerCenters[i].numOfPts
				<< ", dX=" << std::setw(9) << std::fixed << std::setprecision(1) << dX
				<< ", dY=" << std::setw(9) << std::fixed << std::setprecision(1) << dY
				<< ", dZ=" << std::setw(7) << std::fixed << std::setprecision(2) << dZ
				<< ", dist=" << std::setw(10) << std::fixed << std::setprecision(1) << dist
				<< std::endl;

		}
	}
}

int GetTowerCenter(std::vector<cLasPOINT>& towerPoints, float minHei, std::vector<clusterCENTER>& towerCenters, const bool debug)
{
	std::cout << "Detecting Tower using DBScan from " << towerPoints.size() << " points..." << std::endl;

	{
		auto startTime = std::chrono::high_resolution_clock::now();

		std::vector<clusterCENTER> towerCenters0;
		std::vector<std::vector<int> > vClusterPointIndices;
		int nTowers = PCGrid_2D(towerPoints, 3.0, towerCenters0, vClusterPointIndices);  // use 3 as grid size
		
		{
			auto PCGrid_2D_Time = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(PCGrid_2D_Time - startTime);
			std::cout << "  "
				<< "PCGrid_2D Time taken: "
				<< std::fixed << std::setprecision(2)
				<< duration.count()
				<< " seconds"
				<< std::endl;
		}		

		towerCenters.reserve(nTowers);

		std::cout << "  " << towerCenters0.size() << " towers center detected by DBScan, abnormal Towers:" << std::endl;

		// check the z distribution of tower points
		int histZ[400];		// suppose the height of tower is less than 400m
		for (int i = nTowers - 1; i >= 0; --i)
		{
			if (towerCenters0[i].getHeight() < minHei)
			{
				if ( debug ) 
					std::cout << "\tTower Center " << std::setw(5) << i
						<< ": X=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters0[i].pt.x()
						<< ", Y=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters0[i].pt.y()
						<< ", Z=" << std::setw(8) << std::fixed << std::setprecision(2) << towerCenters0[i].pt.z()
						<< ", hei=" << std::setw(5) << std::fixed << std::setprecision(1) << towerCenters0[i].getHeight()
						<< ", nPts=" << std::setw(6) << towerCenters0[i].numOfPts
						<< std::endl;
			}
			else
			{
				float minZ = towerCenters0[i].minZ;
				float maxZ = towerCenters0[i].maxZ;
				for (int j = 0; j < vClusterPointIndices[i].size(); ++j)
				{
					int index = vClusterPointIndices[i][j];
					if (minZ > towerPoints[index].pt.z())
						minZ = towerPoints[index].pt.z();
					if (maxZ < towerPoints[index].pt.z())
						maxZ = towerPoints[index].pt.z();
				}

				memset(histZ, 0, sizeof(int) * 400);
				for(int j = 0; j < vClusterPointIndices[i].size(); ++j)
				{
					int index = vClusterPointIndices[i][j];
					float z = towerPoints[index].pt.z() - minZ; 
					if (z >= 0 && z < 400)
						histZ[(int)z]++;
					else {
						std::cout << "\tWarning: z-zmin=" << z << " out of range [0, 400)" << std::endl;
					}
				}
				// count the number of bins with zero count
				int nBins = maxZ - minZ + 0.5;
				int nZeroBins = 0;
				for (int j = 0; j < nBins; ++j)
				{
					if (0 == histZ[j])
						++nZeroBins;
				}

				if (3 * nZeroBins > nBins) {
					std::cout << "\tTower Center " << std::setw(5) << i
						<< ": X=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters0[i].pt.x()
						<< ", Y=" << std::setw(10) << std::fixed << std::setprecision(2) << towerCenters0[i].pt.y()
						<< ", Z=" << std::setw(8) << std::fixed << std::setprecision(2) << towerCenters0[i].pt.z()
						<< ", hei=" << std::setw(5) << std::fixed << std::setprecision(1) << towerCenters0[i].getHeight()
						<< ", nPts=" << std::setw(6) << towerCenters0[i].numOfPts
						<< ", Z distribution abnormal" << std::endl;

				}
				else 
					towerCenters.push_back(towerCenters0[i]);
			}
		}

		//delete histZ;
	}

	int nTowers = (int)towerCenters.size();

	std::cout << "  " << nTowers << " towers after filtering" << std::endl;

	return nTowers;	
}



// 对塔中心点进行排序，按距离第一个点的距离从小到大排序
void SortTowerCentersWithDistance(std::vector<clusterCENTER>& towerCenters, Eigen::Vector3f curCenter)
{
	// 按距离排序，起点为点云的第一个点
	// Eigen::Vector3f curCenter = towerPoints[0].pt;

	float minDist = calculateDistanceSquare_2D(curCenter, towerCenters[0].pt);
	int jMin = 0;

	int nTowers = (int)towerCenters.size();

	std::vector<float> vDists;
	vDists.resize(nTowers);
	for( int i=0; i<nTowers-1; ++i)
	{
		for (int j = i+1; j < nTowers; j++)
		{
			float dist = calculateDistanceSquare_2D(curCenter, towerCenters[j].pt);
			if ( minDist > dist )
			{
				minDist = dist;
				jMin = j;
			}
		}

		// 与前一个点的距离
		vDists[i] = minDist;

		// 如果最小距离塔与当前基准塔不同，则交换位置
		if (jMin > i) {
			// 交换位置
			clusterCENTER temp = towerCenters[i];
			towerCenters[i] = towerCenters[jMin];
			towerCenters[jMin] = temp;
		}
		
		// 更新最小距离的塔为当前基准塔
		curCenter = towerCenters[i].pt;

		// 更新jMin为下一个默认的最近距离塔
		jMin = i + 1; 
		minDist = calculateDistanceSquare_2D(curCenter, towerCenters[jMin].pt);
	}

	std::vector<clusterCENTER> sortedTowerCenters;

	int iNeedResort = -1;
	for (int i = 1; i < nTowers; i++)
	{
		float dist = calculateDistanceSquare_2D(towerCenters[i].pt, towerCenters[0].pt);

		// 找到第一个反转点
		// 如果与前一个塔的距离大于与第一个塔的距离，则说明该点应该链接到第一个塔之前
		if (dist < vDists[i])
		{
			sortedTowerCenters.resize(nTowers);

			// 将i前面的点反转
			for (int j = 0; j < i; ++j)
			{
				sortedTowerCenters[j] = towerCenters[i-1 - j];
			}

			// 将剩余的点复制到后面
			for (int j = i; j < nTowers; ++j)
			{
				sortedTowerCenters[j] = towerCenters[j];
			}

			iNeedResort = i;
			break;
		}
	}

	// 对i点之后的点，重新进行排序
	if ( iNeedResort > 0 )	{
		towerCenters = std::move(sortedTowerCenters);

		// 设定当前基准塔为iNeedResort-1的塔
		curCenter = towerCenters[iNeedResort-1].pt;

		// 更新jMin为下一个默认的最近距离塔
		jMin = iNeedResort;
		minDist = calculateDistanceSquare_2D(curCenter, towerCenters[jMin].pt);

		for (int i = iNeedResort; i < nTowers - 1; ++i)
		{
			for (int j = i + 1; j < nTowers; j++)
			{
				float dist = calculateDistanceSquare_2D(curCenter, towerCenters[j].pt);
				if (minDist > dist)
				{
					minDist = dist;
					jMin = j;
				}
			}

			// 如果最小距离塔与当前基准塔不同，则交换位置
			if (jMin > i) {
				// 交换位置
				clusterCENTER temp = towerCenters[i];
				towerCenters[i] = towerCenters[jMin];
				towerCenters[jMin] = temp;
			}

			// 更新最小距离的塔为当前基准塔
			curCenter = towerCenters[i].pt;

			// 更新jMin为下一个默认的最近距离塔
			jMin = i + 1;
			minDist = calculateDistanceSquare_2D(curCenter, towerCenters[jMin].pt);
		}
	}

}






void InitLineSegs(const std::vector<clusterCENTER>& towerCenters, std::vector<lineSEG>& lineSegs)
{
	int nTowers = towerCenters.size();
	lineSegs.resize(nTowers);

	for (int i = 0; i < nTowers - 1; i++)
		lineSegs[i].InitLineSeg(towerCenters[i].pt, towerCenters[i + 1].pt);

	int i = nTowers - 1;

	// 最后一个塔之后没有连线
	lineSegs[i].X0 = towerCenters[i].pt.x();
	lineSegs[i].Y0 = towerCenters[i].pt.y();
	lineSegs[i].cosA = 1;	lineSegs[i].sinA = 0;	
	lineSegs[i].len = 0;
}


// 电力线分档
int GroupPowerLinePoints(const std::vector<cLasPOINT>& lineLasPoints,
	const std::vector<clusterCENTER>& towerCenters,		// 杆塔位置
	const std::vector<lineSEG> &lineSegs,				// 杆塔连线
	std::vector<std::vector<int> >& groupedLinePts)
{
	using PointCloud_T = PointCloudAdaptor<clusterCENTER>;
	// 建一个KD树来存储塔中心点，加快查找速度
	PointCloud_T pointCloud_T(towerCenters);

	typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<float, PointCloud_T>,
		PointCloud_T,
		3 /* dim */
	> KDTree;

	KDTree index(3, pointCloud_T, { 10 /* max leaf */ });
	index.buildIndex();

	int nTowers = (int)towerCenters.size();
	groupedLinePts.resize(nTowers);

	int nLinePts = lineLasPoints.size();
	for (int i = 0; i < nLinePts; ++i)
	{
		auto& pt = lineLasPoints[i].pt;
		std::vector<int> neighborPts;

		uint32_t ret_matches[1];
		float  ret_dists[1];

		index.knnSearch(&lineLasPoints[i].pt.x(), 1, ret_matches, ret_dists);

		int indexTower = ret_matches[0];

		float offsetX, dist;
		if ( lineSegs[indexTower].IsPointInSeg(pt, &offsetX, &dist)) 
		{

			float r0 = towerCenters[indexTower].get2DRadius();
			float r1 = towerCenters[indexTower+1].get2DRadius();

			float r = r0 + (r1 - r0) * offsetX / lineSegs[indexTower].len + 1;

			if (fabs(dist) < r)
				groupedLinePts[indexTower].push_back(i);
		}
		else if (indexTower > 0 && lineSegs[indexTower - 1].IsPointInSeg(pt, &offsetX, &dist))
		{
			float r0 = towerCenters[indexTower-1].get2DRadius();
			float r1 = towerCenters[indexTower].get2DRadius();

			float r = r0 + (r1 - r0) * offsetX / lineSegs[indexTower-1].len + 1;

			if (fabs(dist) < r)
				groupedLinePts[indexTower-1].push_back(i);
		}
		else if (indexTower < nTowers - 2 && lineSegs[indexTower + 1].IsPointInSeg(pt, &offsetX, &dist))
		{
			float r0 = towerCenters[indexTower + 1].get2DRadius();
			float r1 = towerCenters[indexTower + 2].get2DRadius();

			float r = r0 + (r1 - r0) * offsetX / lineSegs[indexTower + 1].len + 1;

			if (fabs(dist) < r)
				groupedLinePts[indexTower+1].push_back(i);
		}
	}

	return 0;
}

