#ifndef _DENSE_CLUSTER_H_
#define _DENSE_CLUSTER_H_

#include <iostream>
#include <vector>
#include <map>

#include "nanoflann/nanoflann.hpp"
#include "base.h"





// Function to calculate the Euclidean distance between two points.
inline float calculateDistanceSquare_3D(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2)
{
	return std::pow(p1.x() - p2.x(), 2) + std::pow(p1.y() - p2.y(), 2) + std::pow(p1.z() - p2.z(), 2);
}

inline float calculateDistanceSquare_2D(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2)
{
	return std::pow(p1.x() - p2.x(), 2) + std::pow(p1.y() - p2.y(), 2);
}



struct GridCell
{
public:
	float Xw, Yw, Zw;		// weight center of the grid cell
	int count;				// Number of points in this grid cell
	int lastPtIndex;		// Index of the last point added to this cell, used for clustering

	float minZ, maxZ;		// Minimum and maximum Z values in this cell

	int clusterId;			// 所属的聚类ID, -1表示未分配。 在电力线和杆塔一致性验证总，此ID代表杆塔连线的ID
	int nClusterSamples;	// 所属的聚类ID的采样数。 在电力线和杆塔一致性验证总，此处相当于最近杆塔连线的采样数

public:
	GridCell() : Xw(0),Yw(0),Zw(0), count(0), lastPtIndex(-1),
		minZ(9999.0), maxZ(-9999), 
		clusterId(-1), nClusterSamples(0) {} // Default constructor

	void RegisterPoint(int ptIndex, cLasPOINT &cPt)
	{
		if (0 == count) {
			minZ = cPt.pt.z();
			maxZ = cPt.pt.z();

			lastPtIndex = ptIndex;
			count = 1;
		}
		else {
			cPt.nextPt = lastPtIndex;	// Link to the last point in this cell
			lastPtIndex = ptIndex;		// Update to the last point index

			if ( minZ > cPt.pt.z())
				minZ = cPt.pt.z();
			if ( maxZ < cPt.pt.z())
				maxZ = cPt.pt.z();

			++count;
		}
	}

	void Stat(const std::vector<cLasPOINT>& cPoints)
	{
		int ptIndex = lastPtIndex;

		// 统计每个格网单元的重心
		double XwSum = 0, YwSum = 0, ZwSum = 0;
		double Xw2Sum = 0, Yw2Sum = 0, Zw2Sum = 0;

		for( int i=0; i<count; ++i)
		{
			XwSum += cPoints[ptIndex].pt.x();
			YwSum += cPoints[ptIndex].pt.y();
			ZwSum += cPoints[ptIndex].pt.z();

			Xw2Sum += (double)cPoints[ptIndex].pt.x()* cPoints[ptIndex].pt.x();
			Yw2Sum += (double)cPoints[ptIndex].pt.y()* cPoints[ptIndex].pt.y();
			Zw2Sum += (double)cPoints[ptIndex].pt.z()* cPoints[ptIndex].pt.z();
		}

		Xw = XwSum / count;
		Yw = YwSum / count; 
		Zw = ZwSum / count;
	}
};


class PointGrid2D
{
private:
	float m_X0, m_Y0;		// Origin of the grid
	int64_t m_nX, m_nY;		// Number of grid cells in each dimension

	float m_gridSize;

	std::map<int64_t, GridCell> m_cells; // Map of grid cells

public:
	PointGrid2D() : m_X0(0), m_Y0(0), m_nX(0), m_nY(0), m_gridSize(0)
	{
	}

	PointGrid2D(float x0, float y0, int nX, int nY, float gridSize)
		: m_X0(x0), m_Y0(y0), m_nX(nX), m_nY(nY), m_gridSize(gridSize)
	{
		m_cells.clear();
	}

	void Init(float x0, float y0, int nX, int nY, float gridSize)
	{
		m_X0 = x0;
		m_Y0 = y0;
		m_nX = nX;
		m_nY = nY;
		m_cells.clear();
	}

	inline int64_t Ixy2Id(int ix, int iy) const
	{
		return ix + iy * m_nX;
	}

	inline void Id2IxIy(int64_t cellId, int& ix, int& iy) const
	{
		iy = cellId / m_nX;
		ix = cellId - (int64_t)iy * m_nX;
	}

	inline GridCell* getCell(int64_t cellId)
	{
		auto it = m_cells.find(cellId);

		if (it != m_cells.end())
			return &it->second;

		return NULL;
	}

	inline GridCell* getCell(int ix, int iy)
	{
		return getCell(ix + iy * m_nX);
	}

	inline int GetCount(int64_t cellId) const
	{
		auto it = m_cells.find(cellId);

		if (it != m_cells.end())
			return it->second.count;
		else
			return 0; // No points in this cell
	}

	inline int GetCount(int ix, int iy) const
	{
		return GetCount(ix + iy * m_nX);
	}

	void RegisterPoints( std::vector<cLasPOINT>& cPoints)
	{
		int nPoints = cPoints.size();
		for (int i = 0; i < nPoints; ++i)
		{
			int ix = (int)((cPoints[i].pt.x() - m_X0) / m_gridSize);
			int iy = (int)((cPoints[i].pt.y() - m_Y0) / m_gridSize);

			if (ix < 0 || ix > m_nX - 1 || iy < 0 || iy > m_nY - 1 )
				continue;

			int64_t cellId = ix + iy * m_nX;

			auto cell = m_cells.find(cellId);
			if (m_cells.find(cellId) == m_cells.end()) {
				GridCell cell0;

				cell0.RegisterPoint(i, cPoints[i]);
				m_cells[cellId] = cell0;
			}
			else {
				cell->second.RegisterPoint(i, cPoints[i]);
			}
		}

		// 统计每个格网单元的重心
		for( auto& cell : m_cells)
		{
			cell.second.Stat(cPoints);
		}
	}

	void MarkPointsCluster(std::vector<cLasPOINT>& cPoints)
	{
		int nPoints = cPoints.size();
		for (int i = 0; i < nPoints; ++i)
		{
			int ix = (int)((cPoints[i].pt.x() - m_X0) / m_gridSize);
			int iy = (int)((cPoints[i].pt.y() - m_Y0) / m_gridSize);

			int64_t cellId = ix + iy * m_nX;

			auto cell = m_cells.find(cellId);
			if (cell != m_cells.end())
			{
				cPoints[i].clusterId = cell->second.clusterId;
			}
		}
	}

	// get neighbors of a cell and the count of points in those neighbors
	int GetNeighborPointCount(int ix, int iy) const
	{
		int count = 0;
		for (int dx = -1; dx <= 1; ++dx)
		{
			if (ix + dx < 0 || ix + dx >= m_nX - 1)
				continue;
			for (int dy = -1; dy <= 1; ++dy)
			{
				if (0 == dx && 0 == dy)
					continue;
				if ( iy + dy < 0 || iy + dy >= m_nY-1 )
					continue;

				int64_t cellId = ix + dx + (iy + dy) * m_nX;

				auto nb = m_cells.find(cellId);

				if (nb != m_cells.end()) 
					count += nb->second.count;
			}
		}
		return count;
	}


	int GetNeighbors(int ix, int iy, std::vector<int64_t>& neighbors) const
	{
		int count = 0;
		for (int dx = -1; dx <= 1; ++dx)
		{
			if (ix + dx < 0 || ix + dx >= m_nX - 1)
				continue;

			for (int dy = -1; dy <= 1; ++dy)
			{
				if (0 == dx && 0 == dy)
					continue;

				if ( iy + dy < 0 || iy + dy >= m_nY - 1)
					continue;

				int64_t cellId = ix + dx + (iy + dy) * m_nX;

				auto nb = m_cells.find(cellId);

				if (nb != m_cells.end()) {
					count += nb->second.count;
					neighbors.push_back(cellId);
				}
			}
		}

		return count;
	};

	int GetNeighbors(int64_t cellId, std::vector<int64_t>& neighbors) const
	{
		int ix, iy;
		Id2IxIy(cellId, ix, iy);
		return GetNeighbors(ix, iy, neighbors);
	}

	void GetSortedCells(std::multimap<int, int64_t>& sortedCells) const
	{
		sortedCells.clear();
		for (const auto& cell : m_cells)
		{
			// Insert the point count and index into the sorted map
			sortedCells.insert(std::make_pair(cell.second.count, cell.first));
		}
	}

	int GetAverageCount(float* std, int *minCount, int *maxCount) const
	{
		if (m_cells.empty())
			return 0;

		*minCount = std::numeric_limits<int>::max();
		*maxCount = std::numeric_limits<int>::min();
		double nSum = 0;
		double nnSum = 0;
		for (const auto& cell : m_cells)
		{
			*minCount = std::min(*minCount, cell.second.count);
			*maxCount = std::max(*maxCount, cell.second.count);

			nSum += cell.second.count;
			nnSum += (double)cell.second.count * cell.second.count;
		}

		double mean = nSum / m_cells.size();
		*std = std::sqrt(nnSum / m_cells.size() - mean * mean);

		return (int)mean;
	}

	// get geometric centers of the grid cells
	int GetCellCenters(std::vector<clusterCENTER>& vCellCenter) const
	{
		vCellCenter.resize(m_cells.size());

		int i = 0;
		for (const auto& cell : m_cells)
		{
			int64_t cellId = cell.first;
			int ix, iy;

			Id2IxIy(cellId, ix, iy);

			vCellCenter[i].pt.x() = m_X0 + (ix + 0.5f) * m_gridSize;
			vCellCenter[i].pt.y() = m_Y0 + (iy + 0.5f) * m_gridSize;
			vCellCenter[i].pt.z() = 0;

			vCellCenter[i].numOfPts = cell.second.count;

			// Initialize radius to zero

			i++;
		}

		return vCellCenter.size();
	}

	// get the average centers of points in each grid cell
	int GetWeightedCellCenters(std::vector<clusterCENTER>& vCellCenter) const
	{
		vCellCenter.resize(m_cells.size());

		int i = 0;
		for (const auto& cell : m_cells)
		{
			vCellCenter[i].pt.x() = cell.second.Xw;
			vCellCenter[i].pt.y() = cell.second.Yw;
			vCellCenter[i].pt.z() = 0.0f;

			vCellCenter[i].numOfPts = cell.second.count;

			// Initialize radius to zero

			i++;
		}

		return vCellCenter.size();
	}

	// get the average centers of points in each grid cell
	int GetCellIds(std::vector<int> &vCellIds) const
	{
		vCellIds.resize(m_cells.size());

		int i = 0;
		for (const auto& cell : m_cells)
		{
			vCellIds[i] = cell.first;
			i++;
		}

		return vCellIds.size();
	}


	// 根据邻域格网内的点数进行噪声滤除
	void FilerNoise()
	{
		std::multimap<int, int64_t> sortedCells;
		GetSortedCells(sortedCells);

		// Calculate the average point count in the grid cells	
		int iHalf = sortedCells.size()/2;
		auto cellHalf = sortedCells.begin();
		for(int i=0; i<iHalf; ++i)
			++cellHalf;

		int countTh = cellHalf->first;

		for ( auto iter = m_cells.begin(); iter != m_cells.end(); iter++)
		{
			int64_t cellId = iter->first;
			int ix, iy;
			Id2IxIy(cellId, ix, iy);

			int count = iter->second.count;
			count += GetNeighborPointCount(ix, iy);	
			if (count < countTh)
				iter = m_cells.erase( iter );
		}
	}
};


class PointGrid3D
{
private:
	float m_X0, m_Y0, m_Z0;			// Origin of the grid
	int64_t m_nX, m_nY, m_nZ;
	int64_t m_nXY;					// Number of grid cells in each dimension

	float m_gridSize;

	std::map<int64_t, GridCell> m_cells; // Map of grid cells

public:
	PointGrid3D() : m_X0(0), m_Y0(0), m_Z0(0), m_nX(0), m_nY(0), m_nZ(0), m_gridSize(0)
	{
		m_nXY = 0;
	}

	PointGrid3D(float x0, float y0, float z0, int nX, int nY, int nZ, float gridSize)
		: m_X0(x0), m_Y0(y0), m_Z0(z0), m_nX(nX), m_nY(nY), m_nZ(nZ), m_gridSize(gridSize)
	{
		m_nXY = m_nX * m_nY; // Precompute the number of cells in XY plane
		m_cells.clear();
	}

	void Init(float x0, float y0, float z0, int nX, int nY, int nZ, float gridSize)
	{
		m_X0 = x0;
		m_Y0 = y0;
		m_Z0 = z0;
		m_nX = nX;
		m_nY = nY;
		m_nZ = nZ;
		m_nXY = m_nX * m_nY; // Precompute the number of cells in XY plane

		m_cells.clear();
	}

	inline int64_t Ixyz2Index(int ix, int iy, int iz) const
	{
		return ix + iy * m_nX + iz * m_nXY;
	}

	inline void Id2IxIyIz(int64_t cellId, int& ix, int& iy, int& iz) const
	{
		iz = cellId / m_nXY;

		int64_t cellId2D = cellId - iz * m_nXY;

		iy = cellId2D / m_nX;
		ix = cellId2D - iy * m_nX;
	}

	inline GridCell* getCell(int64_t cellId)
	{
		auto it = m_cells.find(cellId);

		if (it != m_cells.end())
			return &it->second;

		return NULL;
	}

	inline GridCell* getCell(int ix, int iy, int iz)
	{
		return getCell(ix + iy * m_nX + iz * m_nXY);
	}

	inline int GetCount(int64_t cellId)  const
	{
		auto it = m_cells.find(cellId);
		if (it != m_cells.end())
			return it->second.count;
		else
			return 0; // No points in this cell
	}

	inline int GetCount(int ix, int iy, int iz)  const
	{
		return GetCount(ix + iy * m_nX + iz * m_nXY);
	}

	int GetAverageCount(float *std) const
	{
		if (m_cells.empty())
			return 0;

		double nSum = 0;
		double nnSum = 0;
		for (const auto& cell : m_cells)
		{
			nSum += cell.second.count;
			nnSum += (double)cell.second.count * cell.second.count;
		}

		double mean = nSum / m_cells.size();
		*std = std::sqrt( nnSum/m_cells.size() - mean * mean );

		return (int)mean;
	}

	void RegisterPoints(std::vector<cLasPOINT>& cPoints)
	{
		GridCell cell;
		cell.count = 1;

		int nPoints = cPoints.size();
		for (int i = 0; i < nPoints; ++i)
		{
			int ix = (int)((cPoints[i].pt.x() - m_X0) / m_gridSize);
			int iy = (int)((cPoints[i].pt.y() - m_Y0) / m_gridSize);
			int iz = (int)((cPoints[i].pt.z() - m_Z0) / m_gridSize);

			if( ix < 0 || ix > m_nX-1 || iy < 0 || iy > m_nY - 1 || iz < 0 || iz > m_nZ - 1)
				continue;

			int64_t cellId = ix + iy * m_nX + iz * m_nXY;
			auto cell = m_cells.find(cellId);

			if (m_cells.find(cellId) == m_cells.end()) {
				GridCell cell0;

				cell0.RegisterPoint(i, cPoints[i]);
				m_cells[cellId] = cell0;
			}
			else {
				cell->second.RegisterPoint(i, cPoints[i]);
			}
		}

		// 统计每个格网单元的重心
		for (auto& cell : m_cells)
		{
			cell.second.Stat(cPoints);
		}
	}

	void MarkPointsCluster(std::vector<cLasPOINT>& cPoints)
	{
		int nPoints = cPoints.size();
		for (int i = 0; i < nPoints; ++i)
		{
			int ix = (int)((cPoints[i].pt.x() - m_X0) / m_gridSize);
			int iy = (int)((cPoints[i].pt.y() - m_Y0) / m_gridSize);
			int iz = (int)((cPoints[i].pt.z() - m_Z0) / m_gridSize);

			int64_t cellId = ix + iy * m_nX + iz * m_nXY;

			auto cell = m_cells.find(cellId);
			if (cell != m_cells.end())
			{
				cPoints[i].clusterId = cell->second.clusterId;
			}
		}
	}


	// get neighbors of a cell and the count of points in those neighbors
	int GetNeighborPointCount(int ix, int iy, int iz)  const
	{
		int count = 0;
		for (int dx = -1; dx <= 1; ++dx)
		{
			if (ix + dx < 0 || ix + dx >= m_nX-1)
				continue;

			for (int dy = -1; dy <= 1; ++dy)
			{
				if (iy + dy < 0 || iy + dy >= m_nY - 1)
					continue;

				for (int dz = -1; dz <= 1; ++dz)
				{
					if (0 == dx && 0 == dy && 0 == dz)
						continue;

					if ( iz+dz <0 || iz+dz >= m_nZ-1 )
						continue;

					int64_t cellId = ix + dx + (iy + dy) * m_nX + (iz + dz) * m_nXY;

					auto nb = m_cells.find(cellId);

					if (nb != m_cells.end()) 
						count += nb->second.count;
				}
			}
		}

		return count;
	};

	// get neighbors of a cell and the count of points in those neighbors
	int GetNeighbors(int ix, int iy, int iz, std::vector<int64_t>& neighbors)  const
	{
		int count = 0;
		for (int dx = -1; dx <= 1; ++dx)
		{
			if (ix + dx < 0 || ix + dx >= m_nX - 1)
				continue;
			for (int dy = -1; dy <= 1; ++dy)
			{
				if (iy + dy < 0 || iy + dy >= m_nY - 1)
					continue;
				for (int dz = -1; dz <= 1; ++dz)
				{
					if (0 == dx && 0 == dy && 0 == dz)
						continue;

					 if( iz + dz < 0 || iz + dz >= m_nZ-1)
						continue;

					int64_t cellId = ix + dx + (iy + dy) * m_nX + (iz + dz) * m_nXY;

					auto nb = m_cells.find(cellId);

					if (nb != m_cells.end()) {
						count += nb->second.count;
						neighbors.push_back(cellId);
					}
				}
			}
		}

		return count;
	};

	int GetNeighbors(int64_t cellId, std::vector<int64_t>& neighbors) const
	{
		int ix, iy, iz;
		Id2IxIyIz(cellId, ix, iy, iz);

		return GetNeighbors(ix, iy, iz, neighbors);
	}

	// 
	void GetSortedCells(std::multimap<int, int>& sortedCells) const
	{
		sortedCells.clear();
		for (const auto& cell : m_cells)
		{
			// Insert the point count and index into the sorted map
			sortedCells.insert(std::make_pair(cell.second.count, cell.first));
		}
	}

	int GetCellCenters(std::vector<clusterCENTER> &vCellCenter) const
	{
		vCellCenter.resize(m_cells.size());
		
		int i = 0;
		for (const auto& cell : m_cells)
		{
			int64_t cellId = cell.first;
			int ix, iy, iz;

			Id2IxIyIz(cellId, ix, iy, iz);

			vCellCenter[i].pt.x() = m_X0 + (ix + 0.5f) * m_gridSize;
			vCellCenter[i].pt.y() = m_Y0 + (iy + 0.5f) * m_gridSize;
			vCellCenter[i].pt.z() = m_Z0 + (iz + 0.5f) * m_gridSize;

			vCellCenter[i].numOfPts = cell.second.count;

			// Initialize radius to zero

			i++;
		}

		return vCellCenter.size();
	}

	// get the average centers of points in each grid cell
	int GetWeightedCellCenters(std::vector<clusterCENTER>& vCellCenter) const
	{
		vCellCenter.resize(m_cells.size());

		int i = 0;
		for (const auto& cell : m_cells)
		{
			vCellCenter[i].pt.x() = cell.second.Xw;
			vCellCenter[i].pt.y() = cell.second.Yw;
			vCellCenter[i].pt.z() = cell.second.Zw;

			vCellCenter[i].numOfPts = cell.second.count;

			i++;
		}

		return vCellCenter.size();
	}

	// get all cells
	int GetCellIds(std::vector<int64_t>& vCellIds) const
	{
		vCellIds.resize(m_cells.size());

		int i = 0;
		for (const auto& cell : m_cells)
		{
			vCellIds[i] = cell.first;
			i++;
		}

		return vCellIds.size();
	}

	// get all cells
	int GetCellCounts(std::vector<int>& vCellCounts) const
	{
		vCellCounts.resize(m_cells.size());

		int i = 0;
		for (const auto& cell : m_cells)
		{
			vCellCounts[i] = cell.second.count;
			i++;
		}

		return vCellCounts.size();
	}

	// 根据邻域格网内的点数进行噪声滤除
	void FilerNoise()
	{
		std::multimap<int, int> sortedCells;
		GetSortedCells(sortedCells);

		// Calculate the average point count in the grid cells	
		int iHalf = sortedCells.size() / 2;
		auto cellHalf = sortedCells.begin();
		for (int i = 0; i < iHalf; ++i)
			++cellHalf;

		int countTh = cellHalf->first;

		for (auto iter = m_cells.begin(); iter != m_cells.end(); iter++)
		{
			int64_t cellId = iter->first;
			int ix, iy,iz;
			Id2IxIyIz(cellId, ix, iy,iz);

			int count = iter->second.count;
			count += GetNeighborPointCount(ix, iy,iz);
			if (count < countTh)
				iter = m_cells.erase(iter);
		}
	}
};




int PCGrid_2D(std::vector<cLasPOINT>& cPoints, double eps, std::vector<clusterCENTER>& dbScanCenters, std::vector<std::vector<int> >& vvClusterPoints);
int PCGrid_3D(std::vector<cLasPOINT>& cPoints, double eps, std::vector<clusterCENTER>& dbScanCenters, std::vector<std::vector<int> >& vvClusterPoints);

void FilterLinePointsWithPCGrid(std::vector<cLasPOINT>& lineLasPoints, double eps, int minPts, std::vector<clusterCENTER>& clusterCenters, std::vector<cLasPOINT>& mainLinePoints);
void CheckTowersWithLinePoints(std::vector<cLasPOINT>& lineLasPoints, float gridSize, std::vector<clusterCENTER>& towerCenters, const std::string& outDir);


#endif // _DENSE_CLUSTER_H_