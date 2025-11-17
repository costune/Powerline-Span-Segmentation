#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "LasFileIo.h"
#include "DenseCluster.h"
#include "PointCloutFeature.h"

extern bool g_bSaveFile;

void PrintTowerList(const std::vector<clusterCENTER>& towerCenters);

void SortTowerCentersWithDistance(std::vector<clusterCENTER>& towerCenters, Eigen::Vector3f curCenter);

// Function to calculate the point count in each grid cell
// Point count in grids based clustering algorithm
int PCGrid_3D(std::vector<cLasPOINT> &cPoints, double eps, std::vector<clusterCENTER>& dbScanCenters, std::vector<std::vector<int> >& vvClusterPoints)
{
	int nPoints = cPoints.size();
	if (nPoints <= 0 || cPoints.empty()) {
		fprintf(stderr, "No points to process.\n");
		return 0;
	}

	// Determine the grid size
	double gridSize = eps; // sqrt(3); // or eps / 2; depends on how strict you want to be

	// Find the bounding box of the data
	double minX = std::numeric_limits<double>::max();
	double minY = std::numeric_limits<double>::max();
	double minZ = std::numeric_limits<double>::max();
	double maxX = std::numeric_limits<double>::lowest();
	double maxY = std::numeric_limits<double>::lowest();
	double maxZ = std::numeric_limits<double>::lowest();

	for (const auto& point : cPoints) {
		minX = std::min(minX, (double)point.pt.x());
		minY = std::min(minY, (double)point.pt.y());
		minZ = std::min(minZ, (double)point.pt.z());
		maxX = std::max(maxX, (double)point.pt.x());
		maxY = std::max(maxY, (double)point.pt.y());
		maxZ = std::max(maxZ, (double)point.pt.z());
	}

	// Calculate the number of grid cells in each dimension
	int nX = (int)std::ceil((maxX - minX) / gridSize);
	int nY = (int)std::ceil((maxY - minY) / gridSize);
	int nZ = (int)std::ceil((maxZ - minZ) / gridSize);

	PointGrid3D grid(minX, minY, minZ, nX, nY, nZ, gridSize);
	
	grid.RegisterPoints(cPoints);

	// too slow to use std::multimap, use a vector instead
	// std::multimap <int, int> sortedCells;
	// grid.GetSortedCells(sortedCells);

	// minimum point count threshold
	int minPts;

	if( false ) {
		float stdGridPts;
		int meanGridPts = grid.GetAverageCount(&stdGridPts);

		if (meanGridPts <= 0) {
			fprintf(stderr, "No points to process for tower centers.\n");
			return 0;
		}

		minPts = meanGridPts - (int)stdGridPts;
		if (minPts < 15) minPts = 15;
	}
	else
	{
		std::vector<int> vCellCounts;
		grid.GetCellCounts( vCellCounts );
		std::sort(vCellCounts.begin(), vCellCounts.end());

		int nCells = vCellCounts.size();
		minPts = vCellCounts[nCells / 9]; 

		std::cout << "  PCGrid_3D has " << nCells << " sorted cells, minPts = " << minPts << std::endl;
	}

	if (minPts < 15) minPts = 15;

	// clustering
	std::vector<int64_t> vCellIds;
	grid.GetCellIds(vCellIds);

	int clusterId = 0;
	for (auto &cellId: vCellIds)
	{
		GridCell *cell = grid.getCell( cellId );

		// has been visited;
		if (-1 != cell->clusterId)
			continue;

		std::vector<int64_t> neighbors;
		int neighborCount = grid.GetNeighbors(cellId, neighbors);

		// relax the noise cluster restrictions
		if ( cell->count + neighborCount < minPts ) {
			cell->clusterId = 0; // Mark as noise
			continue; // Not enough neighbors to form a cluster
		}

		//if (neighborCount < minPts) {
		//	cell->clusterId = 0; // Mark as noise
		//	continue; // Not enough neighbors to form a cluster
		//

		// Start a new cluster
		clusterId++;
		cell->clusterId = clusterId; 
		std::vector<int64_t> seedList;
		for (int i = 0; i < neighbors.size(); i++)
		{
			auto *nbCell = grid.getCell(neighbors[i]);

			// regain the noise point as border point, but not put in the seedList
			if (nbCell->clusterId == 0) {
				nbCell->clusterId = clusterId; // regain the noise point as border point
				continue;
			}

			// only unvisited points can be seed points
			else if (-1 == nbCell->clusterId) {
				seedList.push_back(neighbors[i]);
				nbCell->clusterId = -2;	// in the seedList, but not in the cluster yet
			}
		}

		// Expand the cluster
		// consuming all points in the seedList
		for (int k = seedList.size() - 1; k >= 0; --k)
		{
			auto* seedCell = grid.getCell(seedList[k]);

			// Already visited?
			if (seedCell->clusterId > -1)
				continue;

			// Mark the point as part of the cluster
			seedCell->clusterId = clusterId;
			std::vector<int64_t> seedNeighbors;
			int seedNeighborsCount = grid.GetNeighbors(seedList[k], seedNeighbors);

			seedList.pop_back();
			if ( seedCell->count + seedNeighborsCount >= minPts )
			{
				// only add new neighbors if they are not already in the seedList
				for (int i = 0; i < seedNeighbors.size(); i++)
				{
					auto* seedNbCell = grid.getCell(seedNeighbors[i]);

					if (seedNbCell->clusterId == 0) {
						seedNbCell->clusterId = clusterId; // regain the noise point as border point
						continue;
					}
					// only unvisited points can be seed points
					else if (-1 == seedNbCell->clusterId) {
						seedList.push_back(seedNeighbors[i]);
						seedNbCell->clusterId = -2;			// in the seedList, but not in the cluster yet
					}
				}

				k = seedList.size();
			}
		}		
	}

	// Mark the points to their clusters according to the grid cells
	grid.MarkPointsCluster(cPoints);

	// Calculate the center for each cluster
	std::map<int, std::vector<int> > mClusterPointIndicess;
	for (int i = 0; i < nPoints; i++)
	{
		int id = cPoints[i].clusterId;
		if (id > 0) {
			mClusterPointIndicess[id].push_back(i);
		}
	}

	for (const auto& cluster : mClusterPointIndicess)
	{
		GetCenterAndRadius(cPoints, cluster.second, dbScanCenters);

		vvClusterPoints.emplace_back(std::move(cluster.second));
	}

	return vvClusterPoints.size();
}


// Point count in grids based clustering algorithm
int PCGrid_2D(std::vector<cLasPOINT>& cPoints, double eps, std::vector<clusterCENTER>& dbScanCenters, std::vector<std::vector<int> >& vvClusterPoints)
{
	int nPoints = cPoints.size();
	if (nPoints <= 0 || cPoints.empty()) {
		std::cout << "No points to process." << std::endl;
		return 0;
	}

	// 1. Determine the grid size
	double gridSize = eps; // / sqrt(2); // or eps / 2; depends on how strict you want to be

	// 2. Find the bounding box of the data
	double minX = std::numeric_limits<double>::max();
	double minY = std::numeric_limits<double>::max();
	double maxX = std::numeric_limits<double>::lowest();
	double maxY = std::numeric_limits<double>::lowest();

	for (const auto& point : cPoints) {
		minX = std::min(minX, (double)point.pt.x());
		minY = std::min(minY, (double)point.pt.y());
		maxX = std::max(maxX, (double)point.pt.x());
		maxY = std::max(maxY, (double)point.pt.y());
	}

	// 3. Calculate the number of grid cells in each dimension
	int nX = (int)std::ceil((maxX - minX) / gridSize);
	int nY = (int)std::ceil((maxY - minY) / gridSize);

	// 4. Init PointGrid2D and register points
	PointGrid2D grid(minX, minY, nX, nY, gridSize);
	grid.RegisterPoints(cPoints);

	// 5. sort cells based on point count in each cell
	std::multimap<int, int64_t> sortedCells;  // <point count, cell id>
	grid.GetSortedCells(sortedCells);

	// 6. set point count noise threshold
	int minPts = 15;

	if( false ) {	
		int minCount, maxCount;
		float stdGridPts;
		int meanGridPts = grid.GetAverageCount( &stdGridPts, &minCount, &maxCount);

		fprintf(stderr, "PCGrid_2D GridPts mean = %d, std = %.1f, min = %d, max = %d\n", meanGridPts, stdGridPts, minCount, maxCount);

		minPts = meanGridPts - (int)stdGridPts; 
		if (minPts < 15)
			minPts = 15;
	}
	else {
		auto iter = sortedCells.begin();
		int n = sortedCells.size() / 5;
		for (int i = 0; i < n; i++)
			++iter;
		
		minPts = iter->first;

		std::cout << "  " << "PCGrid_2D has " << (static_cast<int>(sortedCells.size()) - n) << " sorted cells, minPts = " << minPts << std::endl;
	}

	// 7. clustering every cell using region growing
	int clusterId = 0;
	for (auto iter = sortedCells.rbegin(); iter != sortedCells.rend(); ++iter)
	{
		int count = iter->first;
		int64_t cellId = iter->second;
		GridCell *cell = grid.getCell(cellId);

		// has been visited;
		if (-1 != cell->clusterId)
			continue;

		// ignore cells whose point count is less than minPts
		if ( count < minPts) {
			cell->clusterId = 0; // Mark as noise
			continue; // Not enough neighbors to form a cluster
		}

		// looking for neighbors of cellId
		std::vector<int64_t> neighbors;
		int neighborCount = grid.GetNeighbors(cellId, neighbors);

		// anchor cell
		clusterId++;
		cell->clusterId = clusterId;

		// construct seedList
		std::vector<int64_t> seedList;
		for (int64_t nbId : neighbors)
		{
			auto *nbCell = grid.getCell(nbId);

			// regain the noise point as border point, but not put in the seedList
			if (nbCell->clusterId == 0 || nbCell->count < minPts ) {
				nbCell->clusterId = clusterId; // regain the noise point as border point
				continue;
			}
			// only unvisited points can be seed points
			else if (-1 == nbCell->clusterId) {
				seedList.push_back(nbId);
				nbCell->clusterId = -2;	// in the seedList, but not in the cluster yet
			}
		}

		// Expand the cluster
		// consuming all points in the seedList
		while (!seedList.empty()) {
        	int64_t seedId = seedList.back();
        	seedList.pop_back();
        	auto* seedCell = grid.getCell(seedId);

			// continue if already visited
			if (seedCell->clusterId > -1) continue;

			// update clusterId of seedCell
			seedCell->clusterId = clusterId;

        	// find neighbors of seedCell and put in seedLists
			std::vector<int64_t> seedNeighbors;
			grid.GetNeighbors(seedId, seedNeighbors);

        	for (int64_t seedNbId : seedNeighbors) {
            	auto* seedNbCell = grid.getCell(seedNbId);

				if (seedNbCell->clusterId == 0 || seedNbCell->count < minPts) {
					seedNbCell->clusterId = clusterId; // regain the noise point as border point
					continue;
				}
				
				// only unvisited points can be seed points
				else if (-1 == seedNbCell->clusterId) {
					seedList.push_back(seedNbId);
					seedNbCell->clusterId = -2;			// in the seedList, but not in the cluster yet
				}
			}			
		}
	}

	// Mark the points to their clusters according to the grid cells
	grid.MarkPointsCluster(cPoints);

	// Calculate the center for each cluster
	std::map<int, std::vector<int> > mClusterPointIndicess;
	for (int i = 0; i < nPoints; i++)
	{
		int id = cPoints[i].clusterId;
		if (id > 0) {
			mClusterPointIndicess[id].push_back(i);
		}
	}

	for (const auto& cluster : mClusterPointIndicess)
	{
		GetCenterAndRadius(cPoints, cluster.second, dbScanCenters);

		vvClusterPoints.emplace_back(std::move(cluster.second));
	}

	return dbScanCenters.size();
}


struct ConnectTowerPair
{
	int	  nInsideLineCells;		// 
	int	  nSamples;				// number of samples in the connection, related to the distance between two towers

	float X0, Y0;				// start point
	float cosA, sinA;			// direction 

	ConnectTowerPair() : nInsideLineCells(0), nSamples(0), X0(0), Y0(0), cosA(0), sinA(0) {};

	void Init(const Eigen::Vector3f& p0, const Eigen::Vector3f& p1, int nSamples)
	{
		X0 = p0.x();
		Y0 = p0.y();
		double dX = p1.x() - X0;
		double dY = p1.y() - Y0;

		double len = sqrt(dX * dX + dY * dY);

		cosA = dX / len;
		sinA = dY / len;
	}

	void Point2Line(float X, float Y, float* offset, float* dist)
	{
		double dX = X - X0;
		double dY = Y - Y0;

		*offset =  cosA * dX + sinA * dY;
		*dist   = -sinA * dX + cosA * dY;
	}
};



float EstimateMaxConnectionLength (std::vector<clusterCENTER>& towerCenters )
{
	int nTowers = towerCenters.size();
	int n1 = nTowers - 1;

	double lenSum = 0, len2Sum = 0;
	for( int i = 0; i < n1; ++i)
	{
		// Calculate the distance between two consecutive towers
		double dx = towerCenters[i + 1].pt.x() - towerCenters[i].pt.x();
		double dy = towerCenters[i + 1].pt.y() - towerCenters[i].pt.y();
		double s = dx * dx + dy * dy;

		lenSum += sqrt(s);
		len2Sum += s;
	}

	double avLen = lenSum / n1;
	float stdLen = sqrt(len2Sum / n1 - avLen*avLen);

	// Estimate the maximum connection length as 5 times the standard deviation
	return 5*stdLen;
}

extern bool g_bPrintDetails;
void PrintConnectionMap(const std::vector<ConnectTowerPair>& connectMap, int nTowers)
{
	if (g_bPrintDetails) {
		auto iter = connectMap.begin();
		for (int iTower = 0; iTower < nTowers; ++iTower)
		{
			printf("line %3d: ", iTower);
			for (int jTower = 0; jTower < nTowers; ++jTower)
			{
				if (iTower == jTower)
					printf("*");
				else {
					if( iter->nInsideLineCells > 0 )
						printf("8", iter->nInsideLineCells);
					else
						printf(" ");

					//printf("%2d ", 100 * iter->nInsideLineCells / (iter->nSamples - 2));
				}
				iter++;
			}
			printf("\n");
		}
	}
}

void FindLongestConnection(std::vector<ConnectTowerPair>& connectMap, std::vector<clusterCENTER>& towerCenters)
{
	int nTowers = towerCenters.size();
	if (nTowers <= 0 || connectMap.empty()) {
		std::cerr << "No towers or connections to process." << std::endl;
		return;
	}

	// construct the symmetric connection map
	for (int iTower = 0; iTower < nTowers; ++iTower)
	{
		for (int jTower = iTower + 1; jTower < nTowers; ++jTower)
		{
			int cntIndex = iTower * nTowers + jTower;
			int cntIndexSym = jTower * nTowers + iTower;

			connectMap[cntIndexSym] = connectMap[cntIndex];
		}
	}

	std::cout << "Cnnection Map: initial" << std::endl;
	PrintConnectionMap(connectMap, nTowers);

	// find the longest connection path
	std::multimap<int, std::vector<int> > mPaths;
	for (int startTower = 0; startTower < nTowers; ++startTower)
	{
		for (int nextTower = startTower + 1; nextTower < nTowers; ++nextTower)
		{
			std::vector<int> towers;
			towers.push_back(startTower);

			int traceTower = startTower;
			int candidateTower = nextTower;

			while (candidateTower < nTowers)
			{
				int cntIndex = traceTower * nTowers + candidateTower;
				if ( connectMap[cntIndex].nInsideLineCells > 0 )
				{
					towers.push_back(candidateTower);
					traceTower = candidateTower;
				}
				candidateTower++;
			}

			if (towers.size() > 1) mPaths.insert(std::make_pair(towers.size(), towers));
		}		
	}
	
	std::vector<int> &longestPath = mPaths.rbegin()->second;
	std::vector<clusterCENTER> tempTowerCenter;

	tempTowerCenter.reserve(longestPath.size());
	for (int i = 0; i < longestPath.size(); ++i)
	{
		int index = longestPath[i];
		tempTowerCenter.push_back(towerCenters[index]);
	}

	towerCenters = std::move(tempTowerCenter);
}


//
// check if the towers are in the specified power line corridor
// 1.Register power line points to the grid cells 
// 2.build the kdtree for the grid cells
// 3.filter the power towers with the power line points around the towers, there are 3 cases:
//  (1) the tower is the first or last tower of a power line points
//  (2) the tower is in the middle of a power line points, but the power line points are not connected to the tower
//  (3) the tower is outside the power line corridor, the power tower are at one side of power line points
// 
// 4.For each pair of towers, check their connection is consistent with power line points and update the grid cells id with nearest power connection
//  (1) indexing the line pcgrid point near the mid-line of the tower pair
//  (2) compute the distance and offset from the grid cell to the connection line
//  (3) check the inlier rate to determine if the tower pair is connected
//
void CheckTowersWithLinePoints( std::vector<cLasPOINT>& lineLasPoints, float gridSize, std::vector<clusterCENTER>& towerCenters, const std::string& outDir )
{
	std::cout << "Checking towers with line points..." << std::endl;
	
	// if mainline points num is less than 100, raise error
	int nLinePoints = lineLasPoints.size();
	if (nLinePoints <= 100 * towerCenters.size()) {
		fprintf(stderr, "No points to process.\n");
		return;
	}

	// 1.Register power line points to the grid cells 
	// 1.1 Find the bounding box of the data
	double minX = std::numeric_limits<double>::max();
	double minY = std::numeric_limits<double>::max();
	double minZ = std::numeric_limits<double>::max();
	double maxX = std::numeric_limits<double>::lowest();
	double maxY = std::numeric_limits<double>::lowest();
	double maxZ = std::numeric_limits<double>::lowest();

	{
		int n = lineLasPoints.size();
		double avZ = 0;
		double sumZ2 = 0;
		for (const auto& point : lineLasPoints) {
			minX = std::min(minX, (double)point.pt.x());
			minY = std::min(minY, (double)point.pt.y());
			minZ = std::min(minZ, (double)point.pt.z());

			maxX = std::max(maxX, (double)point.pt.x());
			maxY = std::max(maxY, (double)point.pt.y());
			maxZ = std::max(maxZ, (double)point.pt.z());

			avZ += (double)point.pt.z() / n;
			sumZ2 += (double)point.pt.z() * point.pt.z() / n;
		}
		double stdZ = sqrt(sumZ2 - avZ * avZ);

		// 	prevent outlers
		if (minZ < avZ - 3.5 * stdZ)
			minZ = avZ - 3.5 * stdZ;
		if (maxZ > avZ + 3.5 * stdZ)
			maxZ = avZ + 3.5 * stdZ;
	}

	// 1.2. Calculate the number of grid cells in each dimension
	int nX = (int)std::ceil((maxX - minX) / gridSize);
	int nY = (int)std::ceil((maxY - minY) / gridSize);
	int nZ = (int)std::ceil((maxZ - minZ) / gridSize);


	PointGrid2D grid(minX, minY, nX, nY, gridSize);
	grid.RegisterPoints(lineLasPoints);

	std::vector<int> vCellIds;
	int nCells = grid.GetCellIds(vCellIds);

	std::cout << "  Found " << nCells << " cells centers for line grids" << std::endl;

	if (nCells < towerCenters.size())
	{
		std::cerr << "Not enough grid cells to process." << std::endl;
		return;
	}

	// 2.build the kdtree for the grid cells
	std::vector<cLasPOINT> cCellCenterPoints;

	{
		cCellCenterPoints.resize(nCells);
		for (int i = 0; i < nCells; i++)
		{
			GridCell *cellPtr = grid.getCell( vCellIds[i] );
			cCellCenterPoints[i].pt << cellPtr->Xw, cellPtr->Yw, 0.0f;  // Eigen::Vector3f init
			cCellCenterPoints[i].clusterId = -1; // unvisited
		}
	}

	using PointCloud_C = PointCloudAdaptor<cLasPOINT>;
	PointCloud_C pointCloud_C(cCellCenterPoints);
	typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<float, PointCloud_C >,
		PointCloud_C,
		3 /* dim */
	> KDTree;

	KDTree kdIndex(3, pointCloud_C, { 10 /* max leaf */ });
	kdIndex.buildIndex();

	// 3.filter the power towers with the power line points around the towers, there are 3 cases:
	//  (1) the tower is the first or last tower of a power line points
	//  (2) the tower is in the middle of a power line points, but the power line points are not connected to the tower
	//  (3) the tower is outside the power line corridor, the power tower are at one side of power line points
	{
		SortTowerCentersWithDistance(towerCenters, lineLasPoints[0].pt);
		std::cout << "  " << towerCenters.size() << " towers are sorted based on tower distance" << std::endl;
		PrintTowerList(towerCenters);

		nanoflann::SearchParameters params(0, false);
		double eps2 = 6 * 6;

		std::vector<clusterCENTER> towerCentersInsideCorridor;
		towerCentersInsideCorridor.reserve(towerCenters.size());
		for (int iTower = 0; iTower < towerCenters.size(); ++iTower)
		{
			float xyz[3];

			xyz[0] = towerCenters[iTower].pt.x();
			xyz[1] = towerCenters[iTower].pt.y();
			xyz[2] = 0;	// Z0Max + k * dZ;

			// Find neighbors using nanoflann
			std::vector<nanoflann::ResultItem<unsigned int, float>> ret_matches;

			const size_t nMatches = kdIndex.radiusSearch( xyz, eps2, ret_matches, params );

			std::vector<cLasPOINT> cPoints;
			std::vector<int> cPointIndices;

			cPoints.resize(nMatches);
			cPointIndices.resize(nMatches);
			for (size_t i = 0; i < nMatches; ++i)
			{
				int cellId = vCellIds[ret_matches[i].first];
				GridCell* cell = grid.getCell(cellId);
				
				cPoints[i].pt.x() = cell->Xw;
				cPoints[i].pt.y() = cell->Yw;
				cPoints[i].pt.z() = cell->Zw;
				cPointIndices[i] = i;
			}

			// if there are nearby line grid cells
			if (nMatches >= 1 ) {
				if (iTower > 0 && iTower < towerCenters.size() - 1) {
					Eigen::Vector3f center;
					double eigenValue[2], axis[2];
					float axesLen[2];

					// computePrincipalAxes_2D(cPoints, cPointIndices, center, eigenValue, axis);
					computeAxesLen_2D(cPoints, cPointIndices, center, eigenValue, axis, axesLen);

					if (center.z() < towerCenters[iTower].maxZ)
					{
						double dX = xyz[0] - center.x();
						double dY = xyz[1] - center.y();

						float offset = axis[0] * dX + axis[1] * dY; // offset to the main axis
						float dist = -axis[1] * dX + axis[0] * dY;  // distance to the main axis

						if (fabs(dist) < axesLen[1])
						{
							towerCentersInsideCorridor.push_back(towerCenters[iTower]);
						}
					}
				}

				// the first and late tower will be kept
				else {
					towerCentersInsideCorridor.push_back(towerCenters[iTower]);
				}
			}
		}

		std::cout << "  Found " << towerCentersInsideCorridor.size() << " towers inside the corridor" << std::endl;

		towerCenters = std::move(towerCentersInsideCorridor);
		PrintTowerList(towerCenters);

		if ( g_bSaveFile ) {
			std::string lasCenterInsideFile = outDir + "/InsideTowerCenters.las";
			SaveCenters2LasFile(lasCenterInsideFile.c_str(), towerCenters);
		}
	}


	// 4.For each pair of towers, check their connection is consistent with power line points and update the grid cells id with nearest power connection
	// determine the connection length threshold
	float cntLenTh = 3 * EstimateMaxConnectionLength(towerCenters);

	// We don't need sorted results for radius search
	// nanoflann::SearchParameters params(0, false);
	uint32_t ret_matches[1];
	float  ret_dist2[1];

	int nTowers = towerCenters.size();
	std::vector<ConnectTowerPair> connectMap;
	connectMap.resize(nTowers * nTowers);

	for (int iTower = 0; iTower < nTowers; ++iTower)
	{
		for (int jTower = iTower + 1; jTower < nTowers; ++jTower)
		{
			int cntIndex = iTower * nTowers + jTower;

			// sample the connection line between two towers
			// skip the first and last tower points
			double dX = (towerCenters[jTower].pt.x() - towerCenters[iTower].pt.x());
			double dY = (towerCenters[jTower].pt.y() - towerCenters[iTower].pt.y());

			float  Z0Max = towerCenters[iTower].getMaxZ() + 3;		// add 3 meters tolerance
			double dZ = towerCenters[jTower].getMaxZ() - Z0Max;

			double len = sqrt(dX * dX + dY * dY);
			if ( len > 2*cntLenTh )
				continue;

			int nSamples = len / gridSize;

			nSamples = std::max(nSamples, 5);

			connectMap[cntIndex].Init(towerCenters[iTower].pt, towerCenters[jTower].pt, nSamples);

			dX /= nSamples;
			dY /= nSamples;
			dZ /= nSamples;

			double radius_i = towerCenters[iTower].get2DRadius();
			double radius_j = towerCenters[jTower].get2DRadius();

			double dRadius = (radius_j - radius_i) / nSamples;
			double dOffset = sqrt(dX * dX + dY * dY);

			// continue if the two towers are too close
			if (len < 2 * radius_i || len < 2 * radius_j) {
				continue;
			}

			std::vector<float> insideCellElevations;  // line grid cell height 
			insideCellElevations.reserve(nSamples - 2);
			{
				float xyz[3];

				for (int k = 1; k < nSamples - 1; ++k)
				{
					xyz[0] = towerCenters[iTower].pt.x() + k * dX;
					xyz[1] = towerCenters[iTower].pt.y() + k * dY;
					xyz[2] = 0;	// Z0Max + k * dZ;

					float offset_cnt = k * dOffset;
					float radius = radius_i + k * dRadius;

					// seaech for line cell
					kdIndex.knnSearch(xyz, 1, ret_matches, ret_dist2);

					int cellId = vCellIds[ret_matches[0]];
					GridCell* cell = grid.getCell(cellId);

					// line is under the max height of the towers
					if( (Z0Max + k * dZ) < cell->Zw )
						continue;

					// get the offset and distance to the connection line with respect to tower i
					float dx_cell, dy_cell;
					connectMap[cntIndex].Point2Line(cell->Xw, cell->Yw, &dx_cell, &dy_cell);

					// in the radius
					if ( fabs(dy_cell) < radius) 
					{
						dx_cell -= offset_cnt;
						// within the nearby offset
						if ( fabs(dx_cell) < 2*gridSize )
						{
							insideCellElevations.push_back( cell->maxZ );
						}
					}
				}
			}

			int nInsideCells = insideCellElevations.size();
			float insideRate = (float)nInsideCells / (nSamples - 2);

			// require more than 90% of the sampled cells are inside the line corridor
			if (insideRate > 0.9 ) 
			{
				connectMap[cntIndex].nInsideLineCells = nInsideCells;
				connectMap[cntIndex].nSamples = nSamples;
			}
		}
	}

	FindLongestConnection(connectMap, towerCenters);

}


void FilterLinePointsWithPCGrid(std::vector<cLasPOINT>& lineLasPoints, double eps, int minPts, std::vector<clusterCENTER>& clusterCenters, std::vector<cLasPOINT>& mainLinePoints)
{

	std::cout << "Detecting Mainline Powerline using DBScan from " << int(lineLasPoints.size()) << " points..."<< std::endl;

	// point indexing for each cluster
	std::vector<std::vector<int> > vvClusterPointIndices;
	int nLineClusters = PCGrid_3D(lineLasPoints, 3.0, clusterCenters, vvClusterPointIndices);

	// point count - cluster index
	std::multimap<int, int> mClusterPointIndices;
	for (int i = 0; i < nLineClusters; ++i)
	{
		mClusterPointIndices.insert(std::make_pair(vvClusterPointIndices[i].size(), i));
	};

	std::vector<int> vNumOfClusterPoints;
	{
		for (auto iter = mClusterPointIndices.rbegin(); iter != mClusterPointIndices.rend(); ++iter)
		{
			//fprintf(stderr, "points of line cluster %d: %d\n", n, iter->first);
			vNumOfClusterPoints.push_back(iter->first);
		}
	}

	std::cout << "  " << vNumOfClusterPoints.size() << " power line groups are found" << std::endl;

	int n = 0;
	for (int i = 0; i < vNumOfClusterPoints.size() - 1; i++)
	{
		n += vNumOfClusterPoints[i];
		if(i == vNumOfClusterPoints.size() - 2 || vNumOfClusterPoints[i] > 3 * vNumOfClusterPoints[i + 1] || vNumOfClusterPoints[i + 1] < 100 )
			break;
	}

	mainLinePoints.resize( n );

	n = 0;
	int i = 0;
	for ( auto iter = mClusterPointIndices.rbegin(); iter != mClusterPointIndices.rend(); ++iter, ++i)
	{
		int iCluster = iter->second;
		std::vector<int>& vPointIndices = vvClusterPointIndices[iCluster];

		for (const auto& index : vPointIndices)
		{
			mainLinePoints[n].pt = lineLasPoints[index].pt;
			n++;
		}

		std::cout << "\t" << vNumOfClusterPoints[i] << " power line points are appended" << std::endl;
		if ( i == vNumOfClusterPoints.size()-2 || vNumOfClusterPoints[i] > 3 * vNumOfClusterPoints[i + 1] || vNumOfClusterPoints[i + 1] < 100)
		{
			break;
		}		
	}
}
