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


void PrintTowerList(const std::vector<clusterCENTER>& towerCenters);

int GetTowerCenter(std::vector<cLasPOINT>& towerPoints, float minHei, std::vector<clusterCENTER>& towerCenters, const bool debug);

void SortTowerCentersWithDistance(std::vector<clusterCENTER>& towerCenters, Eigen::Vector3f curCenter);

void InitLineSegs(const std::vector<clusterCENTER>& towerCenters, std::vector<lineSEG>& lineSegs);

int GroupPowerLinePoints(const std::vector<cLasPOINT>& lineLasPoints, const std::vector<clusterCENTER>& towerCenters, const std::vector<lineSEG> &lineSegs, std::vector<std::vector<int> >& groupedLinePts);

