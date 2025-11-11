#ifndef ORS_POINT_CLOUD_FEATURE_H
#define ORS_POINT_CLOUD_FEATURE_H

#include "DenseCluster.h"

void ComputeCovarianceMatrix_2x2(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double QQ[4] );
void ComputeCovarianceMatrix_3x3(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double QQ[9]);

void computePrincipalAxes_2D(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double eigenValue[2], double axis[2]);

void computeAxesLen_2D(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double eigenValue[2], double axis[2], float axisLen[2]);

void computePrincipalAxes_3D(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double eigenValue[3], double axis[9]);

void GetCenterAndRadius(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, std::vector<clusterCENTER>& clusterCenters);



#endif
