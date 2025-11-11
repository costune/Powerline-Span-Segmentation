#ifndef LAS_FILE_IO_H
#define LAS_FILE_IO_H

#include <string>
#include <iostream>
#include <liblas/liblas.hpp>
#include "base.h"

int ReadLasFile(const std::string& lasFile, std::vector<cLasPOINT>& lasPoints);
int EraseRepeatedPoints(std::vector<cLasPOINT>& lasPoints);
// void WriteLasFile(const std::string& lasFile, const std::vector<Eigen::Vector3f>& points);
void SaveCenters2LasFile(const std::string& lasFile, const std::vector<clusterCENTER>& lasPoints);
void SavePoints2LasFile(const std::string& lasFile, const std::vector<cLasPOINT>& lasPoints);
void SavePoints2LasFile(const std::string& lasFile, const std::vector<Eigen::Vector3f>& lasPoints);
void SaveGroupPoints2LasFile(const std::string& outDir, const std::vector<cLasPOINT>& lineLasPoints, const std::vector<std::vector<int>>& groupedLinePts);

#endif // LAS_FILE_IO_H