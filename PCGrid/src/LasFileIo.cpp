#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include <liblas/liblas.hpp>
#include "base.h"


static bool bOffsetSet = false;
static double offsetX = 0;
static double offsetY = 0;
static double offsetZ = 0;

static double lasCoordScales[3];

int ReadLasFile(const std::string& lasFile, std::vector<cLasPOINT>& lasPoints)
{
    if (lasFile.empty()) {
        throw std::invalid_argument("Invalid LAS file path: " + lasFile);
    }

    std::ifstream ifs(lasFile, std::ios::in | std::ios::binary);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open LAS file: " + lasFile);
    }

    liblas::ReaderFactory f;
    liblas::Reader reader = f.CreateWithStream(ifs);
    const liblas::Header& header = reader.GetHeader();

    std::cout << "Reading LAS file: " << lasFile << std::endl;
    std::cout << "  Number of points: " << header.GetPointRecordsCount() << std::endl;

    lasCoordScales[0] = header.GetScaleX();
    lasCoordScales[1] = header.GetScaleY();
    lasCoordScales[2] = header.GetScaleZ();

    offsetX = header.GetOffsetX();
    offsetY = header.GetOffsetY();
    offsetZ = header.GetOffsetZ();

    lasPoints.clear();
    lasPoints.reserve(header.GetPointRecordsCount());

    while (reader.ReadNextPoint()) {
        const liblas::Point& p = reader.GetPoint();

        cLasPOINT point;
        point.pt.x() = p.GetX() - offsetX;
        point.pt.y() = p.GetY() - offsetY;
        point.pt.z() = p.GetZ();

        lasPoints.push_back(point);
    }

    std::cout << "  Total points loaded: " << lasPoints.size() << std::endl;
    return static_cast<int>(lasPoints.size());
}

void WriteLasFile(const std::string& lasFile, const std::vector<Eigen::Vector3f>& points)
{
    if (lasFile.empty()) {
        std::cerr << "Invalid LAS file path: " << lasFile << std::endl;
        return;
    }

    std::ofstream ofs(lasFile, std::ios::out | std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "Failed to open LAS file for writing: " << lasFile << std::endl;
        return;
    }

    liblas::Header header;
    header.SetDataFormatId(liblas::ePointFormat0);
    header.SetScale(lasCoordScales[0], lasCoordScales[1], lasCoordScales[2]);
    header.SetOffset(offsetX, offsetY, offsetZ);
    header.SetPointRecordsCount(points.size());

    liblas::Writer writer(ofs, header);

    for (const auto& p : points) {
        liblas::Point point(&header);
        point.SetCoordinates((double)p.x() + offsetX, (double)p.y() + offsetY, p.z());
        writer.WritePoint(point);
    }

    ofs.close();
    std::cout << "Saved LAS file: " << lasFile << " (" << points.size() << " points)" << std::endl;
}


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


void SavePoints2LasFile(const std::string& lasFile, const std::vector<cLasPOINT>& lasPoints)
{
    std::vector<Eigen::Vector3f> pts;
    pts.reserve(lasPoints.size());
    for (const auto& p : lasPoints)
        pts.push_back(p.pt);

    WriteLasFile(lasFile, pts);
}

void SavePoints2LasFile(const std::string& lasFile, const std::vector<Eigen::Vector3f>& lasPoints)
{
    WriteLasFile(lasFile, lasPoints);
}

void SaveCenters2LasFile(const std::string& lasFile, const std::vector<clusterCENTER>& lasPoints)
{
    std::vector<Eigen::Vector3f> pts;
    pts.reserve(lasPoints.size());
    for (const auto& c : lasPoints)
        pts.emplace_back(c.pt.x(), c.pt.y(), c.maxZ);

    WriteLasFile(lasFile, pts);
}

void SaveGroupPoints2LasFile(const std::string& outDir,
                             const std::vector<cLasPOINT>& lineLasPoints,
                             const std::vector<std::vector<int>>& groupedLinePts)
{
    if (outDir.empty()) {
        std::cerr << "Invalid output directory." << std::endl;
        return;
    }

    for (size_t i = 0; i < groupedLinePts.size(); ++i)
    {
        std::string groupFileName = outDir + "/group_" + std::to_string(i) + ".las";

        std::vector<Eigen::Vector3f> pts;
        pts.reserve(groupedLinePts[i].size());

        for (int idx : groupedLinePts[i]) {
            if (idx >= 0 && idx < (int)lineLasPoints.size())
                pts.push_back(lineLasPoints[idx].pt);
        }

        if (pts.empty()) {
            std::cout << "No points in group " << i << ", skipping file: " << groupFileName << std::endl;
            continue;
        }

        WriteLasFile(groupFileName, pts);
    }
}
