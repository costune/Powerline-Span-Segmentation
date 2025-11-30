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
#include "LasFileIo.h"
#include "PCGrid.h"
#include "Las2PowerLine.h"

#ifdef WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif



bool createDirectory(const std::string& path) {
#ifdef WIN32
    // Windows
    if (CreateDirectory(path.c_str(), NULL) || GetLastError() == ERROR_ALREADY_EXISTS) {
        return true;
    } else {
        return false;
    }
#else
    // Linux / macOS
    if (mkdir(path.c_str(), 0777) == 0 || errno == EEXIST) {
        return true;
    } else {
        return false;
    }
#endif
}


static void usage()
{
	printf("usage: Las2PowerLine lineLasFile towerLasFile outDir [--verbose]\n");
	exit(1);
}


int main(int argc, char **argv){

	if (argc != 4 && argc != 5)
	{
		usage();
		return 0;
	}

	const bool debug = (argc == 5) && (std::string(argv[4]) == "--verbose");

	std::string lineLasFile = argv[1];
	std::string towerLasFile = argv[2];
	std::string outDir = argv[3];

	{

		// Create output Dir
		if (createDirectory(outDir)) {
			std::cout << "Directory " + outDir + " is created successfully." << std::endl;
		} else {
			std::cerr << "Fail to create " + outDir << std::endl;
		}


		auto startTime = std::chrono::high_resolution_clock::now();

		// Read line and tower pointclouds
		std::vector<cLasPOINT> lineLasPoints;
		std::vector<cLasPOINT> towerLasPoints;

		int nLineLasPoints0 = 0;
		try {
			nLineLasPoints0 = ReadLasFile(lineLasFile, lineLasPoints);
		} catch (const std::exception& e) {
			std::cerr << "Error reading line LAS file: " << e.what() << std::endl;
			return -1;
		}
		int nLineLasPoints1 = EraseRepeatedPoints(lineLasPoints);
		std::cout << "  "
			<< nLineLasPoints1
        	<< " valid line points are read, "
			<< std::fixed << std::setprecision(1)
			<< 100.0 * (nLineLasPoints0 - nLineLasPoints1) / nLineLasPoints0
			<< "% points are erased"
			<< std::endl;

		int nTowerLasPoints0 = 0;
		try {
			nTowerLasPoints0 = ReadLasFile(towerLasFile, towerLasPoints);
		} catch (const std::exception& e) {
			std::cerr << "Error reading tower LAS file: " << e.what() << std::endl;
			return -1;
		}
		int nTowerLasPoints1 = EraseRepeatedPoints(towerLasPoints);
		std::cout << "  "
			<< nTowerLasPoints1
        	<< " valid tower points are read, "
			<< std::fixed << std::setprecision(1)
			<< 100.0 * (nTowerLasPoints0 - nTowerLasPoints1) / nTowerLasPoints0
			<< "% points are erased"
			<< std::endl;

		auto endTime_reading = std::chrono::high_resolution_clock::now();
		{
			auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(endTime_reading - startTime);
			std::cout << "Time taken to read las point: "
				<< std::fixed << std::setprecision(2)
				<< duration.count()
				<< " seconds"
				<< std::endl;
		}

		// Main procedure
		if (nLineLasPoints1 > 100 && nTowerLasPoints1 > 100)
		{
			// 1.1 Tower center extraction using PCGrid 2D
			std::vector<clusterCENTER> towerCenters;
			int nTowers = GetTowerCenter(towerLasPoints, 3.0, towerCenters, debug);  // set min height of tower to 3.0
			if (nTowers < 2)
			{
				std::cout << "Not enough towers detected, please check the tower LAS file." << std::endl;
				goto end_l;
			}

			std::cout << "Found " << (int)towerCenters.size() << " tower centers with 2D grid based DBScan." << std::endl;
			PrintTowerList(towerCenters);

			if (debug) {
				std::string lasCenterFile = outDir + "/InitialTowerCenters.las";
				SaveCenters2LasFile(lasCenterFile, towerCenters);
			}

			// 1.2 Mainline powerline extraction using PCGrid 3D
			std::vector<cLasPOINT> mainLineLasPoints;
			{
				std::vector<clusterCENTER> mainLineCenters;
				FilterLinePointsWithPCGrid(lineLasPoints, 3.0, 15, mainLineCenters, mainLineLasPoints);

				if (debug) {
					std::string mainLineLasFile = outDir + "/mainLines.las";
					SavePoints2LasFile(mainLineLasFile, mainLineLasPoints);
				}
			}

			// 1.3 The towers are paired and verified using power lines to check if there are sufficient power line points between them.
			CheckTowersWithLinePoints(mainLineLasPoints, 1.0, towerCenters, outDir, debug);
			std::cout << towerCenters.size() << " towers passed checking" << std::endl;
			PrintTowerList(towerCenters);

			if (debug) {
				std::string CheckedTowerCenterFile = outDir + "/CheckedTowerCenters.las";
				SaveCenters2LasFile(CheckedTowerCenterFile, towerCenters);
			}

			// 2. Power line segmentation between towers
			// 2.1 Init span according to tower centers
			std::vector<lineSEG> lineSegs;
			InitLineSegs(towerCenters, lineSegs);

			std::vector<std::vector<int>> vvSegLinePtIndices;

			// 2.2 根据点杆塔连线的垂直距离，对电力线进行分档
			int nGroups = GroupPowerLinePoints(mainLineLasPoints, towerCenters, lineSegs, vvSegLinePtIndices);
			
			if( debug )
				SaveGroupPoints2LasFile(outDir, mainLineLasPoints, vvSegLinePtIndices);

			{
				// Record end time
				auto endTime = std::chrono::high_resolution_clock::now();

				// Calculate the duration
				auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(endTime_reading - startTime);

				std::cout << "Power Line Points Segmentaion, time taken: " << duration.count() << " seconds" << std::endl;
			}
		}
	}
end_l:
	return 0;
}

