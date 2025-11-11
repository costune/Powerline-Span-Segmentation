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
#include "DenseCluster.h"


void process_points(const std::string& lineLasFile, const std::string& towerLasFile, const std::string& outDir);