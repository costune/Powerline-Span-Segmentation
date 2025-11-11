#ifndef BASE_H
#define BASE_H

#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <Eigen/Dense>


// clusterId:
//  -2: Has been added to the neighborhood point list
//	-1: Unvisited
//   0：Noise points or boundary points
// > 1: True cluster ID
struct cLasPOINT
{
    Eigen::Vector3f pt = Eigen::Vector3f::Zero();
    int clusterId = -1;
    int nextPt = -1;

    cLasPOINT() = default;
    explicit cLasPOINT(const Eigen::Vector3f& c) : pt(c) {}
};


struct clusterCENTER
{
    Eigen::Vector3f pt = Eigen::Vector3f::Zero();
    int numOfPts = 0;
    float minZ = 0.0f;
    float maxZ = 0.0f;
    double eigenValues[3]   = {0.0, 0.0, 0.0};
    double eigenVectors[9]  = {0.0};
    float  axesLen[3]       = {0.0f, 0.0f, 0.0f};
    float  direction = 0.0f;

    clusterCENTER() = default;
    explicit clusterCENTER(const Eigen::Vector3f& c) : pt(c) {}

    float get2DRadius() const noexcept { return axesLen[0] * 0.5f; }
    float get3DRadius() const noexcept { return axesLen[0] * 0.5f; }
    float getHeight() const noexcept { return maxZ - minZ; }
	float getMaxZ() const noexcept { return maxZ; }
};


// custom data structure for nanoflann
template <typename T>
struct PointCloudAdaptor
{
    const std::vector<T>& points;

    explicit PointCloudAdaptor(const std::vector<T>& pts) : points(pts) {}

    // 返回点数量
    size_t kdtree_get_point_count() const noexcept { return points.size(); }

    // 返回单个点在指定维度的值
    double kdtree_get_pt(const size_t idx, const size_t dim) const noexcept
    {
        switch (dim)
        {
            case 0: return points[idx].pt.x();
            case 1: return points[idx].pt.y();
            default: return points[idx].pt.z();
        }
    }

    template <typename BBOX>
	bool kdtree_get_bbox(BBOX& bbox) const
	{
		if (points.empty()) return false;

		bbox[0].low = bbox[0].high = points[0].pt.x();
		bbox[1].low = bbox[1].high = points[0].pt.y();
		bbox[2].low = bbox[2].high = points[0].pt.z();
		for (size_t i = 1; i < points.size(); ++i)
		{
			bbox[0].low = bbox[0].low < points[i].pt.x() ? bbox[0].low : points[i].pt.x();
			bbox[0].high = bbox[0].high > points[i].pt.x() ? bbox[0].high : points[i].pt.x();
			bbox[1].low = bbox[1].low < points[i].pt.y() ? bbox[1].low : points[i].pt.y();
			bbox[1].high = bbox[1].high > points[i].pt.y() ? bbox[1].high : points[i].pt.y();
			bbox[2].low = bbox[2].low < points[i].pt.z() ? bbox[2].low : points[i].pt.z();
			bbox[2].high = bbox[2].high > points[i].pt.z() ? bbox[2].high : points[i].pt.z();
		}
		return true;
	}
};

struct  lineSEG
{
	float X0, Y0;

	// 方向
	float cosA, sinA;
	float len;

public:
	void InitLineSeg(const Eigen::Vector3f& from, const Eigen::Vector3f& to)
	{
		X0 = from.x();	Y0 = from.y();

		double dX = to.x() - from.x();
		double dY = to.y() - from.y();

		len = sqrt(dX * dX + dY * dY);
		cosA = dX / len;
		sinA = dY / len;
	}


	bool ProjecPoint2Line(const Eigen::Vector3f& pt, float* offsetX, float* dist)  const
	{
		if (len <= 0)
		{
			*offsetX = 0;
			*dist = 0;
			return false;
		}

		double dX = pt.x() - X0;
		double dY = pt.y() - Y0;

		*offsetX = cosA * dX + sinA * dY;
		*dist = -sinA * dX + cosA * dY;

		return true;
	}

	bool PointFromProject(float offsetX, float dist, float *X, float *Y) const
	{
		if (len <= 0)
		{
			*X = 0;
			*Y = 0;
			return false;
		}

		*X = X0 + offsetX * cosA - dist * sinA;
		*Y = Y0 + offsetX * sinA + dist * cosA;

		return true;
	}

	bool IsPointInSeg(const Eigen::Vector3f& pt, float* offsetX, float *dist) const
	{
		if( !ProjecPoint2Line(pt, offsetX, dist) )
			return false;

		if (*offsetX >= 0 && *offsetX < len)
			return true;
		
		return false;
	}
};





#endif
