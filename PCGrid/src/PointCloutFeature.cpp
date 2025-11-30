#include "base.h"
#include "PointCloutFeature.h"


void computeCenter_Std(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, Eigen::Vector3f& std)
{
	double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
	double sumX2 = 0.0, sumY2 = 0.0, sumZ2 = 0.0;

	for (const auto& index : clusterPointIndices)
	{
		auto& p = cPoints[index].pt;
		sumX += p.x(); sumY += p.y();	sumZ += p.z();
		sumX2 += (double)p.x() * p.x();
		sumY2 += (double)p.y() * p.y();
		sumZ2 += (double)p.z() * p.z();
	}

	int n = clusterPointIndices.size();
	double Xc = sumX / n, Yc = sumY / n, Zc = sumZ / n;

	center.x() = Xc;	center.y() = Yc;	center.z() = Zc;

	std.x() = sqrt(sumX2 / n - Xc * Xc);
	std.y() = sqrt(sumY2 / n - Yc * Yc);
	std.z() = sqrt(sumZ2 / n - Zc * Zc);
}


void ComputeCovarianceMatrix_2x2(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double QQ[4])
{
	double sumX = 0.0, sumY = 0.0, sumZ = 0.0;

	for (const auto& index : clusterPointIndices)
	{
		auto& p = cPoints[index].pt;
		sumX += p.x(); sumY += p.y();	sumZ += p.z();
	}

	int n = clusterPointIndices.size();
	double Xc = sumX / n, Yc = sumY / n, Zc = sumZ/n;

	for (int i = 0; i < 4; i++)
		QQ[i] = 0;

	for (const auto& index : clusterPointIndices)
	{
		auto& p = cPoints[index].pt;
		double dX = p.x() - Xc;
		double dY = p.y() - Yc;

		QQ[0] += dX * dX;	QQ[1] += dX * dY;
							QQ[3] += dY * dY;
	}
	QQ[2] = QQ[1];

	for (int i = 0; i < 4; i++)
		QQ[i] /= n;

	center.x() = Xc; center.y() = Yc; center.z() = Zc;
}


void ComputeCovarianceMatrix_3x3(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double QQ[9])
{
	double sumX = 0.0, sumY = 0.0, sumZ = 0.0;

	for (const auto& index : clusterPointIndices)
	{
		auto& p = cPoints[index].pt;
		sumX += p.x(); sumY += p.y();	sumZ += p.z();
	}

	int n = clusterPointIndices.size();
	double Xc = sumX / n, Yc = sumY / n, Zc = sumZ / n;

	for (int i = 0; i < 9; i++)
		QQ[i] = 0;
	
	for (const auto& index : clusterPointIndices)
	{
		auto& p = cPoints[index].pt;
		double dX = p.x() - Xc;
		double dY = p.y() - Yc;
		double dZ = p.z() - Zc;

		QQ[0] += dX * dX;	QQ[1] += dX * dY;	QQ[2] += dX * dZ;
							QQ[4] += dY * dY;	QQ[5] += dY * dZ;
												QQ[8] += dZ * dZ;
	}
	QQ[3] = QQ[1];	
	QQ[6] = QQ[2];	QQ[7] = QQ[5];

	for (int i = 0; i < 9; i++)
		QQ[i] /= n;

	center.x() = Xc; center.y() = Yc; center.z() = Zc;
}

void computePrincipalAxes_2D(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f &center, double eigenValue[2], double axis[2])
{
	double QQ[4];
	
	ComputeCovarianceMatrix_2x2(cPoints, clusterPointIndices, center, QQ);

	// trace of matrix QQ
	double tr = QQ[0] + QQ[3];

	// det of matrix QQ
	double det = QQ[0] * QQ[3] - QQ[1] * QQ[2];

	// equation:  lambda*lambda - trace*lambda + det = 0;
	// 

	double m = tr / 2;	// (lambda1 + lambda2)/2  = tr/2
	double p = det;		// lambda1 * lambda2 = det 
	
	double sqrt_delta = sqrt(m * m - p);

	eigenValue[0] = m + sqrt_delta;
	eigenValue[1] = m - sqrt_delta;

	// eigen vector
	double a0 = QQ[0] - eigenValue[0],	a1 = QQ[1];
	double a2 = QQ[1], a3 = QQ[3] - eigenValue[0];

	axis[0] = 1;
	axis[1] = -a0/a1*eigenValue[0];

	// normalize
	double s = sqrt( axis[0] * axis[0] + axis[1] * axis[1] );

	axis[0] /= s;	axis[1] /= s;
}


void computeAxesLen_2D(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f& center, double eigenValue[2], double axis[2], float axesLen[2])
{
	computePrincipalAxes_2D(cPoints, clusterPointIndices, center, eigenValue, axis);

	float minOffset = 9999.0, maxOffset = -9999.0;
	float minDist = 9999.0, maxDist = -9999.0;

	// 易受噪声干扰，需要先过滤噪声
	// Calculate the lengths of the axes based on eigenvector
	for (const auto& index : clusterPointIndices)
	{
		auto& p = cPoints[index].pt;

		double dX = p.x() - center.x();
		double dY = p.y() - center.y();

		// Project the point onto the principle eigenvector
		float offset = (float)( dX * axis[0] + dY * axis[1]);
		float dist   = (float)(-dX * axis[1] + dY * axis[0]);

		if( offset < minOffset)
			minOffset = offset;

		if (offset > maxOffset)
			maxOffset = offset;

		if ( dist < minDist )
			minDist = dist;

		if (dist > maxDist)
			maxDist = dist;
	}

	axesLen[0] = (maxOffset - minOffset);
	axesLen[1] = (maxDist - minDist);
}

// Calculate eigenvalues and eigenvectors using the Jacobi iteration method
void Eigen3x3_Jacobi(double QQ[9], double eigenValue[3], double axis[9])
{
	// Initialize eigenvectors to identity matrix
	axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
	axis[3] = 0.0; axis[4] = 1.0; axis[5] = 0.0;
	axis[6] = 0.0; axis[7] = 0.0; axis[8] = 1.0;

	// Copy covariance matrix to a local variable that will be modified
	double A[9];
	for (int i = 0; i < 9; ++i) {
		A[i] = QQ[i];
	}

	// Perform Jacobi iterations
	for (int iter = 0; iter < 50; ++iter) // Iterate until convergence or max iterations reached
	{
		// Find the largest off-diagonal element
		int p = 0, q = 1;
		for (int i = 0; i < 3; ++i) {
			for (int j = i + 1; j < 3; ++j) {
				if (fabs(A[i * 3 + j]) > fabs(A[p * 3 + q])) {
					p = i;
					q = j;
				}
			}
		}

		// Calculate the rotation angle
		double theta;
		if (A[p * 3 + q] == 0.0) {
			theta = 0.0;
		}
		else {
			theta = 0.5 * atan2(2.0 * A[p * 3 + q], A[q * 3 + q] - A[p * 3 + p]);
		}

		// Calculate the rotation matrix
		double c = cos(theta);
		double s = sin(theta);

		// Update the covariance matrix
		for (int i = 0; i < 3; ++i)
		{
			double t1 = A[i * 3 + p];
			double t2 = A[i * 3 + q];
			A[i * 3 + p] = c * t1 - s * t2;
			A[i * 3 + q] = s * t1 + c * t2;
		}

		for (int i = 0; i < 3; ++i) {
			double t1 = A[p * 3 + i];
			double t2 = A[q * 3 + i];
			A[p * 3 + i] = c * t1 - s * t2;
			A[q * 3 + i] = s * t1 + c * t2;
		}

		// Update the eigenvector matrix
		for (int i = 0; i < 3; ++i) {
			double t1 = axis[i * 3 + p];
			double t2 = axis[i * 3 + q];
			axis[i * 3 + p] = c * t1 - s * t2;
			axis[i * 3 + q] = s * t1 + c * t2;
		}
	}

	// The eigenvalues are now on the diagonal of A
	eigenValue[0] = A[0];
	eigenValue[1] = A[4];
	eigenValue[2] = A[8];
}


// traditional Cardano method
void Eigen3x3_Cardano(double QQ[9], double eigenValue[3], double axis[9])
{
	double a = 1.0; // Coefficient of x^3 (always 1 for a normalized cubic equation)
	double b = -(QQ[0] + QQ[4] + QQ[8]); // From trace
	double c = (QQ[0] * QQ[4] + QQ[0] * QQ[8] + QQ[4] * QQ[8] - QQ[1] * QQ[1] - QQ[2] * QQ[2] - QQ[5] * QQ[5]);
	double d = -(QQ[0] * (QQ[4] * QQ[8] - QQ[5] * QQ[5]) - QQ[1] * (QQ[1] * QQ[8] - QQ[2] * QQ[5]) + QQ[2] * (QQ[1] * QQ[5] - QQ[2] * QQ[4]));

	// Convert to depressed cubic x^3 + px + q = 0
	double p = (3.0 * c - b * b) / 3.0;
	double q = (2.0 * b * b * b - 9.0 * b * c + 27.0 * d) / 27.0;

	double delta = (q * q / 4.0 + p * p * p / 27.0);

	std::complex<double> u, v;

	if (delta >= 0) {
		// One real root and two complex conjugate roots
		u = std::pow((-q / 2.0 + std::sqrt(delta)), 1.0 / 3.0);
		v = std::pow((-q / 2.0 - std::sqrt(delta)), 1.0 / 3.0);
	}
	else {
		// Three real roots
		double rho = std::sqrt(-p * p * p / 27.0);
		double theta = std::acos(-q / (2.0 * rho));
		u = std::pow(rho, 1.0 / 3.0) * std::complex<double>(std::cos(theta / 3.0), std::sin(theta / 3.0));
		v = std::conj(u);
	}

	// The roots of the depressed cubic
	std::complex<double> x1 = u + v;
	std::complex<double> x2 = -0.5 * (u + v) + std::complex<double>(0, std::sqrt(3.0) / 2.0) * (u - v);
	std::complex<double> x3 = -0.5 * (u + v) - std::complex<double>(0, std::sqrt(3.0) / 2.0) * (u - v);

	// The roots of the original cubic
	eigenValue[0] = (x1 - b / 3.0).real();
	eigenValue[1] = (x2 - b / 3.0).real();
	eigenValue[2] = (x3 - b / 3.0).real();

	// Eigenvectors (simplified for symmetric matrix)
	// For each eigenvalue, solve (QQ - lambda*I)v = 0
	// Since it's a 3x3, this involves solving a system of linear equations.
	// This is a simplified approach that may not be numerically stable.
	for (int i = 0; i < 3; ++i) {
		double lambda = eigenValue[i];
		// Simplified eigenvector calculation (may need more robust method)
		axis[i * 3 + 0] = 1.0;
		axis[i * 3 + 1] = (lambda - QQ[0]) == 0 ? 1.0 : QQ[1] / (lambda - QQ[0]);
		axis[i * 3 + 2] = (lambda - QQ[0]) == 0 ? 1.0 : QQ[2] / (lambda - QQ[0]);

		// Normalize eigenvector
		double norm = std::sqrt(axis[i * 3 + 0] * axis[i * 3 + 0] + axis[i * 3 + 1] * axis[i * 3 + 1] + axis[i * 3 + 2] * axis[i * 3 + 2]);
		axis[i * 3 + 0] /= norm;
		axis[i * 3 + 1] /= norm;
		axis[i * 3 + 2] /= norm;
	}
}

// On Closed-Form Formulas for the 3D Nearest Rotation Matrix Problem
// Simplified Cardano solution
void Eigen3x3_Cardano_Simplified(double QQ[9], double eigenValue[3], double axis[9])
{
	double m = (QQ[0] + QQ[4] + QQ[8]) / 3; // From trace

	double A[9];
	{
		for (int i = 0; i < 9; ++i)
			A[i] = QQ[i];
		A[0] -= m;	A[4] -= m;	A[8] -= m;
	}

	double detA = A[0] * (A[4] * A[8] - A[7] * A[5]) - A[3] * (A[1] * A[8] - A[7] * A[2]) + A[6] * (A[1] * A[5] - A[4] * A[2]);
	double q = detA / 2;

	// Frobenius norms
	double p = 0;
	{
		for (int i = 0; i < 9; ++i)
			p += A[i] * A[i];
		p /= 6;
	}

	double delta = p * p * p - q * q;

	if (delta > 0) {
		double theta = atan2(delta, q);
		double cosT = cos(theta);
		double sinT = sin(theta);

		eigenValue[0] = m + 2 * sqrt(p) * cosT;
		eigenValue[1] = m - 2 * sqrt(p) * (cosT + 3 * sinT);
		eigenValue[1] = m - 2 * sqrt(p) * (cosT - 3 * sinT);
	}
	else {
		eigenValue[0] = m;
		eigenValue[1] = m;
		eigenValue[1] = m;
	}

	// Eigenvectors (simplified for symmetric matrix)
	// For each eigenvalue, solve (QQ - lambda*I)v = 0
	// Since it's a 3x3, this involves solving a system of linear equations.
	// This is a simplified approach that may not be numerically stable.
	for (int i = 0; i < 3; ++i) {
		double lambda = eigenValue[i];
		// Simplified eigenvector calculation (may need more robust method)
		axis[i * 3 + 0] = 1.0;
		axis[i * 3 + 1] = (lambda - QQ[0]) == 0 ? 1.0 : QQ[1] / (lambda - QQ[0]);
		axis[i * 3 + 2] = (lambda - QQ[0]) == 0 ? 1.0 : QQ[2] / (lambda - QQ[0]);

		// Normalize eigenvector
		double norm = std::sqrt(axis[i * 3 + 0] * axis[i * 3 + 0] + axis[i * 3 + 1] * axis[i * 3 + 1] + axis[i * 3 + 2] * axis[i * 3 + 2]);
		axis[i * 3 + 0] /= norm;
		axis[i * 3 + 1] /= norm;
		axis[i * 3 + 2] /= norm;
	}
}

void computePrincipalAxes_3D(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices, Eigen::Vector3f &center, double eigenValue[3], double axis[9])
{
	double QQ[9];

	ComputeCovarianceMatrix_3x3(cPoints, clusterPointIndices, center, QQ);

	Eigen3x3_Jacobi(QQ, eigenValue, axis);
	
	double eigenValue0[3], axis0[9];
	Eigen3x3_Cardano(QQ, eigenValue, axis);

	double eigenValue1[3], axis1[9];
	Eigen3x3_Cardano_Simplified(QQ, eigenValue, axis);

	{
		double dEigenValue0[3];
		dEigenValue0[0] = fabs(eigenValue0[0] - eigenValue[0]);
		dEigenValue0[1] = fabs(eigenValue0[1] - eigenValue[1]);
		dEigenValue0[2] = fabs(eigenValue0[2] - eigenValue[2]);
	
		double dEigenValue1[3];
		dEigenValue1[0] = fabs(eigenValue1[0] - eigenValue[0]);
		dEigenValue1[1] = fabs(eigenValue1[1] - eigenValue[1]);
		dEigenValue1[2] = fabs(eigenValue1[2] - eigenValue[2]);
	}

}


void GetCenterAndRadius(const std::vector<cLasPOINT>& cPoints, const std::vector<int>& clusterPointIndices0, std::vector<clusterCENTER>& dbScanCenters)
{
	std::vector<int>& clusterPointIndices = const_cast<std::vector<int>&>(clusterPointIndices0);

	int n = clusterPointIndices.size();
	float minZ = 9999, maxZ = -9999;

	// robust min and max Z values
	{
		std::multiset <float> zValues;
		for (const auto& index : clusterPointIndices)
		{
			zValues.insert(cPoints[index].pt.z());
		}

		int nNoise = (int)(n * 0.001); // 0.1% noise
		if (nNoise < 1)
			nNoise = 1;

		if (zValues.size() < 3)
			nNoise = 0;

		if (!zValues.empty()) {
			auto iter0 = zValues.begin();
			auto iter1 = zValues.rbegin();
			for (int i = 0; i < nNoise; i++) {
				++iter0;
				++iter1;
			}

			minZ = *iter0;
			maxZ = *iter1;
		}

		minZ -= 0.2f;
		maxZ += 0.2f;
	}

	// 标记噪声
	{
		Eigen::Vector3f center, std;
		computeCenter_Std(cPoints, clusterPointIndices, center, std);

		int n = clusterPointIndices.size();
		for (int i = 0; i < n; i++)
		{
			int index = clusterPointIndices[i];
			double dX = cPoints[index].pt.x() - center.x();
			double dY = cPoints[index].pt.y() - center.y();

			if ( fabs(dX) > 3.75 * std.x() || fabs(dY) > 3.75 * std.y() || cPoints[index].pt.z() < minZ || cPoints[index].pt.z() > maxZ) 
			{
				// replace the point from the cluster with the last point
				clusterPointIndices[i] = clusterPointIndices[n - 1];
				--i;
				--n;
			}
		}

		// 调整大小
		clusterPointIndices.resize(n);
	}

	clusterCENTER center;

	double eigenValue[2];
	double principleAxes[2];
	float axesLen[2];

	computeAxesLen_2D(cPoints, clusterPointIndices, center.pt, eigenValue, principleAxes, axesLen);

	center.numOfPts = n;
	center.minZ = minZ;
	center.maxZ = maxZ;

	center.eigenValues[0] = eigenValue[0];
	center.eigenValues[1] = eigenValue[1];

	center.eigenVectors[0] = principleAxes[0]; 
	center.eigenVectors[1] = principleAxes[1];

	center.direction = atan2(principleAxes[1], principleAxes[0]) * 180.0 / 3.1415926;

	center.axesLen[0] = axesLen[0];
	center.axesLen[1] = axesLen[1];

	dbScanCenters.push_back(center);
}
