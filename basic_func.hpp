#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <chrono>
#include <algorithm>
#include <iterator>
#include <map>
#include <eigen3/Eigen/Core>
//#include <Eigen/Core>

/* Define structs -------------------------- */
struct CoordInt { int x, y, z; };

struct CoordFloat{ float x, y, z; };

//struct CoordId {float x,y,z; int t,id;};

struct CoordId {int x,y,z,t,id;};

struct CoordPixel
{
	float x,y,z;
	int greyval, pxnum, label;
};

struct CoordIntGrey
{
	std::vector<CoordInt> raypts;
	int greyval;
};

struct PointCloud
{
	struct Point
	{
		size_t id {0};
		float x, y, z, vx, vy, vz, ax, ay, az;
	};
	Point spt;
	std::vector<Point> pts;
	inline unsigned int kdtree_get_point_count() const { return pts.size(); }
	inline float kdtree_get_pt(const unsigned int idx, const unsigned int dim) const {
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }
};

/* Cross product ---------------------- */
CoordFloat cross_pdt( CoordFloat& a, CoordFloat& b )
{
	CoordFloat c;
	c.x = (a.y * b.z) - (a.z * b.y);
	c.y = (a.z * b.x) - (a.x * b.z);
	c.z = (a.x * b.y) - (a.y * b.x);
	return c;
}

/* Square of norm of vector ------------ */
float normSq( CoordFloat a ){
	return a.x * a.x + a.y * a.y + a.z * a.z;
}

PointCloud readPtCloud(const std::string filename) {
	PointCloud cloud;

	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error: unable to open file " << filename << '\n';
		return cloud;
	}
  
	while (true) {
		file >> cloud.spt.x >> cloud.spt.y >> cloud.spt.z;
		if (file.eof()) break;
		cloud.pts.emplace_back(cloud.spt);
	}

	file.close();
	return cloud;
}

template<typename T>
std::vector<T> readCoord(const std::string filename) {
	std::vector<T> coordinates;
	std::ifstream file(filename);

	if (!file.is_open()) {
		std::cerr << "Error: unable to open file " << filename << '\n';
		return coordinates;
	}
  
	while (true) {
		T coordinate;
		file >> coordinate.x >> coordinate.y >> coordinate.z;
		if (file.eof()) break;
		coordinates.emplace_back(coordinate);
	}

	file.close();
	return coordinates;
}

template<typename T>
std::vector<T> readfile(const std::string filename) {

	std::vector<T> my_array;
	std::ifstream file(filename);
	float num{};

	while (file >> num)
		my_array.emplace_back(num);

	file.close();
	return my_array;
}

template<typename T>
void writeToFile(const std::string filename, std::vector<T>& data) {
	std::ofstream outdata;
	outdata.open(filename, std::ios::app);
	for (T& points : data) {
		outdata << points << '\t';
	}
	outdata.close();
}

template<typename T>
void writeCoord(std::string filename, std::vector<T>& data)
{
	std::ofstream outdata;
	outdata.open(filename, std::ios::app);
	for (T& points : data) {
		outdata << points.x << '\t' << points.y << '\t' << points.z << '\n';
	}
	outdata.close();
}

/* Delete all contents of a file ------------------------------ */
void clearfile(const std::string& filename)
{
	std::ofstream outdata;
	outdata.open(filename);
	outdata.close();
}

void writePtCloud(const std::string filename, int t, PointCloud& data)
{
	std::ofstream outdata;
	outdata.open(filename, std::ios::app);
	for (size_t i=0; i<data.kdtree_get_point_count(); i++) {
		outdata << data.pts[i].x << '\t' << data.pts[i].y << '\t' << data.pts[i].z << '\t' << data.pts[i].id << '\t' << t << '\n';
	}
	outdata.close();
}

void writeCoordId(const std::string filename, std::vector<CoordId> data)
{
	std::ofstream outdata;
	outdata.open(filename, std::ios::app);
	for (size_t i=0; i<data.size(); i++) {
		outdata << data[i].x << '\t' << data[i].y << '\t' << data[i].z << '\t' << data[i].id << '\t' << data[i].t << '\n';
	}
	outdata.close();
}

/* Convert degree to radians -------------- */
std::vector<CoordFloat> DegToRad(std::vector<CoordFloat>& angle_deg)
{
	CoordFloat radang;
	std::vector<CoordFloat> radanglist;

	const float conv {3.14159/180.0};
	for (CoordFloat& ang : angle_deg) {
		radang.x = ang.x*conv; radang.y = ang.y*conv; radang.z = ang.z*conv;
		radanglist.emplace_back(radang);
	}

	return radanglist;
}

/* Calculate 3D rotation matrix, inverts if inv=1 --------- */
Eigen::Matrix3f rotmat(float& theta, float& phi, float& gamma, const bool& inv)
{
	Eigen::Matrix3f Rx; Eigen::Matrix3f Ry; Eigen::Matrix3f Rz;
	Rx << 1,0,0,0,cos(theta),-sin(theta),0,sin(theta),cos(theta);
	Ry << cos(phi),0,sin(phi),0,1,0,-sin(phi),0,cos(phi);
	Rz << cos(gamma),-sin(gamma),0,sin(gamma),cos(gamma),0,0,0,1;

	if (inv) {return (Rx*Ry*Rz).transpose();}
	else {return Rx*Ry*Rz;}
}

/* measure total run time of program -------- */
class Timer 
{
private:
	using Clock = std::chrono::steady_clock;
	using Second = std::chrono::duration<float, std::ratio<1>>;
	std::chrono::time_point<Clock> m_beg{ Clock::now() };
public:
	void reset() { m_beg = Clock::now(); }
	float elapsed() const{
		return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
	}
};
