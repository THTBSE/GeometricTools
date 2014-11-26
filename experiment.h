#ifndef _EXPERIMENT_H_
#define _EXPERIMENT_H_
#include "../Mathematics/GteIntrLine3Cone3.h"
#include "../Mathematics/GteIntrSegment3Cone3.h"
#include "../Mathematics/GteIntrRay3Cone3.h"
#include "../Mathematics/GteMatrix4x4.h"
#include "../Mathematics/GteMatrix3x3.h"
#include "../Mathematics/GteRotation.h"
#include "../Mathematics/GteIntrSegment3Triangle3.h"
#include "../Mathematics/GteIntrLine2Triangle2.h"
#include "../trimesh/include/TriMesh.h"
#include "../Mathematics/GteIntrLine3Triangle3.h"
#include <array>
#include <vector>
#include <map>

class IntrPair
{
public:
	IntrPair(double _1st = 0.0, double _2nd = 0.0, 
		size_t index = std::numeric_limits<size_t>::max()) :TriIndex(index) 
	{
		first = std::min(_1st, _2nd);
		second = std::max(_1st, _2nd);
	}

	double first, second;
	size_t TriIndex;

	friend bool operator< (const IntrPair& lhs, const IntrPair& rhs)
	{
		return (lhs.first + lhs.second) < (rhs.first + rhs.second);
	}
};

template<int N>
class CompOneDim
{
public:
	CompOneDim(size_t i) :dim(i){}

	bool operator() (const gte::Vector<N,double>& lhs, const gte::Vector<N,double>& rhs)
	{
		return lhs[dim] < rhs[dim];
	}
private:
	size_t dim;
};

class experiment
{
	typedef gte::Vector4<double> gteVec4;
	typedef gte::Vector3<double> gteVec3;
	//crdSystem represent a coordinate system
	typedef std::array<gteVec3, 4> crdSystem;
public:
	//transform from one local coordinate system to another.
	//orig and tar are the origin local coordinate system and target system respectively.
	//orig[3] is the base point.
	//orig[0],[1],[2] are x,y,z axis respectively,and they must be unit-length.
	//atom is the point need to transform.
	void interstellar(std::array<gteVec4, 4>& orig, std::array<gteVec4, 4>& tar, gteVec3& atom);

	//build a coordinate system from a vertex and a vector;
	//the vector then be treated as z-axis.
	crdSystem getCrdSystem(const gte::Vector3<double>& Vec, const gte::Vector3<double>& Vertex);

	//make a coordinate to homogeneous.
	std::array<gteVec4, 4>
		makeHomogeneous(const std::array<gteVec3, 4>& crdSys);

	void draw()
	{
		drawPoints();
		drawVectors();
		drawCrdSystem();
	}
	//draw coordinate system by OpenGL
	void drawCrdSystem();

	//draw Vectors
	void drawVectors();

	//experiment data
	std::vector<std::pair<gteVec3, gteVec3>> Vectors;
	std::vector<crdSystem> crdSystems;

	/*           Cone Intersection with Cone              */
	std::vector<gte::Vector3<double>>
		IntrConeToCone(const gte::Cone3<double>& c1, const gte::Cone3<double>& c2);
	std::vector<gte::Vector3<double>>
		IntrConeToConeInfinite(const gte::Cone3<double>& c1, const gte::Cone3<double>& c2);
	//get cone bottom disk points;
	std::vector<gte::Vector3<double>>
		GetConeDisk(const gte::Cone3<double>& c1);

	void IntrConeExperiment();
	//experiment data
	std::vector<gteVec3> points;
	void drawPoints();

	/*        Cone Intersection with Triangle Mesh        */
	std::vector<gte::Vector3<double>>
		IntrConeToTriangle(const gte::Cone3<double>& cone, std::vector<gte::Triangle3<double>>& TList);
	std::vector<gte::Triangle3<double>>
		GetTriangleList(TriMesh* mesh);

	//experiment data
	std::vector<gte::Triangle3<double>> TriList;

	/*  Get Sampling points  */

	//Getting support points from specified triangles list.
	std::vector<gte::Vector3<double>>
		GetSupportPoint(const std::vector<int>& IndexList, 
		const std::vector<gte::Triangle3<double>>& TriList, double resolution = 0.25);

	/*  Compute overhang triangles ,return their index  */
	std::vector<int>
		GetOverhangTriangle(const std::vector<gte::Triangle3<double>>& TList, double alphcC = GTE_C_QUARTER_PI);
};



#endif