#ifndef _BASICGEOMETRY_H_
#define _BASICGEOMETRY_H_
#include "include\Vec.h"
#include "material.h"
#include <memory>
#include <vector>
using std::vector;
using std::shared_ptr;

class geometry;
typedef Vec<3, double> Vector3;

struct IntersectResult
{
	bool intersect;
	Vector3 position;
	Vector3 normal;
	double distance;
	shared_ptr<geometry> object;
};

class Ray3
{
public:
	Ray3(const Vector3 &orig, const Vector3 &dir);
	Vector3 getPoint(double t) const;
	IntersectResult intersectNearest(const vector<shared_ptr<geometry>> &objects);

	Vector3 origin, direction;
};

class geometry
{
public:
	geometry(shared_ptr<material> m) :mat(m){}
	virtual IntersectResult intersectTest(const Ray3 &ray) const 
	{ 
		return IntersectResult(); 
	}
	virtual ~geometry() {}

	shared_ptr<material> mat;
};

class sphere :public geometry
{
public:
	sphere(const Vector3 &cen, double rad, const phongMaterial &m);
	virtual IntersectResult intersectTest(const Ray3 &ray) const;


	Vector3 center;
	double radius;
	double sqrRadius;
};

class plane :public geometry
{
public:
	plane(const Vector3 &nor, double dis, const checkerMaterial &m);
	virtual IntersectResult intersectTest(const Ray3 &ray) const;

	Vector3 normal, position;
	double distance;
};

#endif