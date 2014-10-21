#ifndef _GEOMETRYPROC_H_
#define _GEOMETRYPROC_H_
#include "Vec.h"


namespace CGGToolKit
{
	//The instance of Line3D is immutable,it can not be modified after construction.
	class Line3D
	{
	public:
		Line3D():d(vec(0.0f,0.0f,1.0f)){}
		Line3D(const point& linePoint, const vec& dirOrPoint, bool isPoint) :p(linePoint)
		{
			if (isPoint)
				d = dirOrPoint - linePoint;
			else
				d = dirOrPoint;
		}

		const point& origin() const { return p; }
		const vec& direction() const { return d; }
	protected:
		point p;
		vec   d;
	};

	/*
	a plane is represented with formula 'ax + by + cz + d = 0'.
	vector (a,b,c) is the unit normal vector of the plane.
	d is the minimum distance from base point(0,0,0) to the plane.
	(x,y,z) is the point in the plane.
	so the dot product of n(a,b,c) and v(x,y,z) is '-d'.
	Here the Plane is also immutable.
	*/
	class Plane
	{
	public:
		//default constructor build a xoy plane 
		Plane() :nor(vec(0.0f, 0.0f, 1.0f)), distance(0.0f){}
		//PointA,B,C is counterclockwise
		Plane(const point& pointA, const point& pointB, const point& pointC);
		float a() const { return nor[0]; }
		float b() const { return nor[1]; }
		float c() const { return nor[2]; }
		float d() const { return distance; }
		const vec& normal() const { return nor; }

	private:
		vec nor;
		float distance;
	};


}




#endif 