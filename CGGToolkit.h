#ifndef _CGGTOOLKIT_H_
#define _CGGTOOLKIT_H_
/*
This class collected methods in <Geometric Tools for Computer Graphics>,
all of methods are 'public' and 'static'.

Implemented by Chilin Fu.
*/
#include "Vec.h"
#include "GeometryFigure.h"
#include <vector>
#include <tuple>
using std::vector;

namespace CGGToolKit
{
	const float epsilon = 1e-3f;
	const int CGG_PARALLEL = 0;
	const int CGG_INTERSECT = 1;
	const int CGG_COINCIDE = 2;
	//3 dimension intersection
	class intersection3D
	{
	public:
		static std::tuple<int,float,point> lineIntersectPlane(const Line3D& line, 
			const Plane& plane);


	};
}


#endif