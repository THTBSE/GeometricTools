#include "GeometryFigure.h"

using namespace CGGToolKit;


Plane::Plane(const point& pointA, const point& pointB, const point& pointC)
{
	vec edge1 = pointB - pointA;
	vec edge2 = pointC - pointA;
	nor = edge1 CROSS edge2;
	normalize(nor);
	vec v(pointA);
	distance = -(v DOT nor);
}