#include "CGGToolkit.h"
using namespace CGGToolKit;

std::tuple<int,float,point> intersection3D::lineIntersectPlane(const Line3D& line,
	const Plane& plane)
{
	float t = std::numeric_limits<float>::min();
	point intersection;
	float denomiator = line.direction() DOT plane.normal();
	if (fabs(denomiator) < epsilon)
	{
		if (fabs((line.origin() DOT plane.normal()) + plane.d()) < epsilon)
		{
			t = 0.0f;
			return std::make_tuple(CGG_COINCIDE,t,intersection);
		}
		else
			return std::make_tuple(CGG_PARALLEL,t,intersection);
	}


	t = -((plane.normal() DOT line.origin()) + plane.d()) / denomiator;
	intersection = line.origin() + t * line.direction();

	return std::make_tuple(CGG_INTERSECT,t,intersection);
}
