#include "basicGeometry.h"

Ray3::Ray3(const Vector3 &orig, const Vector3 &dir) :origin(orig), direction(dir)
{
}

Vector3 Ray3::getPoint(double t) const
{
	return origin + t * direction;
}

IntersectResult Ray3::intersectNearest(const vector<shared_ptr<geometry>> &objects) const
{
	double minDistance = std::numeric_limits<double>::max();
	IntersectResult ret;
	ret.intersect = false;
	for (auto object : objects)
	{
		auto result = object->intersectTest(*this);
		if (result.intersect)
		{
			if (result.distance < minDistance)
			{
				ret = result;
				ret.object = object;
				minDistance = result.distance;
			}
		}
	}
	return ret;
}

sphere::sphere(const Vector3 &cen, double rad, const phongMaterial &m) :
geometry(std::make_shared<phongMaterial>(m)), center(cen), radius(rad)
{
	sqrRadius = radius * radius;
}

IntersectResult
sphere::intersectTest(const Ray3 &ray) const
{
	auto v = ray.origin - center;
	auto a0 = len2(v) - sqrRadius;
	auto Ddotv = ray.direction DOT v;

	IntersectResult result;
	if (Ddotv <= 0)
	{
		auto discr = Ddotv * Ddotv - a0;
		if (discr >= 0)
		{
			result.intersect = true;
			result.distance = -Ddotv - sqrt(discr);
			result.position = ray.getPoint(result.distance);
			result.normal = result.position - center;
			normalize(result.normal);
			return result;
		}
	}

	result.intersect = false;
	return result;
}

plane::plane(const Vector3 &nor, double dis, const checkerMaterial &m) :
normal(nor), distance(dis), geometry(std::make_shared<checkerMaterial>(m))
{
	position = distance * normal;
}

IntersectResult
plane::intersectTest(const Ray3 &ray) const
{
	IntersectResult result;
	auto a = ray.direction DOT normal;
	if (a > 0 || fabs(a) < 1e-6)
	{
		result.intersect = false;
		return result;
	}

	auto b = normal DOT (ray.origin - position);
	result.intersect = true;
	result.distance = -b / a;
	result.position = ray.getPoint(result.distance);
	result.normal = normal;
	return result;
}
