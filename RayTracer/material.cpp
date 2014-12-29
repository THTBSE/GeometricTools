#include "material.h"
#include "basicGeometry.h"

Color::Color(double red, double green, double blue) :r(red), g(green), b(blue)
{

}

Color Color::add(const Color &c) const
{
	return Color(r + c.r, g + c.g, b + c.b);
}

Color Color::multiply(double s) const
{
	return Color(r*s, g*s, b*s);
}

Color Color::modulate(const Color &c) const
{
	return Color(r*c.r, g*c.g, b*c.b);
}

Color
phongMaterial::sample(const Ray3 &ray, const Vector3 &position,
const Vector3 &normal)
{
	auto NdotL = normal DOT lightDir;
	auto H = lightDir - ray.direction;
	normalize(H);
	auto NdotH = normal DOT H;
	auto diffuseTerm = diffuse.multiply(std::max(NdotL, 0.0));
	auto specularTerm = specular.multiply(pow(std::max(NdotH, 0.0), shininess));
	return COLOR_WHITE.modulate(diffuseTerm.add(specularTerm));
}

checkerMaterial::checkerMaterial(double sca, double ref) :material(ref), scale(sca)
{
}

Color
checkerMaterial::sample(const Ray3 &ray, const Vector3 &position, const Vector3 &normal)
{
	return abs((int)(position[0] * scale) + (int)(position[2] * scale)) % 2 < 1 ? Color(0.0, 0.0, 0.0) :
		Color(1.0, 1.0, 1.0);
}