#ifndef _FCL_MATERIAL_H_
#define _FCL_MATERIAL_H_
#include "include\Vec.h"


class Ray3;
typedef Vec<3, double> Vector3;

class Color
{
public:
	Color(double red, double green, double blue);
	Color add(const Color &c) const;
	Color multiply(double s) const;
	Color modulate(const Color &c) const;

	double r, g, b;
};

const Color COLOR_RED(1.0, 0.0, 0.0);
const Color COLOR_WHITE(1.0, 1.0, 1.0);
const Color COLOR_BLUE(0.0, 0.0, 1.0);
const Color COLOR_GREEN(0.0, 1.0, 0.0);

class material
{
public:
	virtual Color sample(const Ray3 &ray, const Vector3 &position, const Vector3 &normal)
	{
		return Color(1.0, 1.0, 1.0);
	};
	virtual ~material() {}
};

class phongMaterial :public material
{
public:
	phongMaterial(const Color &diff,const Color &spec, double shin, double refl, 
		const Vector3 &ld) :diffuse(diff), specular(spec), 
		shininess(shin), reflectiveness(refl), lightDir(ld)
	{
		normalize(lightDir);
	}
	virtual Color sample(const Ray3 &ray, const Vector3 &position, const Vector3 &normal);
private:
	Color diffuse, specular;
	double shininess, reflectiveness;
	Vector3 lightDir;
};

class checkerMaterial :public material
{
public:
	checkerMaterial(double sca, double ref = 0.0);
	virtual Color sample(const Ray3 &ray, const Vector3 &position, const Vector3 &normal);

	double scale, reflectiveness;
};

#endif