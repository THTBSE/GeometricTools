#ifndef _RAYTRACER_H_
#define _RAYTRACER_H_
#include "include\GLFW\glfw3.h"
#include "basicGeometry.h"
#include "material.h"
#include <vector>
using namespace std;
typedef shared_ptr<geometry> GPtr;

void GenerateObjects(vector<GPtr> &objects);

class rayTracer
{
public:
	rayTracer(int w, int h, const Vector3 &e);
	void renderScene();
	void simpleExp(int width, int height);


	void rayTrace(const vector<GPtr> &objects);
	
	int width, height, lw, rw, lh, rh;
	Vector3 eye;
	vector<Color> colors;
};




#endif