#include "raytracer.h"

//the center of screen is (0,0,0) 
rayTracer::rayTracer(int w, int h, const Vector3 &e) :eye(e)
{
	width = w;
	height = h;
	lw = -0.5*width, rw = -lw;
	lh = -0.5*height, rh = -lh;
}

void
rayTracer::rayTrace(const vector<GPtr> &objects)
{
	for (int h = lh; h <= rw; h++)
	{
		for (int w = lw; w <= rw; w++)
		{
			auto dir = Vector3((double)w, (double)h, 0.0) - eye;
			normalize(dir);
			Ray3 ray(eye, dir);
			auto result = ray.intersectNearest(objects);
			if (result.intersect)
			{
				auto cr = result.object->mat->sample(ray, result.position, result.normal);
				colors.push_back(cr);
			}
			else
				colors.push_back(Color(0.0, 0.0, 0.0));
		}
	}
}

void rayTracer::renderScene()
{
	int index = 0;
	glBegin(GL_POINTS);
	for (int h = lh; h <= rh; h++)
	{
		for (int w = lw; w <= rw; w++)
		{
			glColor3d(colors[index].r, colors[index].g, colors[index].b);
			glVertex3i(w, h, 0);
			index++;
		}
	}
	glEnd();

}

void GenerateObjects(vector<GPtr> &objects)
{
	Vector3 lightDir(1.0, 1.0, 1.0);
	phongMaterial mat1(COLOR_RED, COLOR_WHITE, 36.0, 20.0, lightDir);
	phongMaterial mat2(COLOR_BLUE, COLOR_WHITE, 36.0, 20.0, lightDir);
	checkerMaterial mat3(0.05);
	objects.push_back(std::make_shared<sphere>(sphere(Vector3(-120, 100, -100), 100, mat1)));
	objects.push_back(std::make_shared<sphere>(sphere(Vector3(120, 100, -100), 100, mat2)));
	plane p(Vector3(0.0, 1.0, 0.0), 0.0, mat3);
	objects.push_back(std::make_shared<plane>(p));
}