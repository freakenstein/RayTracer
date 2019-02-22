#include <iostream>
#include <fstream>
#include <iomanip>
#include "sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include <random>
#include <functional>
#include "drand48.h"

using namespace std;

/*
vec3 color(const ray& r)
{
	vec3 unit_direction = unit_vector(r.direction());
	float t = 0.5*(unit_direction.y() + 1.0);
	return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
*/

/*
bool hit_sphere(const vec3& center, float radius, const ray& r)
{
	float dB = dot(r.direction(), r.direction());
	float dBAC = dot(r.direction(), r.origin() - center);
	float dAC = dot(r.origin() - center, r.origin() - center);
	float delta = pow(2*dBAC,2) - 4 * dB*(dAC - radius * radius); // ÅÐ±ðÊ½(discriminant) b^2 - 4ac;
	return delta > 0;
}

vec3 color(const ray& r)
{
	vec3 unit_direction = unit_vector(r.direction());
	float t = 0.5*(unit_direction.y() + 1.0);
	vec3 ret = (1.0 - t)*vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
	bool hit = hit_sphere(vec3(0.0, 0.0, -1.0),0.5,r);
	if (hit)
	{
		ret = vec3(1.0, 0.0, 0.0);
	}
	return ret;
}
*/

/*
float hit_sphere(const vec3& center, float radius, const ray& r)
{
	float a = dot(r.direction(), r.direction());
	float b = 2 * dot(r.direction(), r.origin() - center);
	float c = dot(r.origin() - center, r.origin() - center) - radius * radius;
	float delta = pow(b, 2) - 4*a*c;
	if (delta < 0)
	{
		return -1.0;
	}
	else
	{
		return (-b - sqrt(delta)) / (2 * a);
	}
}

vec3 color(const ray& r)
{
	vec3 unit_direction = unit_vector(r.direction());
	float t = 0.5*(unit_direction.y() + 1.0);
	vec3 ret = (1.0 - t)*vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);

	vec3 center = vec3(0.0, 0.0, -1.0);
	float hitT = hit_sphere(center, 0.5, r);
	if (hitT != -1.0)
	{
		vec3 normal = r.point_at_parameter(hitT) - center;
		normal.make_unit_vector();
		ret = 0.5*vec3(normal.r() + 1.0, normal.g() + 1.0, normal.b() + 1.0);
	}
	return ret;
}
*/

vec3 random_in_unit_sphere()
{
	vec3 p;
	do
	{
		p = 2.0*vec3(drand48(), drand48(), drand48()) - vec3(1,1,1);
	} while (p.squared_length() >= 1.0);
	return p;
}

vec3 color(const ray& r, hitable *world)
{
	hit_record rec;
	if (world->hit(r, 0.001, FLT_MAX, rec))
	{
		//return 0.5*vec3(rec.normal.r() + 1.0, rec.normal.g() + 1.0, rec.normal.b() + 1.0);
		vec3 target = rec.p + rec.normal + random_in_unit_sphere();
		return 0.5*color(ray(rec.p, target - rec.p), world);
	}
	else
	{
		vec3 unit_direction = unit_vector(r.direction());
		float t = 0.5*(unit_direction.y() + 1.0);
		return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
	}
}

void main()
{
	ofstream f("test.ppm");
	if (!f)
	{
		cout << "File error!\n";
		return;
	}

	int nx = 200;
	int ny = 100;
	int ns = 100;
	f<< "P3\n" << nx << " " << ny << "\n255\n";
	//std::cout << "P3\n" << nx << " " << ny << "\n255\n";
	//vec3 lower_left_corner(-2.0, -1.0, -1.0);
	//vec3 horizontal(4.0, 0.0, 0.0);
	//vec3 vertical(0.0, 2.0, 0.0);
	//vec3 origin(0.0, 0.0, 0.0);
	hitable *list[2];
	list[0] = new sphere(vec3(0, 0, -1), 0.5);
	list[1] = new sphere(vec3(0, -100.5, -1), 100);
	hitable *world = new hitable_list(list, 2);
	camera cam;
	for (int j = ny-1; j >= 0; j--)
	{
		for (int i = 0; i < nx; i++)
		{
			vec3 col(0, 0, 0);
			default_random_engine e;
			uniform_real_distribution<double> urd(0.0, double(1- DBL_MIN));
			auto dice = bind(urd, e);
			for (int s = 0; s < ns; s++)
			{
				//drand48()
				float u = float(i + dice()) / float(nx);
				float v = float(j + dice()) / float(ny);
				ray r = cam.get_ray(u, v);
				col += color(r, world);
			}
			col /= float(ns);

			//float u = float(i) / float(nx);
			//float v = float(j) / float(ny);
			//ray r = cam.get_ray(u, v);
			//col = color(r, world);

			col = vec3(pow(col[0], 0.45), pow(col[1], 0.45), pow(col[2], 0.45));
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);
			f << ir << " " << ig << " " << ib << "\n";
			//std::cout << ir << " " << ig << " " << ib << "\n";
		}
	}
	f.close();
}