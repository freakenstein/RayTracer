#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"

class sphere : public hitable
{
	public:
		sphere() {}
		sphere(vec3 cen, float r, material *m) : center(cen), radius(r), mat_ptr(m) {};
		virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
		vec3 center;
		float radius;
		material *mat_ptr;
};

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const
{
	float a = dot(r.direction(), r.direction());
	float b = 2 * dot(r.direction(), r.origin() - center);
	float c = dot(r.origin() - center, r.origin() - center) - radius * radius;
	float discriminant = pow(b, 2) - 4 * a*c;
	if (discriminant >= 0)
	{
		float temp = (-b - sqrt(discriminant)) / (2 * a);
		if(temp > t_min && temp < t_max)
		{
			rec.t = temp;
			rec.p = r.point_at_parameter(temp);
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mat_ptr;
			return true;
		}
		temp = (-b + sqrt(discriminant)) / (2 * a);
		if (temp > t_min && temp < t_max)
		{
			rec.t = temp;
			rec.p = r.point_at_parameter(temp);
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mat_ptr;
			return true;
		}
	}
	return false;
}

#endif