#ifndef RAYH
#define RAYH

#include "vec3.h"

class ray
{
public:
	ray() {}
	ray(const vec3& a, const vec3& b) { A = a; B = b;}
	vec3 origin() const { return A; } //起始点
	vec3 direction() const { return B; } //起始点作为原点的方向
	vec3 point_at_parameter(float t) const { return A + t * B; } //so that

	vec3 A;
	vec3 B;
};

#endif