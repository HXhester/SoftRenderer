#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "geometry.h"

// Camera with frustum
class Camera {
private:
	enum { Top, Bottom, Left, Right, Near, Far };

public:
	static enum { Inside, Intersect, Outside };
	Vec3f pos;
	Plane plane[6];
	Vec3f ntl, ntr, nbr;
	Vec3f fbr, ftr, ftl;

	float fov, ratio, near, far;
	Camera();
	int pointInFrustum(Vec3f& p);
	void setCameraPram(float fov, float ratio, float near, float far);
};

#endif
