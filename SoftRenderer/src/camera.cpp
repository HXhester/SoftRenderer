#include "camera.h"

#define FRUSTUM_PLANE_COUNT 6

Camera::Camera() {}

int Camera::pointInFrustum(Vec3f& p) {
	for (int i=0; i < FRUSTUM_PLANE_COUNT; i++) {
		if (plane[i].distance(p) < 0)
			return Outside;
	}

	return Inside;
}

void Camera::setCameraPram(float fov, float ratio, float near, float far) {

}