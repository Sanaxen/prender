#ifndef _RAY_H_
#define _RAY_H_

#include "Vector3d.h"

namespace prender {

struct Ray {
	Vector3d org, dir;
	double doppler_factor; // ドップラー効果の係数

	Ray(const Vector3d &org, const Vector3d &dir) : org(org), dir(dir) {}
};

};

#endif
