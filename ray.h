#ifndef _RAY_H_
#define _RAY_H_

#include "Vector3d.h"

namespace prender {

struct Ray {
	Vector3d org, dir;
	Ray(const Vector3d &org, const Vector3d &dir) : org(org), dir(dir) {}
};

};

#endif
