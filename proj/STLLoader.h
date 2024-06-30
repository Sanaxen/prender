#ifndef _STLLOADER_H
#define _STLLOADER_H

#include "ObjLoader.h"

class STLImporter
{
public:
	Obj* mesh;

	STLImporter()
	{
		mesh = new Obj();
		mesh->max[0] = -9999999.0f;
		mesh->max[1] = -9999999.0f;
		mesh->max[2] = -9999999.0f;
		mesh->min[0] = 99999999.0f;
		mesh->min[1] = 99999999.0f;
		mesh->min[2] = 99999999.0f;
		mesh->daiag = 1.0;
		mesh->nv = 0;
		mesh->nf = 0;

	}
	int Load(const char* filename, const FLOAT_DATA_T tol);
	void stor(meshDiff::Vector3d vertex[3], const FLOAT_DATA_T tol_);
};

#endif