#define _CRT_SECURE_NO_WARNINGS
#include <ctype.h>
#include "vector3d.h"
#include "STLLoader.h"


int STLImporter::Load( const char* filename, const FLOAT_DATA_T tol)
{

	FILE* fp = fopen(filename, "r");

	if (fp == NULL)
	{
		return -1;
	}

	char buf[256];

	meshDiff::Vector3d normal;
	meshDiff::Vector3d vertex[3];

	bool stl_start = false;
	while (fgets(buf, 256, fp) != NULL)
	{
		char* line = buf;
		while (isspace(*line)) line++;

		if (strncmp(line, "solid", 5) == 0)
		{
			stl_start = true;
			continue;
		}
		if (strncmp(line, "endsolid", 8) == 0)
		{
			if ( stl_start ) break;
			else return -1;
		}
		if (strncmp(line, "facet normal", 12) == 0)
		{
#if FLOAT_FLAG==1
			sscanf(line, "facet normal %f %f %f", &normal.x, &normal.y, &normal.z);
#else
			sscanf(line, "facet normal %lf %lf %lf", &normal.x, &normal.y, &normal.z);
#endif
			continue;
		}
		if (strncmp(line, "endloop", 7) == 0)
		{
			continue;
		}
		if (strncmp(line, "endfacet", 8) == 0)
		{
			stor(vertex, tol);
			continue;
		}

		if (strncmp(line, "outer loop", 10) == 0)
		{
			for (int i = 0; i < 3; i++)
			{
				fgets(buf, 256, fp);
				char* line = buf;
				while (isspace(*line)) line++;
				if (strncmp(line, "vertex", 6) == 0)
				{
#if FLOAT_FLAG==1
					sscanf(line, "vertex %f %f %f", &(vertex[i].x), &(vertex[i].y), &(vertex[i].z));
#else
					sscanf(line, "vertex %lf %lf %lf", &(vertex[i].x), &(vertex[i].y), &(vertex[i].z));
#endif
					continue;
				}
			}
			continue;
		}
	}
	fclose(fp);
	if (!stl_start) return -3;

	mesh->fnorm.resize(mesh->nf);
	mesh->norm.resize(mesh->nv);
	mesh->norm_id.resize(mesh->nf);
	mesh->CalcNormalVector();

	mesh->dump2();

	return 0;
}


void STLImporter::stor(meshDiff::Vector3d vertex[3], const FLOAT_DATA_T tol_)
{
	const FLOAT_DATA_T tol = tol_*tol_;

	const float v[3][3] = {
		{ vertex[0].x, vertex[0].y, vertex[0].z },
		{ vertex[1].x, vertex[1].y, vertex[1].z },
		{ vertex[2].x, vertex[2].y, vertex[2].z },
	};

	Vector3df vert[3];
	vert[0].g[0] = v[0][0];
	vert[0].g[1] = v[0][1];
	vert[0].g[2] = v[0][2];
	vert[1].g[0] = v[1][0];
	vert[1].g[1] = v[1][1];
	vert[1].g[2] = v[1][2];
	vert[2].g[0] = v[2][0];
	vert[2].g[1] = v[2][1];
	vert[2].g[2] = v[2][2];

	idx vector_idx;

	vector_idx.g[0] = UINT_MAX;
	vector_idx.g[1] = UINT_MAX;
	vector_idx.g[2] = UINT_MAX;

	const int n = mesh->nv;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < n; i++)
	{
		if (vector_idx.g[0] != UINT_MAX && vector_idx.g[1] != UINT_MAX && vector_idx.g[2] != UINT_MAX)
		{
#ifdef _OPENMP
			continue;
#else
			break;
#endif
		}

		const meshDiff::Vector3d v = meshDiff::Vector3d(mesh->vert[i].g);
		const FLOAT_DATA_T len1 = dot((v - vertex[0]), (v - vertex[0]));
		const FLOAT_DATA_T len2 = dot((v - vertex[1]), (v - vertex[1]));
		const FLOAT_DATA_T len3 = dot((v - vertex[2]), (v - vertex[2]));
#pragma omp critical
		{
			if (len1 < tol) vector_idx.g[0] = i;
			if (len2 < tol) vector_idx.g[1] = i;
			if (len3 < tol) vector_idx.g[2] = i;
		}
	}

	if (vector_idx.g[0] == UINT_MAX)
	{
		mesh->vert.push_back(vert[0]);
		vector_idx.g[0] = (int)mesh->vert.size() - 1;
	}
	if (vector_idx.g[1] == UINT_MAX)
	{
		mesh->vert.push_back(vert[1]);
		vector_idx.g[1] = (int)mesh->vert.size() - 1;
	}
	if (vector_idx.g[2] == UINT_MAX)
	{
		mesh->vert.push_back(vert[2]);
		vector_idx.g[2] = (int)mesh->vert.size() - 1;
	}
	mesh->face.push_back(vector_idx);

	mesh->nv = (int)mesh->vert.size();
	mesh->nf = (int)mesh->face.size();
}
