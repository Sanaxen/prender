#ifndef _COLOR_H_
#define _COLOR_H_

#include "Vector3d.h"

namespace prender {

typedef Vector3d Color;

inline float* ColotToFloat(const Color* image, const int width, const int height)
{
	float* d = new float[3 * width*height];
	int idx = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			const int id = (height - i - 1) * width + j;
			d[idx + 0] = image[id].x;
			d[idx + 1] = image[id].y;
			d[idx + 2] = image[id].z;
			idx += 3;
		}
	}
	return d;
}

inline Color* FloatToColor(const float* data, const int width, const int height)
{
	Color* image = new Color[3 * width*height];
	int idx = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			const int id = (height - i - 1) * width + j;
			image[id].x = data[idx + 0];
			image[id].y = data[idx + 1];
			image[id].z = data[idx + 2];
			idx += 3;
		}
	}
	return image;
}

};
#endif