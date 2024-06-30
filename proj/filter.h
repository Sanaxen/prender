#ifndef _FILTER_H
#define FILTER_H

#include "color.h"
#include <vector>

namespace prender {

inline void rgb2YCoCg(double r, double g, double b, double *Y, double *Co, double *Cg) 
{
    *Y = 1/4.0f * r + 1/2.0f * g + 1/4.0f * b;
    *Co= 1/2.0f * r + 0/1.0f * g - 1/2.0f * b;
    *Cg=-1/4.0f * r + 1/2.0f * g - 1/4.0f * b;
}
inline void YCoCg2rgb(double Y, double Co, double Cg, double *r, double *g, double *b) 
{
    *r = Y + Co - Cg;
    *g = Y + 0  + Cg;
    *b = Y - Co - Cg;
}


inline void set_vector(std::vector<Color>& arr, const int x, const int y, const int ch, const int width, const int height, float *vector, int learn_radius) {
    const int size = learn_radius * 2 + 1;

    for (int oy = -learn_radius; oy <= learn_radius; ++oy) {
        for (int ox = -learn_radius; ox <= learn_radius; ++ox) {
            const int nx = ox + x;
            const int ny = oy + y;
            if (0 <= nx && nx < width && 0 <= ny && ny < height) {
                if ( ch == 0 ) vector[(oy + learn_radius) * size + (ox + learn_radius)] = arr[(ny * width + nx)].x;
                if ( ch == 1 ) vector[(oy + learn_radius) * size + (ox + learn_radius)] = arr[(ny * width + nx)].y;
                if ( ch == 2 ) vector[(oy + learn_radius) * size + (ox + learn_radius)] = arr[(ny * width + nx)].z;
            }
        }
    }
}
inline float length__(float *v0, float *v1, int size) 
{
    float sum = 0;
    for (int i = 0; i < size; ++i) {
        const float a = v0[i] - v1[i];
        sum += a * a;
    }
    return sum;
}

inline void FilterNLM(std::vector<Color>& arr0, std::vector<Color>& arr1, const int width, const int height)
{
    const int learn_radius = 3;
    const int compare_raidus = 6;
    const float sigma = 0.3f;

    // •ÏŠ·
    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            const int idx = iy * width + ix;
            const float r = arr0[idx].x;
            const float g = arr0[idx].y;
            const float b = arr0[idx].z;
            rgb2YCoCg(r, g, b, &arr0[idx].x, &arr0[idx].y, &arr0[idx].z);
        }
    }

    // NLM
    const int size = learn_radius * 2 + 1;
    for (int iy = 0; iy < height; ++iy) {
        //std::cout << "Y: " << iy << "    \r";
        #pragma omp parallel for schedule(dynamic, 1)
        for (int ix = 0; ix < width; ++ix) {

            for (int ch = 0; ch < 3; ++ch) {
                float vector0[size * size] = {0};
                set_vector(arr0, ix, iy, ch, width, height, vector0, learn_radius);
                
                const int compare_size = compare_raidus * 2 + 1;
                float weight_map[compare_size * compare_size] = {0};
                float value_map[compare_size * compare_size] = {0};

                // ’Tõ
                for (int oy = -compare_raidus; oy <= compare_raidus; ++oy) {
                    for (int ox = -compare_raidus; ox <= compare_raidus; ++ox) {
                        const int nx = ox + ix;
                        const int ny = oy + iy;
                        const int compare_idx = (oy + compare_raidus) * compare_size + (ox + compare_raidus);
                        if (0 <= nx && nx < width && 0 <= ny && ny < height) {
                            float vector1[size * size] = {0};
                            set_vector(arr0, nx, ny, ch, width, height, vector1, learn_radius);

                            // d‚ÝŒvŽZ
                            if ( ch == 0 ) value_map[compare_idx] = arr0[(ny * width + nx)].x;
                            if ( ch == 1 ) value_map[compare_idx] = arr0[(ny * width + nx)].y;
                            if ( ch == 2 ) value_map[compare_idx] = arr0[(ny * width + nx)].z;
                            weight_map[compare_idx] = length__(vector0, vector1, size * size);
                        } else {
                            weight_map[compare_idx] = -1;
                        }
                    }
                }

                // Œ‹‰ÊŒvŽZ
                float sum = 0;
                float total_weight = 0;
                for (int cy = 0; cy < compare_size; ++cy) {
                    for (int cx = 0; cx < compare_size; ++ cx) {
                        const int compare_idx = cy * compare_size + cx;
                        if (weight_map[compare_idx] < 0)
                            continue;
                        const float weight = exp(-weight_map[compare_idx] / (sigma * sigma));
                        sum += value_map[compare_idx] * weight;
                        total_weight += weight;
                    }
                }
                if (total_weight > 0)
                    sum /= total_weight;
            	if ( ch == 0 ) arr1[(iy * width + ix)].x = sum;
            	if ( ch == 1 ) arr1[(iy * width + ix)].y = sum;
            	if ( ch == 2 ) arr1[(iy * width + ix)].z = sum;
            }
            const int idx = iy * width + ix;
            const float Y = arr1[idx].x;
            const float Co= arr1[idx].y;
            const float Cg= arr1[idx].z;
            YCoCg2rgb(Y, Co, Cg, &arr1[idx].x, &arr1[idx].y, &arr1[idx].z);
        }
    }
}

class Filter
{
	std::vector<Color> tmpImg1;
	std::vector<Color> tmpImg2;
	std::vector<Color> tmpImg3;
	int w;
	int h;
	int size;
	Color* image;
public:
	Filter(Color* image_, const int width, const int height)
	{
		image = image_;
		w = width;
		h = height;
		size = w*h;
		tmpImg1.resize(size);
		tmpImg2.resize(size);
		tmpImg3.resize(size);
#pragma omp parallel for schedule(dynamic, 1)
		for (int i = 0; i < size; i++)
		{
			tmpImg1[i] = image[i];
			tmpImg3[i] = image[i];
		}
	}
	~Filter()
	{
		tmpImg1.resize(0);
		tmpImg2.resize(0);
		tmpImg3.resize(0);
		image = NULL;
	}

	void execute()
	{
		FilterNLM(tmpImg1, tmpImg2, w, h);
#pragma omp parallel for schedule(dynamic, 1)
		for (int i = 0; i < size; i++) image[i] = tmpImg2[i];
	}
	void restor()
	{
#pragma omp parallel for schedule(dynamic, 1)
		for (int i = 0; i < size; i++) image[i] = tmpImg3[i];
	}
};
	
};

#endif