#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

#include "ArHosekSkyModel.h"

#include "bitmap.h"

typedef struct RawImage_ {
  double r, g, b;
} RawImage;

void conv(RawImage* pixel, unsigned char* out)  {
  float d;
  int e;
  d = pixel->r > pixel->g ? pixel->r : pixel->g;
  if (pixel->b > d) d = pixel->b;
  
  if (d <= 1e-32) {
    out[0] = out[1] = out[2] = 0;
    out[3] = 0;
    return;
  }
  
  float m = frexp(d, &e); // d = m * 2^e
  d = m * 256.0 / d;
  
  out[0] = pixel->r * d;
  out[1] = pixel->g * d;
  out[2] = pixel->b * d;
  out[3] = (e + 128);
}


void saveHDRfile(char *filename, RawImage *imgdata, int width, int height) {
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) {
    printf("Error: %s", filename);
    return;
  }

  unsigned char ret = 0x0a;
  // Header
  fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
  fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
  fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
  fprintf(fp, "EXPOSURE=          1.0000000000000%c%c", ret, ret);

  // Data
  fprintf(fp, "-Y %d +X %d%c", height, width, ret);
  for (int i = height - 1; i >= 0; i --) {

	  unsigned char *buf = new unsigned char[width * 4];
    for (int j = 0; j < width; j ++) {
      conv(&imgdata[i * width + j], &buf[j * 4]);
    }

    fprintf(fp, "%c%c", 0x02, 0x02);
    fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);

    for (int j = 0; j < 4; j ++) {
      int cursor = 0;
      for (;;) {
	int w = width - cursor;
	if (w >= 128)
	  w = 127;
	fprintf(fp, "%c", w);
	for (int idx = cursor; idx < cursor + w; idx ++) 
	    fprintf(fp, "%c", buf[idx * 4 + j]); 
	cursor += w;
	if (cursor >= width)
	  break;
      }
    }
	delete [] buf;
  }
	
  fclose(fp);
}

// tekito converter
void xyz2rgb(double* xyz, double* rgb) {
  double X = xyz[0];
  double Y = xyz[1];
  double Z = xyz[2];
  rgb[0] = 3.2410*X - 1.5374*Y - 0.4986*Z;
  rgb[1] = -0.9692*X + 1.8760*Y + 0.0416*Z;
  rgb[2] = 0.0556*X - 0.2040*Y + 1.5070*Z;
}

inline double clamp(double x){
	if (x < 0.0)
		return 0.0;
	if (x > 1.0)
		return 1.0;
	return x;
}

inline int to_int(double x){
	return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

// 
// skydome filename width height turbidity albedo solarElevation
//
int main(int argc, char **argv) {  
	ArHosekSkyModelState* skymodel_state[3];
  if (argc < 7)
    return 0;

  int width = atof(argv[2]);
  int height = atof(argv[3]);
  double turbidity = atof(argv[4]);
  double albedo = atof(argv[5]);
  double solarElevation = atof(argv[6])*3.14159 / 180.0;
  RawImage *img = (RawImage*)malloc(sizeof(RawImage) * width * height);

  if ( turbidity < 1 ) turbidity = 1.0;
  if ( turbidity > 10 ) turbidity = 10.0;
  if ( albedo < 0.0 ) albedo = 0.0;
  if ( albedo > 1.0 ) albedo = 1.0;

#if 10
 skymodel_state[0] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo, solarElevation);
 skymodel_state[1] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo, solarElevation);
 skymodel_state[2] = arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo, solarElevation);
#else
 skymodel_state[0] = arhosek_xyz_skymodelstate_alloc_init(turbidity, albedo, solarElevation);
 skymodel_state[1] = arhosek_xyz_skymodelstate_alloc_init(turbidity, albedo, solarElevation);
 skymodel_state[2] = arhosek_xyz_skymodelstate_alloc_init(turbidity, albedo, solarElevation);
#endif

	const double PI = 3.14159265358979323846264338327950288;

#if 0
  for (int y = 0; y < height; y ++) {
    for (int x = 0; x < width; x ++) {

      int cx = width / 2;
      int cy = height / 2;
      int rx = (x - cx);
      int ry = (y - cy);

      double nr = sqrt( (double)(rx*rx + ry*ry)) / (width / 2.0);
	  double u = rx / (width / 2.0);
	  double v = ry / (width / 2.0);
	  double theta = atan2(v, u);
	  double ph = PI*sqrt(u*u+v*v);

	  if (nr <= 1.0) 
	  {
			double xyz[3];
			double rgb[3];

			for (int i = 0; i < 3; i ++) 
			{
#if 10
				rgb[i] = arhosek_tristim_skymodel_radiance(skymodel_state[i], theta, ph, i);
#else
				xyz[i] = arhosek_tristim_skymodel_radiance(skymodel_state[i], theta, ph, i);
#endif
			}
#if 10
			/* empty */
#else
			xyz2rgb(xyz, rgb);
#endif
			img[y * width + x].r = rgb[0] *0.02;
			img[y * width + x].g = rgb[1] *0.02;
			img[y * width + x].b = rgb[2] *0.02;
	   } else {
			img[y * width + x].r = 0;
			img[y * width + x].g = 0;
			img[y * width + x].b = 0;
      }
    }
  }
#else

  for (int y = 0; y < height; y ++) {
    for (int x = 0; x < width; x ++) {

      double cx = width / 2.0;
      double cy = height / 2.0;
      double rx = (x - cx);
      double ry = (y - cy);

      double nr = sqrt( (double)(rx*rx + ry*ry)) / ((width+width*0.01) / 2.0);
      double th = nr * 0.5 * PI;
	  double ph = atan2((double)rx, (double)ry);
      
	  double gamma = acos(cos(solarElevation) * sin(th) * sin(ph) + sin(solarElevation) * cos(th));
      double theta = th;

	  if (nr < 1.0) 
	  {
			double xyz[3];
			double rgb[3];

			for (int i = 0; i < 3; i ++) 
			{
#if 10
				rgb[i] = arhosek_tristim_skymodel_radiance(skymodel_state[i], theta, gamma, i);
#else
				xyz[i] = arhosek_tristim_skymodel_radiance(skymodel_state[i], theta, gamma, i);
#endif
			}
#if 10
			/* empty */
#else
			xyz2rgb(xyz, rgb);
#endif
			img[y * width + x].r = rgb[0] *0.02;
			img[y * width + x].g = rgb[1] *0.02;
			img[y * width + x].b = rgb[2] *0.02;
	   } else {
			img[y * width + x].r = 0;
			img[y * width + x].g = 0;
			img[y * width + x].b = 0;
      }
    }
  }
#endif

  arhosekskymodelstate_free(skymodel_state[0]);
  arhosekskymodelstate_free(skymodel_state[1]);
  arhosekskymodelstate_free(skymodel_state[2]);
  


  BitMap bmp;

  bmp.Create(width, height);
  for (int y = 0; y < height; y++) {
	  for (int x = 0; x < width; x++) {
		  bmp.cell(y, x).r = to_int(img[y * width + x].r);
		  bmp.cell(y, x).g = to_int(img[y * width + x].g);
		  bmp.cell(y, x).b = to_int(img[y * width + x].b);
	  }
  }
  bmp.Write("ccc.bmp");

  // save
  saveHDRfile(argv[1], img, width, height);

 free(img);

  return 0;
}
