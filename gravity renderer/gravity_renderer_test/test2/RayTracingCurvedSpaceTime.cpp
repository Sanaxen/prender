//#define _GNU_SOURCE
#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <err.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <unistd.h>
#include <errno.h>

#include "random.h"
#include "bitmap.h"

#include "WormHole.h"
#include <omp.h>
using namespace prender;

static int sizeX = 640;
static int sizeY = 480;


static BitMap* backgroundTexture1;
static BitMap* backgroundTexture2;


inline double luminance(const Rgb &color) {
	return 0.2126*color.r + 0.7152*color.g + 0.0722*color.b;
}

float GetBrightness(Rgb& color)
{
	//return luminance(color);
	float r = color.r / 255.0f;
	float g = color.g / 255.0f;
	float b = color.b / 255.0f;

	float max, min;

	max = r; min = r;

	if (g > max) max = g;
	if (b > max) max = b;

	if (g < min) min = g;
	if (b < min) min = b;

	return (max + min) / 2;
}

int Cap(int x, int max)
{
	return x > max ? max : x;
}

int CapMin(int x, int min)
{
	return x < min ? min : x;
}

Rgb AddColor(Rgb hitColor, Rgb tintColor)
{
	//return tintColor;

	//return Rgb((hitColor.r + tintColor.r)*0.5, (hitColor.g + tintColor.g)*0.5, (hitColor.b + tintColor.b)*0.5);

	float brightness = GetBrightness(tintColor);
	Rgb result = Rgb(
		(int)Cap((int)((1 - brightness) * hitColor.r) + CapMin(tintColor.r - 20, 0) * 255 / 205, 255),
		(int)Cap((int)((1 - brightness) * hitColor.g) + CapMin(tintColor.g - 20, 0) * 255 / 205, 255),
		(int)Cap((int)((1 - brightness) * hitColor.b) + CapMin(tintColor.b - 20, 0) * 255 / 205, 255)
		);
	return result;
}




//std::vector<Spherical> points;


class DiscMapping
{
private:
	double rMin;
	double rMax;
	int sizeX;
	int sizeY;

public:

	DiscMapping(double rMin_, double rMax_, int sizex_, int sizey_)
	{
		rMax = rMax_;
		rMin = rMin_;
		sizeX = sizex_;
		sizeY = sizey_;
	}


	void Map(double r, double theta, double phi, int& x, int& y)
	{
		if (r < rMin || r > rMax)
		{
			x = 0;
			y = sizeY;
		}

		x = (int)(phi / (2 * M_PI) * sizeX) % sizeX;
		if (x < 0) x = sizeX + x;
		y = (int)((r - rMin) / (rMax - rMin) * sizeY);
		if (y > sizeY - 1)
			y = sizeY - 1;

		//y = sizeY - y - 1;

	}
};

class SphericalMapping
{
	int SizeX;
	int SizeY;

public:
	SphericalMapping(int sizex_, int sizey_)
	{
		SizeX = sizex_;
		SizeY = sizey_;
	}

	void Map(double r, double theta, double phi,  int& x, int& y)
	{
		// do mapping of texture image
		double textureScale = 1.0;

		x = (int)(((phi * textureScale) / (2 * M_PI)) * this->SizeX) % this->SizeX;
		y = (int)((theta * textureScale / M_PI) * this->SizeY) % this->SizeY;

		if (x < 0) x = this->SizeX + x;
		if (y < 0) y = this->SizeY + y;

		//y = SizeY - y - 1;

	}
};
SphericalMapping* backgroundMap1;
SphericalMapping* backgroundMap2;


#define SGN(x)	( (x) < 0 ? -1 : 1)

XorShift rnd(1);
/*
r ≥ 0
0° ≤ θ ≤ 180° (π rad)
0° ≤ φ < 360° (2π rad)
*/
static void fire_ray(WormHole& wormhole, Rgb& color, double x1, double y1)
{
	Rgb hitPixel(0, 0, 0);

	double htry = 0.1, escal = 1e11, hdid = 0.0, hnext = 0.0;


	double y[WORMHOLE_ODE_N], dydx[WORMHOLE_ODE_N], yscal[WORMHOLE_ODE_N], ylaststep[WORMHOLE_ODE_N];

	int side;

	Vector3d tnv;

	Spherical ray;
	ray.th = y1 * M_PI/ (sizeY - 1);
	ray.ph = x1 * 2.0*M_PI/ (sizeX - 1);
	ray.r = 1.0;


	wormhole.initial(y, dydx, ray, tnv);

	//y[1] = ray.th;
	//y[2] = ray.ph;

	double t = 0.0;
	int loopCnt = 0;

	Vector3d tnv_pre;

	double yinit = y[0];
	while (10)
	{
		wormhole.geodesic( y, dydx);

		hnext = htry;

		t += hnext;

		hnext = wormhole.NextPosition( y, dydx, y, dydx, hnext, &tnv, &hdid);

		if ( hnext < 1.0e-14 )
		{
			hnext = 1.0e-14;
		}
		//fprintf(stderr, "l = %f\n", y[0]);

		if (wormhole.r_wh(y[0]) > 1.0e10 || fabs(y[0]) > 1.0e10)
		{
			//fprintf(stderr, "%d #l = %f r = %f\n", loopCnt, y[0], wormhole.r_wh(y[0]));
			break;
			char buf[256];
			fgets(buf, 256, stdin);
		}
		loopCnt++;
		if (loopCnt > wormhole.geodesics_max_length)
		{
			//fprintf(stderr, "t %f\n", t);
			fprintf(stderr, "l = %f r = %f\n", y[0], wormhole.r_wh(y[0]));
			break;
		}
		tnv_pre = tnv;
		htry = hnext;
	}

	double r = wormhole.r_wh(y[0]);

	Spherical s = Cartesian(tnv.x, tnv.y, tnv.z).ToSpherical();
	y[1] = s.th;
	y[2] = s.ph;

	y[1] = fmod(y[1], M_PI);
	y[2] = fmod(y[2], 2.0*M_PI);

	//背景に到達
	if (y[0] >= 0)
	{
		int xPos, yPos;

		backgroundMap1->Map(r, y[1], y[2], xPos, yPos);
		color = backgroundTexture1->cell(xPos, yPos);
		//fprintf(stderr, "+count %d\n", loopCnt);
		return;
	}

	if (y[0] < 0)
	{
		int xPos, yPos;

		backgroundMap2->Map(r, y[1], y[2], xPos, yPos);
		color = backgroundTexture2->cell(xPos, yPos);
		//fprintf(stderr, "-count %d\n", loopCnt);
		return;
	}

}


void PointsDump()
{
	//FILE* fp = fopen("points.obj", "w");

	//for (int i = 0; i < points.size(); i++)
	//{
	//	Cartesian p = points[i].ToCartesian();
	//	fprintf(fp, "v %f %f %f\n", p.x, p.y, p.z);
	//}
	//fclose(fp);
}

int main(int argc, char** argv)
{
	double r0 = 3.0;
	double theta0 = (M_PI / 180.0) * 90;
	double phi0 = (M_PI / 180.0) * 0;

	char* background1 = "InterstellarWormhole_Fig6b.jpg";
	char* background2 = "InterstellarWormhole_Fig6a.jpg";

	int geodesics_max_length = 10000;
	for ( int i = 1; i < argc; i++ )
	{
		if ( strcmp(argv[i], "-r" ) == 0 ) r0 = atof(argv[i+1]);
		if ( strcmp(argv[i], "-theta" ) == 0 ) theta0 = (M_PI / 180.0) * atof(argv[i+1]);
		if ( strcmp(argv[i], "-phi" ) == 0 ) phi0 = (M_PI / 180.0) * atof(argv[i+1]);
		if ( strcmp(argv[i], "-X" ) == 0 ) sizeX = atoi(argv[i+1]);
		if ( strcmp(argv[i], "-Y" ) == 0 ) sizeY = atoi(argv[i+1]);
		if (strcmp(argv[i], "-back1") == 0) background1 = argv[i + 1];
		if (strcmp(argv[i], "-back2") == 0) background2 = argv[i + 1];
		if (strcmp(argv[i], "-t") == 0) geodesics_max_length = atoi(argv[i + 1]);
		if (strcmp(argv[i], "-thread") == 0) omp_set_num_threads(atoi(argv[i + 1]));
	}

	Spherical cameraPos(r0, theta0, phi0);

	Cartesian cameraPosC = cameraPos.ToCartesian();
	

	WormHole wormhole(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 1.3, 0.5, 0.05, cameraPosC.ToVector(), Vector3d(0,0,1));

	wormhole.geodesics_max_length = geodesics_max_length;
	backgroundTexture1 = new BitMap;
	backgroundTexture2 = new BitMap;

	backgroundTexture1->Read(background1);
	backgroundTexture2->Read(background2);


	backgroundMap1 = new SphericalMapping(backgroundTexture1->W(), backgroundTexture1->H());
	backgroundMap2 = new SphericalMapping(backgroundTexture2->W(), backgroundTexture2->H());

	sizeX = backgroundTexture1->W();
	sizeY = backgroundTexture1->H();

	BitMap renderImage;
	renderImage.Create(sizeX, sizeY);

	double supersamples = 1.0;
	double rate = (1.0 / supersamples);

#pragma omp parallel for
	for (int j = 0; j < sizeY; j++)
	{
		for (int i = 0; i < sizeX; i++)
		{

			const int n = supersamples;
			for (int k = 0; k < n; k++)
			{
				Rgb color;
				WormHole prm = wormhole;
				fire_ray( prm, color, i + (k + 0.5)*rate, j + (k + 0.5)*rate);

				//int xPos, yPos;

				//backgroundMap1->Map(0, j*M_PI/ sizeY, i*2*M_PI/ sizeX, xPos, yPos);
				//color = backgroundTexture1->cell(xPos, yPos);

				//renderImage.cell(i, sizeY - j - 1) = color;
				renderImage.cell(i, j) = color;
			}

			renderImage.cell(i, sizeY - j - 1).r /= n;
			renderImage.cell(i, sizeY - j - 1).g /= n;
			renderImage.cell(i, sizeY - j - 1).b /= n;
		}
		printf("%.2f%%\n", 100.0*(float)j/ (float)sizeY);
		//renderImage.Write("imageTmp.bmp");
	}
	//PointsDump();

	renderImage.Write("image.bmp");

	return 0;
}

