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

#include "KerrBlackHole.h"

static int sizeX = 640;
static int sizeY = 480;


static BitMap* backgroundTexture;
static BitMap* discTexture;

static double  cameraTilt = 0.0;

#define RAY_TEST	10

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
DiscMapping *discMap;
SphericalMapping* backgroundMap;




XorShift rnd(1);
/*
r ≥ 0
0° ≤ θ ≤ 180° (π rad)
0° ≤ φ < 360° (2π rad)
*/
static void fire_ray(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, Rgb& color, double x1, double y1)
{
	Rgb pixel(-1, 0, 0);
	Rgb hitPixel(0, 0, 0);

	double htry = 0.5, escal = 1e11, hdid = 0.0, hnext = 0.0;

	double range = 1.0 / (sizeX - 1); //0.0025 * Rdisk / (sizeX - 1.0);

	double y[KERRBLACKHOLE_ODE_N], dydx[KERRBLACKHOLE_ODE_N], yscal[KERRBLACKHOLE_ODE_N], ylaststep[KERRBLACKHOLE_ODE_N];

	int side;


#if 0
	kerrBlackHole.initial(y, dydx, (x1 - (sizeX + 1.0) / 2) * range, (y1 - (sizeY + 1.0) / 2) * range);
#else

	double tiltSin = sin((cameraTilt / 180) * M_PI);
	double tiltCos = cos((cameraTilt / 180) * M_PI);

	double xRot = x1 - (sizeX + 1.0) / 2;
	double yRot = y1 - (sizeY + 1.0) / 2;

	Vector3d v(0, 0, 0);
	kerrBlackHole.initial(prm, y, dydx,
		(xRot * tiltCos - yRot * tiltSin) * range,
		(yRot * tiltCos + xRot * tiltSin) * range,
		v);

#if RAY_TEST
	//計算されたRayベクトルと一致する微分方程式初期条件を求める
	kerrBlackHole.RayToOrdinaryDifferentialEquationInitial(prm, v, y, dydx);
#endif

#endif


	//points.push_back(Spherical(y[0], y[1], y[2]));

	int loopCnt = 0;
	while (1)
	{
		memcpy(ylaststep, y, KERRBLACKHOLE_ODE_N * sizeof(double));

		kerrBlackHole.geodesic(prm, y, dydx);


		if (y[1] > M_PI / 2) side = 1;
		else if (y[1] < M_PI / 2) side = -1;
		else side = 0;

		//if (kerrBlackHole.M < 1.0e-10) side = 0;

#if RAY_TEST
		Vector3d tnv;
		hnext = kerrBlackHole.NextPosition(prm, y, dydx, y, dydx, 0, &hdid);

#else
		for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
		{
			yscal[i] = fabs(y[i]) + fabs(dydx[i] * htry) + 1.0e-3;
		}
		hnext = rkqs(prm, kerrBlackHole, y, dydx, htry, escal, yscal, &hdid);
#endif
		//points.push_back(Spherical(y[0], y[1], y[2]));

		if ((y[1] - M_PI / 2)*side < 0)
		{
			memcpy(y, ylaststep, KERRBLACKHOLE_ODE_N * sizeof(double));

			//交点を探す
			binarysearch(prm, kerrBlackHole, y, dydx, hdid);

			//降着円盤にヒット
			/* Did we hit the disk? */
			//Ray hits accretion disc?
			if ((y[0] <= kerrBlackHole.Rdisk) && (y[0] >= kerrBlackHole.Rmstable))
			{
				// y[0] - radial position
				// y[2] - phi (horizontal) angular position

				unsigned p1 = 4.0 * (y[0] - kerrBlackHole.Rmstable) / (kerrBlackHole.Rdisk - kerrBlackHole.Rmstable);
				unsigned p2 = floor(y[2] * 6.0 / M_PI);

#if 0
				if ((p1 ^ p2) & 1)
				{
					rgb[0] = 255;
					rgb[1] = 128;
					rgb[2] = 128;
				}
				else
				{
					rgb[0] = 255;
					rgb[1] = 0;
					rgb[2] = 0;
				}
#else
				int xPos, yPos;
				if ( discTexture->GetImage())
				{
					discMap->Map(y[0], y[1], y[2],  xPos,  yPos);
					if (pixel.r >= 0)
					{
						pixel = AddColor(discTexture->cell(xPos, yPos), pixel);
					}
					else
					{
						pixel = discTexture->cell(xPos, yPos);
					}
				}
				// don't return yet, just remember the color to 'tint' the texture later 
#endif
				//return;
			}
		}

		/* Inside the hole, or escaped to infinity */

		//事象の地平面に到達
		if (y[0] < kerrBlackHole.Rhor)
		{
			hitPixel = Rgb(0, 0, 0);

			if (pixel.r >= 0)
			{
				color = AddColor(hitPixel, pixel);
			}
			else
			{
				color = hitPixel;
			}
			return;
		}

		//背景に到達
		if (y[0] > kerrBlackHole.r0)
		{
#if 0
			color = Rgb(128, 128, 128);
#else
			int xPos, yPos;

			backgroundMap->Map(y[0], y[1], y[2], xPos, yPos);
			hitPixel = backgroundTexture->cell(xPos, yPos);
			if (pixel.r >= 0)
			{
				color = AddColor(hitPixel, pixel);
			}
			else
			{
				color = hitPixel;
			}
			return;
#endif
		}
		loopCnt++;
		if (loopCnt > 10000)
		{
			color = Rgb(1, 0, 0);
			return;
		}

		htry = hnext;
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
	double r0 = 50.0;
	double theta0 = (M_PI / 180.0) * 85;
	double phi0 = (M_PI / 180.0) * 0;


	char* background = "InterstellarWormhole_Fig6b.jpg";
	char* disc = "adisk_skewed.png";

	for ( int i = 1; i < argc; i++ )
	{
		if ( strcmp(argv[i], "-r" ) == 0 ) r0 = atof(argv[i+1]);
		if ( strcmp(argv[i], "-theta" ) == 0 ) theta0 = (M_PI / 180.0) * atof(argv[i+1]);
		if ( strcmp(argv[i], "-phi" ) == 0 ) phi0 = (M_PI / 180.0) * atof(argv[i+1]);
		if ( strcmp(argv[i], "-X" ) == 0 ) sizeX = atoi(argv[i+1]);
		if ( strcmp(argv[i], "-Y" ) == 0 ) sizeY = atoi(argv[i+1]);
		if ( strcmp(argv[i], "-back" ) == 0 ) background = argv[i+1];
		if ( strcmp(argv[i], "-disc" ) == 0 ) disc = argv[i+1];
		if ( strcmp(argv[i], "--disc" ) == 0 ) disc = 0;
	}

	Spherical cameraPos(r0, theta0, phi0);

	Cartesian cameraPosC = cameraPos.ToCartesian();
	
	double M = 1.0;
	double a = -0.5;
	double Rdisk = 20.0;

	KerrBlackHole kerrBlackHole(0, 0, 50, Rdisk, M, a, cameraPos);

	backgroundTexture = new BitMap;
	discTexture = new BitMap;

	//backgroundTexture->Read("bgedit.jpg");
	//discTexture->Read("adisk_skewed2.png");

	backgroundTexture->Read(background);
	if ( disc ) discTexture->Read(disc);

	//backgroundTexture->Read("InterstellarWormhole_Fig10.jpg");
	//discTexture->Read("adisk_skewed.png");


	if ( disc ) discMap = new DiscMapping(kerrBlackHole.Rmstable, Rdisk, discTexture->W(), discTexture->H());
	backgroundMap = new SphericalMapping(backgroundTexture->W(), backgroundTexture->H());

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
				KerrBlackHoleWrkParam prm;
				fire_ray( prm, kerrBlackHole, color, i + (k + 0.5)*rate, j + (k + 0.5)*rate);
				renderImage.cell(i, sizeY - j - 1) = color;
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

