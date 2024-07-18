#ifndef _WORMHOLE_
#define _WORMHOLE_

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <string>
#include "constant.h"
#include "vector3d.h"
using namespace prender;

#define BACK_GLOUND2		-997

#define WORMHOLE_ODE_N	5
//inline double Clamp(double val, double low, double high) {
//	if (val < low) return low;
//	else if (val > high) return high;
//	else return val;
//}

//class Spherical;
//
//class Cartesian
//{
//public:
//	double x;
//	double y;
//	double z;
//	double r;
//
//	double xx;
//	double yy;
//	double zz;
//
//	Cartesian(double x_, double y_, double z_) :x(x_), y(y_), z(z_)
//	{
//		xx = yy = zz = 0.0;
//		r = x*x + y*y + z*z;
//		if (r > 1.0e-10)
//		{
//			r = sqrt(r);
//			xx = x / r;
//			yy = y / r;
//			zz = z / r;
//		}
//	}
//
//	inline double SphericalTheta() const
//	{
//		return acos(Clamp(zz, -1.0, 1.0));
//	}
//	inline double SphericalPhi()  const
//	{
//		double p = atan2(yy, xx);
//		return (p < 0.0) ? p + 2.0*M_PI : p;
//	}
//
//	inline Spherical ToSpherical() const;
//
//	inline Vector3d ToVector() const
//	{
//		return Vector3d(x, y, z);
//	}
//	inline Vector3d ToNormalVector() const
//	{
//		return Vector3d(xx, yy, zz);
//	}
//
//};
//
//class Spherical
//{
//public:
//	double r;
//	double th;
//	double ph;
//
//	Spherical() {}
//	Spherical(double r_, double th_, double ph_) :r(r_), th(th_), ph(ph_)
//	{
//	}
//
//	inline Cartesian ToCartesian() const
//	{
//		return Cartesian(
//			r * cos(ph) * sin(th),
//			r * sin(ph) * sin(th),
//			r * cos(th)
//			);
//	}
//	inline Vector3d ToVector() const
//	{
//		Cartesian p = ToCartesian();
//		return Vector3d(p.x, p.y, p.z);
//	}
//	inline Vector3d ToNormalVector() const
//	{
//		Cartesian p = ToCartesian();
//		return Vector3d(p.xx, p.yy, p.zz);
//	}
//
//
//};
//inline Spherical Cartesian::ToSpherical() const
//{
//	return Spherical(sqrt(x*x + y*y + z*z), SphericalTheta(), SphericalPhi());
//}


class WormHole;

inline double rkqs(WormHole& wormhole, double *y, double *dydx, double htry, double escal, double *yscal, double *hdid);
inline void rkstep(WormHole& wormhole, double *y, double *dydx, double h, double *yout, double *yerr);

inline double UnitStep(const double x)
{
	if (x < 0.0) return 0.0;
	return 1.0;
}


class WormHole
{
	
public:
	double l0, theta0, phi0;

	double Rho;	// the radius of the wormhole
	double a;	// the length of the wormhole
	double M;	//a parameter describing the curvature, described in the paper as the "gentleness of the transition from the wormhole's cylindrical interior to its asymptotically flat exterior"
	double W;

	double Nx;
	double Ny;
	double Nz;

	double n_l;
	double n_ph;
	double n_th;


	double traceDir;

	double geodesics_max_length;

	Vector3d camera_pos;
	Vector3d camera_dir;
	Vector3d camera_up;
	Vector3d screen_x;
	Vector3d screen_y;

	Vector3d position0;
	Vector3d position1;

	int sign;
	int usr_set_sign_ = 0;

	std::string background_texture[2];
	double background_texture_map_coef[2][4];
	double background_texture_coef[2];

	double r_wh(double l)
	{
		double x = 2.0*(fabs(l) - a) / (PS_PI*M);

		if (fabs(l) > a)
		{
			return  Rho + M*(x*atan(x) - 0.5*log(1.0 + x*x));
		}
		return Rho;
	}


	double drdl(double l)
	{
		double s = 1;
		if (l < 0) s = -1;
		if (l == 0.0) s = 0.0;

		if (fabs(l) > a)
		{
			return (2.0*atan((2.0*(-a + fabs(l))) / (M*PS_PI))*s) / PS_PI;
		}
		return 0;
	}

	double rtol(double r, double linit)
	{
		double x = 2.0*(fabs(linit) - a) / (PS_PI*M);

		int stat = -1;
		for (int i = 0; i < 30; i++)
		{
			double xi = x;
			x = x - (Rho + M*(x*atan(x) - 0.5*log(1.0 + x*x)) - r) / (M*atan(x));

			if (fabs(x - xi) < 1.0e-6)
			{
				stat = 0;
				break;
			}
			if (!isfinite(x))
			{
				break;
			}
		}
		double l = x*PS_PI*M*0.5 + a;

		if (stat != 0)
		{
			x = a;
			for (int i = 0; i < 30; i++)
			{
				double xi = x;
				x = x - (Rho + M * (x * atan(x) - 0.5 * log(1.0 + x * x)) - r) / (M * atan(x));

				if (fabs(x - xi) < 1.0e-6)
				{
					stat = 0;
					break;
				}
				if (!isfinite(x))
				{
					break;
				}
			}
		}
		if (stat != 0)
		{
			//fprintf(stderr, "error l=%f\n", (linit < 0)? -fabs(l): fabs(l));
		}
		return (linit < 0) ? -fabs(l) : fabs(l);
	}

	inline double zfunction(const double l, const double ll) const
	{
		double p = (UnitStep(fabs(ll) - a)*(2.0 * atan((2.0 * (-a + fabs(ll))) / (M*M_PI))*((l < 0) ? -1.0 : 1.0)) / M_PI);

		return sqrt(1 - p*p);
	}
	//Simpson 法積分
	inline double NIntegrate(const double a, const double b, const double l, const int divnum) const
	{
		double integral; // 積分結果
		double h; // 積分範囲を n 個に分割したときの幅
		double x, dA;
		int i;


		h = (b - a) / (2.0*divnum);
		x = a;
		integral = zfunction(l, x);

		for (i = 1; i<divnum; i++) {
			dA = 4.0*zfunction(l, x + h) + 2.0*zfunction(l, x + 2.0*h);
			integral += dA;
			x += 2.0*h;
		}

		integral += (4.0*zfunction(l, x + h) + zfunction(l, x + 2.0*h));
		integral *= h / 3.0;

		return(integral);
	}

	//NIntegrate[Sqrt[1 - (UnitStep[Abs[ll] -a] (2 ArcTan[(2 (-a + Abs[ll]))/(M π)] Sign[l])/π)^2], {ll, 0, l}]
	inline double distance_through(const double l) const
	{
		return NIntegrate(0.0, l, l, 10);
	}


	WormHole(double x0_, double y0_, double z0_, double x1_, double y1_, double z1_, double rho, double a_, double W_, Vector3d& cameraPos, Vector3d& cameraDir)
	{
		background_texture[0] = "";
		background_texture[1] = "";
		background_texture_map_coef[0][0] = 1.0;
		background_texture_map_coef[0][1] = 0.0;
		background_texture_map_coef[1][0] = 1.0;
		background_texture_map_coef[1][1] = 0.0;

		camera_pos = cameraPos;
		camera_dir = normalize(cameraDir);
		traceDir = 1.0;
		position0 = Vector3d(x0_, y0_, z0_);
		position1 = Vector3d(x1_, y1_, z1_);

		Rho = rho;
		W = W_;
		a = fabs(a_);
		M = W/1.42953;

		sign = 1;
		Cartesian pos0(cameraPos.x, cameraPos.y, cameraPos.z);
		if ( cameraPos.length() > 1000000)
		{
			sign = -1.0;
		}

		pos0.x -= position0.x;
		pos0.y -= position0.y;
		pos0.z -= position0.z;
		Spherical poss0 = pos0.ToSpherical();

		l0 = rtol(poss0.r, sign*a/2);
		theta0 = poss0.th;
		phi0 = poss0.ph;

		//if (fabs(l0) < a)
		//{
		//	sign = -1;
		//	l0 = sign * (a);
		//}

		//if (!isfinite(l0))
		//{
		//	sign = -1;
		//	l0 = sign * (a);
		//}

		//if (usr_set_sign_)
		//{
		//	sign = usr_set_sign_;
		//	l0 = sign * fabs(l0);
		//}
	}

	int Setup( double x0, double y0, double z0, Spherical& ray, double sgn)
	{
		Cartesian pos(x0, y0, z0);

		if ( sgn >= 0 )
		{
			pos.x -= position0.x;
			pos.y -= position0.y;
			pos.z -= position0.z;
		}else
		{
			pos.x -= position1.x;
			pos.y -= position1.y;
			pos.z -= position1.z;
		}
		Spherical poss = pos.ToSpherical();

		Mat.LoadIdentity();
		invMat.LoadIdentity();
#if 10
		Vector3d X = normalize(camera_pos - position0);
		if ( sgn < 0 )
		{
			X = normalize(camera_pos - position1);
		}

		Vector3d Y = screen_y;
		Vector3d Z = screen_x;
		//OrthonormalBasis(X, X, Y, Z);

		Y = normalize(cross(Vector3d(0,0,1), X));
		if ( Y.length() < 1.0e-16 )
		{
			Y = normalize(cross(Vector3d(0,1,0), X));
		}
		Z = normalize(cross(X, Y));

		Matrix4D mat(
			X.x, Y.x, Z.x, 0,
			X.y, Y.y, Z.y, 0,
			X.z, Y.z, Z.z, 0,
			0, 0, 0, 1.0);
		mat = mat.Transpose(mat);
		
		Matrix4D rt;
		rt.LoadIdentity();

		Mat = mat;

		double det;
		bool stat = InvertMatrix_(Mat, invMat, det);

		Spherical possS = Cartesian(Mat*poss.ToVector()).ToSpherical();
		theta0 = possS.th;
		phi0 = possS.ph;
		l0 = rtol( possS.r, sgn*a / 2);


		//l0 = rtol(poss.r, sgn*a / 2);
		//theta0 = poss.th;
		//phi0 = poss.ph;

#else
		l0 = rtol(poss.r, sgn*a / 2);
		theta0 = poss.th;
		phi0 = poss.ph;
#endif
		return 1;
	}

	Matrix4D Mat;
	Matrix4D invMat;

	/* Initial Conditions for Ray */
	void initial( double *y0, double* dydx, Spherical& ray,  const Vector3d& tnv)
	{
#if 0
		Nx = sin(ray.th)*cos(ray.ph);		//(A9a)
		Ny = sin(ray.th)*sin(ray.ph);
		Nz = cos(ray.th);

		n_l = -Nx;							//(A9b)
		n_ph = -Ny;
		n_th = Nz;
#else

		Spherical rayS = Cartesian(Mat*ray.ToCartesian().ToVector()).ToSpherical();
		
		const double sint = sin(rayS.th);
		Nx = sint*cos(rayS.ph);		//(A9a)
		Ny = sint*sin(rayS.ph);
		Nz = cos(rayS.th);

		n_l = -Nx;							//(A9b)
		n_ph = -Ny;
		n_th = Nz;
#endif

		y0[0] = l0;
		y0[1] = theta0;
		y0[2] = phi0;
		y0[3] = n_l;					//(A9c)
		y0[4] = r_wh(l0)*n_th;			//(A9c)
	}

	inline Vector3d WorldAxisSystem(double xx, double yy, double zz, double sgn) const
	{
		if ( sgn >= 0 )
		{
			return invMat*Vector3d(xx, yy, zz) + position0;
		}
		return invMat*Vector3d(xx, yy, zz) + position1;
	}

	inline Vector3d WorldAxisSystem(Spherical& pos, double sgn) const
	{
		Cartesian q = pos.ToCartesian();
		return WorldAxisSystem(q.x, q.y, q.z, sgn);
	}
	inline Vector3d WorldAxisSystem(Cartesian& pos, double sgn) const
	{
		return WorldAxisSystem(pos.x, pos.y, pos.z, sgn);
	}

	inline Vector3d Coordinate_transformation(double xx, double yy, double zz, double sgn) const
	{
		return Coordinate_transformation(Vector3d(xx, yy, zz), sgn);
	}
	inline Vector3d Coordinate_transformation(Vector3d& p, double sgn) const
	{
		if ( sgn >= 0 )
		{
			return Mat*(p - position0);
		}
		return Mat*(p - position1);
	}

	void geodesic(double* y, double* dydx)
	{

		const double r = r_wh(y[0]);
		const double r2 = r*r;

		double sint = sin(y[1]);

		if (fabs(sint) < 1.0e-8)
		{
			if (sint > 0.0) sint = 1.0e-8;
			else sint = -1.0e-8;
		}
		const double sint2 = sint*sint;

		const double r_l0 = r_wh(l0);
		const double b = r_l0*sin(theta0)*n_ph;											//(A9d)
		const double B2 = r_l0*r_l0*(n_th*n_th + n_ph*n_ph);			//(A9d)

		const double invr2 = 1.0 / r2;
		dydx[0] = y[3];										//(A7a)
		dydx[1] = y[4] * invr2;								//(A7b)
		dydx[2] = b*invr2 / sint2;							//(A7c)
		dydx[3] = B2*drdl(y[0])*invr2 /r;					//(A7d)
		dydx[4] = (b*b*invr2)* (cos(y[1]) / (sint*sint2));			//(A7e)
	}

	//微分方程式の次の位置と接ベクトルを計算する
	double NextPosition(double* y, double* dydx,  double* y_next, double* dydx_next, double htry, Vector3d* tnv, double *hdid)
	{
		double escal = 1e11, hnext = 0.0;

		double y1[WORMHOLE_ODE_N], dydx1[WORMHOLE_ODE_N], yscal[WORMHOLE_ODE_N], ylaststep[WORMHOLE_ODE_N], ytemp[WORMHOLE_ODE_N], yerr[WORMHOLE_ODE_N];

		memcpy(y1, y, WORMHOLE_ODE_N*sizeof(double));
		memcpy(dydx1, dydx, WORMHOLE_ODE_N*sizeof(double));

		for (int i = 0; i < WORMHOLE_ODE_N; i++)
		{
			yscal[i] = fabs(y1[i]) + fabs(dydx1[i] * htry) + 1.0e-3;
		}

		double hh = 0.0;
		double* hd = hdid;
		if (hd == 0) hd = &hh;

		hnext = rkqs( *this, y1, dydx1, htry, escal, yscal, hd);
		if (tnv)
		{

			Spherical org( r_wh(y[0]), y[1], y[2]);
			Cartesian c_org = org.ToCartesian();

			Spherical p(r_wh(y1[0]), y1[1], y1[2]);

			Cartesian q = p.ToCartesian();

			tnv->x = q.x - c_org.x;
			tnv->y = q.y - c_org.y;
			tnv->z = q.z - c_org.z;

			double len = tnv->length();
			if ( len > 1.0e-16)
			{
				tnv->x /= len;
				tnv->y /= len;
				tnv->z /= len;
			}
		}
		if (y_next && dydx_next)
		{
			memcpy(y_next, y1, WORMHOLE_ODE_N*sizeof(double));
			memcpy(dydx_next, dydx1, WORMHOLE_ODE_N*sizeof(double));
		}

		return hnext;
	}
};

inline void rkstep(WormHole& wormhole, double *y, double *dydx, double h, double *yout, double *yerr)
{

	double ak[WORMHOLE_ODE_N];

	double ytemp1[WORMHOLE_ODE_N], ytemp2[WORMHOLE_ODE_N], ytemp3[WORMHOLE_ODE_N], ytemp4[WORMHOLE_ODE_N], ytemp5[WORMHOLE_ODE_N];

	h *= -1.0*wormhole.traceDir;
	for (int i = 0; i < WORMHOLE_ODE_N; i++)
	{
		double hdydx = h * dydx[i];
		double yi = y[i];
		ytemp1[i] = yi + 0.2 * hdydx;
		ytemp2[i] = yi + (3.0 / 40.0) * hdydx;
		ytemp3[i] = yi + 0.3 * hdydx;
		ytemp4[i] = yi - (11.0 / 54.0) * hdydx;
		ytemp5[i] = yi + (1631.0 / 55296.0) * hdydx;
		yout[i] = yi + (37.0 / 378.0) * hdydx;
		yerr[i] = ((37.0 / 378.0) - (2825.0 / 27648.0)) * hdydx;
	}

	wormhole.geodesic( ytemp1, ak);

	for (int i = 0; i < WORMHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp2[i] += (9.0 / 40.0) * yt;
		ytemp3[i] -= 0.9 * yt;
		ytemp4[i] += 2.5 * yt;
		ytemp5[i] += (175.0 / 512.0) * yt;
	}

	wormhole.geodesic( ytemp2, ak);

	for (int i = 0; i < WORMHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp3[i] += 1.2 * yt;
		ytemp4[i] -= (70.0 / 27.0) * yt;
		ytemp5[i] += (575.0 / 13824.0) * yt;
		yout[i] += (250.0 / 621.0) * yt;
		yerr[i] += ((250.0 / 621.0) - (18575.0 / 48384.0)) * yt;
	}

	wormhole.geodesic( ytemp3, ak);

	for (int i = 0; i < WORMHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp4[i] += (35.0 / 27.0) * yt;
		ytemp5[i] += (44275.0 / 110592.0) * yt;
		yout[i] += (125.0 / 594.0) * yt;
		yerr[i] += ((125.0 / 594.0) - (13525.0 / 55296.0)) * yt;
	}

	wormhole.geodesic(ytemp4, ak);

	for (int i = 0; i < WORMHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp5[i] += (253.0 / 4096.0) * yt;
		yerr[i] -= (277.0 / 14336.0) * yt;
	}

	wormhole.geodesic( ytemp5, ak);

	for (int i = 0; i < WORMHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		yout[i] += (512.0 / 1771.0) * yt;
		yerr[i] += ((512.0 / 1771.0) - 0.25) * yt;
	}
}

//適応ルンゲ・クッタ積分アルゴリズム。
//実際に実装された方法のバリエーションは、Cash-Karp.として知られています。
//ステップhtry
inline double rkqs(WormHole& wormhole, double *y, double *dydx, double htry, double escal, double *yscal, double *hdid)
{
	int i;

	double hnext;

	double errmax, h = htry, htemp;
	double yerr[WORMHOLE_ODE_N], ytemp[WORMHOLE_ODE_N];

	int loopCnt = 0;
	while (loopCnt < 2000)
	{
		// Run a single step of integration using step h 
		rkstep(wormhole, y, dydx, h, ytemp, yerr);

		//誤差評価
		errmax = 0.0;
		for (i = 0; i < WORMHOLE_ODE_N; i++)
		{
			double temp = fabs(yerr[i] / yscal[i]);
			if (temp > errmax) errmax = temp;
		}

		//誤差Maxが1.0以下ならOK
		errmax *= escal;
		if (errmax <= 1.0) break;

		//ステップをもう少し小さく設定する
		htemp = 0.9 * h / sqrt(sqrt(errmax));

		h *= 0.1;

		if (h >= 0.0)
		{
			if (htemp > h) h = htemp;
		}
		else
		{
			if (htemp < h) h = htemp;
		}
		loopCnt++;
	}
	if (loopCnt == 2000)
	{
		printf("--------------------\n");
	}

	//次回のステップを調整
	if (errmax > 1.89e-4)
	{
		hnext = 0.9 * h * pow(errmax, -0.2);
	}
	else
	{
		hnext = 5.0 * h;
	}

	*hdid = h;

	memcpy(y, ytemp, WORMHOLE_ODE_N * sizeof(double));

	return hnext;
}


#endif
