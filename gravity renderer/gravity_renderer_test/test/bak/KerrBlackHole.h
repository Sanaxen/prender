#ifndef _KERRBLACKHOLE_
#define _KERRBLACKHOLE_

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constant.h"
#include "vector3d.h"

#define KERRBLACKHOLE_ODE_N	6
//inline double Clamp(double val, double low, double high) {
//	if (val < low) return low;
//	else if (val > high) return high;
//	else return val;
//}

class Spherical;

class Cartesian
{
public:
	double x;
	double y;
	double z;
	double r;

	double xx;
	double yy;
	double zz;

	Cartesian(double x_, double y_, double z_) :x(x_), y(y_), z(z_)
	{
		xx = yy = zz = 0.0;
		r = x*x + y*y + z*z;
		if (r > 1.0e-10)
		{
			r = sqrt(r);
			xx = x / r;
			yy = y / r;
			zz = z / r;
		}
	}

	inline double SphericalTheta() const
	{
		return acos(Clamp(zz, -1.0, 1.0));
	}
	inline double SphericalPhi()  const 
	{
		double p = atan2(yy, xx);
		return (p < 0.0) ? p + 2.0*M_PI : p;
	}

	inline Spherical ToSpherical() const;

	inline Vector3d ToVector() const
	{
		return Vector3d(x, y, z);
	}
	inline Vector3d ToNormalVector() const
	{
		return Vector3d(xx, yy, zz);
	}

};

class Spherical
{
public:
	double r;
	double th;
	double ph;

	Spherical(double r_, double th_, double ph_) :r(r_), th(th_), ph(ph_)
	{
	}

	inline Cartesian ToCartesian() const
	{
		return Cartesian(
			r * cos(ph) * sin(th),
			r * sin(ph) * sin(th),
			r * cos(th)
			);
	}
	inline Vector3d ToVector() const
	{
		Cartesian p = ToCartesian();
		return Vector3d(p.x, p.y, p.z);
	}
	inline Vector3d ToNormalVector() const
	{
		Cartesian p = ToCartesian();
		return Vector3d(p.xx, p.yy, p.zz);
	}


};
inline Spherical Cartesian::ToSpherical() const
{
	return Spherical(sqrt(x*x + y*y + z*z), SphericalTheta(), SphericalPhi());
}

inline void sincos(double t, double* sint, double* cost)
{
	*sint = sin(t);
	*cost = cos(t);
}

class KerrBlackHoleWrkParam
{
public:
	double L;		//angular momentum in phi direction
	double kappa;	// Carter's constant's element

};

class KerrBlackHole;

static double rkqs(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, double *y, double *dydx, double htry, double escal, double *yscal, double *hdid);
static void rkstep(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, double *y, double *dydx, double h, double *yout, double *yerr);

class KerrBlackHole
{
public:
	// Initial conditions :
	// Ray starting location in Boyer-Lindquist coordinates
	double r0, theta0, phi0;

	//double L;		//angular momentum in phi direction
	//double kappa;	// Carter's constant's element

	double a;	// angular momentum in theta direction (1.0を超えると裸の特異点が出てしまう）
	double a2;	// a-squared

				// Dimensions of the "accretion disc"
	double Rhor;
	double Rmstable;
	double Rdisk;

	double M;

	double x;
	double y;
	double z;
	KerrBlackHole(double x_, double y_, double z_, double rDisc, double m, double a_, Spherical& cameraPos)
	{
		x = x_;
		y = y_;
		z = z_;

		//a = -0.5; // angular momentum in theta direction (hard-coded)

		a = a_;
		M = m;
		Cartesian pos = cameraPos.ToCartesian();

		pos.x -= x;
		pos.y -= y;
		pos.z -= z;
		Spherical poss = pos.ToSpherical();

		r0 = poss.r;
		theta0 = poss.th;
		phi0 = poss.ph;

		a2 = a * a;

		//Rhor = 1.0 + sqrt(1.0 - a2) + 1e-5;
		Rhor = M + sqrt(M*M - a2) + 1e-5;	//event horizon
		Rdisk = rDisc*M;
		Rmstable = inner_orbit();
	}

	double cbrt(const double x) const
	{
		return pow(x, 1.0 / 3.0);
	}
	double inner_orbit(void) const
	{
		double m2 = M*M;
		double z1, z2;
		if (M > 1.0e-10)
		{
			z1 = 1 + cbrt(1 - a2/m2)*(cbrt(1 + a/m2) + cbrt(1 - a/m2));
			z2 = sqrt(3 * a2/m2 + z1*z1);
		}else
		{ 
			z1 = 1;
			z2 = sqrt(3  + z1*z1);
		}
		return M*(3 + z2 - sqrt((3 - z1)*(3 + z1 + 2 * z2)));
	}

	/* Initial Conditions for Ray */
	// 接ベクトルを返す
	void initial(KerrBlackHoleWrkParam& prm, double *y0, double *ydot0, double x, double y, Vector3d& tnv)
	{
		y0[0] = r0;
		y0[1] = theta0;
		y0[2] = phi0;
		y0[3] = 0;
		y0[4] = 0;
		y0[5] = 0;

		double sintheta, costheta;
		sincos(theta0, &sintheta, &costheta);
		double cos2 = costheta*costheta;
		double sin2 = sintheta*sintheta;

		double rdot0 = cos(y)*cos(x);
		double thetadot0 = sin(y) / r0;

		double r2 = r0 * r0;
		double sigma = r2 + a2*cos2;			//Σ = r^2 + a^2 cos^2(θ)
		double delta = r2 - 2.0 * M*r0 + a2;	//Δ = r^2 - 2Mr + a^2
		double s1 = sigma - 2.0 * r0;

		y0[4] = rdot0*sigma / delta;
		y0[5] = thetadot0*sigma;

		ydot0[0] = rdot0;
		ydot0[1] = thetadot0;
		ydot0[2] = cos(y)*sin(x) / (r0*sin(theta0));

		double phidot0 = ydot0[2];
		double energy2 = s1*(rdot0*rdot0 / delta + thetadot0*thetadot0)
			+ delta*sin2*phidot0*phidot0;

		double energy = sqrt(energy2);

		/* Rescale */
		y0[4] = y0[4] / energy;
		y0[5] = y0[5] / energy;

		/* Angular Momentum with E = 1 */
		prm.L = ((sigma*delta*phidot0 - 2.0*a*r0*energy)*sin2 / s1) / energy;

		prm.kappa = y0[5] * y0[5] + a2*sin2 + prm.L*prm.L / sin2;

		/* Hack - make sure everything is normalized correctly by a call to geodesic */
		geodesic(prm, y0, ydot0);


		if (1) {
			double y1[KERRBLACKHOLE_ODE_N], dydx1[KERRBLACKHOLE_ODE_N];

			memcpy(y1, y0, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx1, ydot0, KERRBLACKHOLE_ODE_N*sizeof(double));

			//接ベクトルを求める
			NextPosition(prm, y1, dydx1, 0, 0, &tnv, 0);
		}
	}
	
	inline Vector3d WorldAxisSystem(double xx, double yy, double zz) const
	{
		return Vector3d(xx + x, yy + y, zz + z);
	}

	inline Vector3d WorldAxisSystem(Spherical& pos) const
	{
		Cartesian q = pos.ToCartesian();
		return WorldAxisSystem(q.x, q.y, q.z);
	}

	//微分方程式の次の位置と接ベクトルを計算する
	double NextPosition(KerrBlackHoleWrkParam& prm, double* y, double* dydx, double* y_next, double* dydx_next, Vector3d* tnv, double *hdid)
	{
		double htry = 0.5, escal = 1e11, hnext = 0.0;

		double y1[KERRBLACKHOLE_ODE_N], dydx1[KERRBLACKHOLE_ODE_N], yscal[KERRBLACKHOLE_ODE_N], ylaststep[KERRBLACKHOLE_ODE_N], ytemp[KERRBLACKHOLE_ODE_N], yerr[KERRBLACKHOLE_ODE_N];

		memcpy(y1, y, KERRBLACKHOLE_ODE_N*sizeof(double));
		memcpy(dydx1, dydx, KERRBLACKHOLE_ODE_N*sizeof(double));

		for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
		{
			yscal[i] = fabs(y1[i]) + fabs(dydx1[i] * htry) + 1.0e-3;
		}

		double hh = 0.0;
		double* hd = hdid;
		if (hd == 0) hd = &hh;

		hnext = rkqs(prm, *this, y1, dydx1, htry, escal, yscal, hd);

		if (tnv)
		{
			Spherical org(y[0], y[1], y[2]);
			Cartesian c_org = org.ToCartesian();

			Spherical p(y1[0], y1[1], y1[2]);
			Cartesian q = p.ToCartesian();

			double xx = q.x - c_org.x;
			double yy = q.y - c_org.y;
			double zz = q.z - c_org.z;

			double r = xx*xx + yy*yy + zz*zz;
			if (r > 1.0e-10)
			{
				r = sqrt(r);

				tnv->x = xx / r;
				tnv->y = yy / r;
				tnv->z = zz / r;
			}
		}
		if (y_next && dydx_next)
		{
			memcpy(y_next, y1, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx_next, dydx1, KERRBLACKHOLE_ODE_N*sizeof(double));
		}

		return hnext;
	}


	//Rayベクトルと一致する微分方程式初期条件を求める
	void RayToOrdinaryDifferentialEquationInitial(KerrBlackHoleWrkParam& prm, Vector3d& dir, double* y, double* dydx)
	{
		//printf("%f %f %f\n", dir.xx, dir.yy, dir.zz);

		double y0[KERRBLACKHOLE_ODE_N], ydot0[KERRBLACKHOLE_ODE_N];


		Spherical org(r0, theta0, phi0);
		Cartesian c_org = org.ToCartesian();

		Vector3d v(0, 0, 0);
		bool lockup = false;
		double xx, yy;


		xx = 0.0;
		yy = 0.0;

		double delta = 3.0;
		double step = 2.0*delta / 5.0;
		double dmin = 999999999.0;
		double tol = 1.0e-8;
		KerrBlackHoleWrkParam prm0;
		double y2[KERRBLACKHOLE_ODE_N];
		double dydx2[KERRBLACKHOLE_ODE_N];

		for (int k = 0; k < 40; k++)
		{
			double xxx, yyy;
			bool lockup_min = false;

			lockup = false;
			for (double x = xx - delta; x <= xx + delta; x += step)
			{
				for (double y = yy - delta; y <= yy + delta; y += step)
				{

					KerrBlackHoleWrkParam prm1;
					initial(prm1, y0, ydot0, x, y, v);

					double l = dot(v, dir);

					l = fabs(1.0 - l);
					if (l < dmin)
					{
						dmin = l;
						xxx = x;
						yyy = y;
						prm0 = prm1;
						memcpy(y2, y0, KERRBLACKHOLE_ODE_N*sizeof(double));
						memcpy(dydx2, ydot0, KERRBLACKHOLE_ODE_N*sizeof(double));
						lockup_min = true;
						//printf("%d LoockUP!! !!%f\n", k, dmin);
					}
				}
			}

			if (dmin < tol)
			{
				//printf("@@@@@@%d LoockUP!! !!%f\n", k, dmin);
				lockup = true;
				prm = prm0;
				memcpy(y, y2, KERRBLACKHOLE_ODE_N*sizeof(double));
				memcpy(dydx, dydx2, KERRBLACKHOLE_ODE_N*sizeof(double));
				break;
			}

			if (lockup_min)
			{
				delta = 2.0*step;
				step = 2.0*delta / 10.0;
				prm = prm0;
				memcpy(y, y2, KERRBLACKHOLE_ODE_N*sizeof(double));
				memcpy(dydx, dydx2, KERRBLACKHOLE_ODE_N*sizeof(double));
				xx = xxx;
				yy = yyy;
			}
			else
			{
				step *= 0.5;
			}
			//printf("dmin %f step %f\n", dmin, step);
		}
		if (dmin > tol)printf("@@@@@@%f\n", dmin);

	}

	/* Coupled differential equations describing motion of photon */
	void geodesic(KerrBlackHoleWrkParam& prm, double *y, double *dydx)
	{
		double r, theta, pr, ptheta;

		r = y[0];
		theta = y[1];
		pr = y[4];
		ptheta = y[5];

		double r2 = r*r;
		double twor = 2.0*M*r;

		double sintheta, costheta;
		sincos(theta, &sintheta, &costheta);
		double cos2 = costheta*costheta;
		double sin2 = sintheta*sintheta;

		double sigma = r2 + a2*cos2;
		double delta = r2 - twor + a2; 
		double sd = sigma*delta;
		double siginv = 1.0 / sigma;
		double bot = 1.0 / sd;

		/* Prevent problems with the axis */
		if (sintheta < 1e-8)
		{
			sintheta = 1e-8;
			sin2 = 1e-16;
		}

		dydx[0] = -pr*delta*siginv;
		dydx[1] = -ptheta*siginv;
		dydx[2] = -(twor*a + (sigma - twor)*prm.L / sin2)*bot;
		dydx[3] = -(1.0 + (twor*(r2 + a2) - twor*a*prm.L)*bot);
		dydx[4] = -(((r - 1.0)*(-prm.kappa) + twor*(r2 + a2) - 2.0*a*prm.L)*bot - 2.0*pr*pr*(r - 1.0)*siginv);
		dydx[5] = -sintheta*costheta*(prm.L*prm.L / (sin2*sin2) - a2)*siginv;
	}
};


// Calculate single step of integration algorithm.
// Run a single step of integration using step h 
inline void rkstep(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, double *y, double *dydx, double h, double *yout, double *yerr)
{

	double ak[KERRBLACKHOLE_ODE_N];

	double ytemp1[KERRBLACKHOLE_ODE_N], ytemp2[KERRBLACKHOLE_ODE_N], ytemp3[KERRBLACKHOLE_ODE_N], ytemp4[KERRBLACKHOLE_ODE_N], ytemp5[KERRBLACKHOLE_ODE_N];

	for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
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

	kerrBlackHole.geodesic(prm, ytemp1, ak);

	for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp2[i] += (9.0 / 40.0) * yt;
		ytemp3[i] -= 0.9 * yt;
		ytemp4[i] += 2.5 * yt;
		ytemp5[i] += (175.0 / 512.0) * yt;
	}

	kerrBlackHole.geodesic(prm, ytemp2, ak);

	for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp3[i] += 1.2 * yt;
		ytemp4[i] -= (70.0 / 27.0) * yt;
		ytemp5[i] += (575.0 / 13824.0) * yt;
		yout[i] += (250.0 / 621.0) * yt;
		yerr[i] += ((250.0 / 621.0) - (18575.0 / 48384.0)) * yt;
	}

	kerrBlackHole.geodesic(prm, ytemp3, ak);

	for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp4[i] += (35.0 / 27.0) * yt;
		ytemp5[i] += (44275.0 / 110592.0) * yt;
		yout[i] += (125.0 / 594.0) * yt;
		yerr[i] += ((125.0 / 594.0) - (13525.0 / 55296.0)) * yt;
	}

	kerrBlackHole.geodesic(prm, ytemp4, ak);

	for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		ytemp5[i] += (253.0 / 4096.0) * yt;
		yerr[i] -= (277.0 / 14336.0) * yt;
	}

	kerrBlackHole.geodesic(prm, ytemp5, ak);

	for (int i = 0; i < KERRBLACKHOLE_ODE_N; i++)
	{
		double yt = h * ak[i];
		yout[i] += (512.0 / 1771.0) * yt;
		yerr[i] += ((512.0 / 1771.0) - 0.25) * yt;
	}
}

//適応ルンゲ・クッタ積分アルゴリズム。
//実際に実装された方法のバリエーションは、Cash-Karp.として知られています。
//ステップhtry
inline double rkqs(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, double *y, double *dydx, double htry, double escal, double *yscal, double *hdid)
{
	int i;

	double hnext;

	double errmax, h = htry, htemp;
	double yerr[KERRBLACKHOLE_ODE_N], ytemp[KERRBLACKHOLE_ODE_N];

	int loopCnt = 0;
	while (loopCnt < 2000)
	{
		// Run a single step of integration using step h 
		rkstep(prm, kerrBlackHole, y, dydx, h, ytemp, yerr);

		//誤差評価
		errmax = 0.0;
		for (i = 0; i < KERRBLACKHOLE_ODE_N; i++)
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

	memcpy(y, ytemp, KERRBLACKHOLE_ODE_N * sizeof(double));

	return hnext;
}

//光線が降着円盤か事象の地平面に当たったときの停止する必要があります。
inline void binarysearch(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, double *y, double *dydx, double hbig)
{
	double hsmall = 0.0;

	int side;
	if (y[1] > M_PI / 2.0)
	{
		side = 1;
	}
	else if (y[1] < M_PI / 2.0)
	{
		side = -1;
	}
	else
	{
		/* Already at the equator */
		return;
	}

	kerrBlackHole.geodesic(prm, y, dydx);

	//事象の地平面の外で背景に到達していない
	while ((y[0] > kerrBlackHole.Rhor) && (y[0] < kerrBlackHole.r0) && (side != 0))
	{
		double yout[KERRBLACKHOLE_ODE_N], yerr[KERRBLACKHOLE_ODE_N];

		double hdiff = hbig - hsmall;

		if (hdiff < 1e-7)
		{
			rkstep(prm, kerrBlackHole, y, dydx, hbig, yout, yerr);

			memcpy(y, yout, KERRBLACKHOLE_ODE_N * sizeof(double));

			return;
		}

		double hmid = (hbig + hsmall) / 2;

		rkstep(prm, kerrBlackHole, y, dydx, hmid, yout, yerr);

		if (side * (yout[1] - M_PI / 2.0) > 0)
		{
			hsmall = hmid;
		}
		else
		{
			hbig = hmid;
		}
	}
}

#endif
