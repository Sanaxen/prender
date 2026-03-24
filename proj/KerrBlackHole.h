#ifndef _KERRBLACKHOLE_
#define _KERRBLACKHOLE_

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "constant.h"
#include "vector3d.h"

extern bool ___debug;

using namespace prender;

#define EVENT_HORIZON	-998
#define BACK_GLOUND		-999

#define KERRBLACKHOLE_ODE_N	6
//inline double Clamp(double val, double low, double high) {
//	if (val < low) return low;
//	else if (val > high) return high;
//	else return val;
//}


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

	double traceDir;
	double rDisc;

	double geodesics_max_length;

	Vector3d camera_pos;
	Vector3d camera_dir;
	Vector3d position;
	Matrix4D Mat;
	Matrix4D invMat;
	int initial_condition;

	double color_doppler_factor_effect = 0;

	std::vector<std::string> accretion_disk_texture;
	std::string background_texture;
	double background_texture_map_coef[4];
	double background_texture_coef;

	int use_accretion_disk_temperature;


	KerrBlackHole(double x_, double y_, double z_, double rDisc_, double m, double a_, Vector3d& cameraPos, Vector3d& cameraDir)
	{
		use_accretion_disk_temperature = 0;
		background_texture = "";

		camera_pos = cameraPos;
		camera_dir = normalize(cameraDir);
		traceDir = 1.0;
		rDisc = rDisc_;
		position = Vector3d(x_, y_, z_);

		//a = -0.5; // angular momentum in theta direction (hard-coded)

		a = a_;
		M = m;
		Cartesian pos(cameraPos.x, cameraPos.y, cameraPos.z);

		pos.x -= position.x;
		pos.y -= position.y;
		pos.z -= position.z;
		//Spherical poss = pos.ToSpherical();
		Spherical poss = pos.ToBoyerLindquist(a);

		r0 = poss.r;
		theta0 = poss.th;
		phi0 = poss.ph;

		a2 = a * a;

		//Rhor = 1.0 + sqrt(1.0 - a2) + 1e-5;
		Rhor = M + sqrt(M*M - a2) + 1e-10;	//event horizon
		Rdisk = rDisc*M;
		Rmstable = inner_orbit();

#if 0
		Mat.LoadIdentity();
		invMat.LoadIdentity();

		Vector3d Z(pos.x, pos.y, pos.z);
		Z = normalize(Z);
		//Vector3d X, Y, ZZ=Z;
		//OrthonormalBasis(ZZ, Z, X, Y);

		Vector3d X, Y;
		X = normalize(cross(Vector3d(0, 1, 0), Z));
		if (X.length() < 1.0e-16)
		{
			Y = normalize(cross(Z, Vector3d(1, 0, 0)));
			X = normalize(cross(Y, Z));
		}
		else
		{
			Y = normalize(cross(Z, X));
		}

		Matrix4D mat(
			X.x, X.y, X.z, 0,
			Y.x, Y.y, Y.z, 0,
			Z.x, Z.y, Z.z, 0,
			0, 0, 0, 1.0);
		Mat = mat;

		//Mat.LoadIdentity();
		invMat = Mat.InvertMatrix();
#endif
	}

	int Setup(double x0, double y0, double z0)
	{
		//a = -0.5; // angular momentum in theta direction (hard-coded)

		Cartesian pos(x0, y0, z0);

		pos.x -= position.x;
		pos.y -= position.y;
		pos.z -= position.z;
		//Spherical poss = pos.ToSpherical();
		Spherical poss = pos.ToBoyerLindquist(a);

		r0 = poss.r;
		theta0 = poss.th;
		phi0 = poss.ph;

#if 0
		Vector3d Z(pos.x, pos.y, pos.z);
		Z = normalize(Z);
		//Vector3d X, Y, ZZ=Z;
		//OrthonormalBasis(ZZ, Z, X, Y);

		Vector3d X, Y;
		X = normalize(cross(Vector3d(0, 1, 0), Z));
		if (X.length() < 1.0e-16)
		{
			Y = normalize(cross(Z, Vector3d(1, 0, 0)));
			X = normalize(cross(Y, Z));
		}
		else
		{
			Y = normalize(cross(Z, X));
		}

		Matrix4D mat(
			X.x, Y.x, Z.x, 0,
			X.y, Y.y, Z.y, 0,
			X.z, Y.z, Z.z, 0,
			0, 0, 0, 1.0);
		Mat = mat;

		//Vector3d ss = Mat*Vector3d(1,0,0);
		//fprintf(stderr, "%f %f %f %f\n", ss.x, ss.y, ss.z, ss.length());
		//mat = Mat.InvertMatrix();
		//ss = mat*ss;
		//fprintf(stderr, "%f %f %f %f\n\n", ss.x, ss.y, ss.z, ss.length());

		double det;
		//Mat.LoadIdentity();
		bool stat = InvertMatrix_(Mat, invMat, det);
		if ( !stat )
		{
			return 0;
		}
#endif
		return 1;
	}

	double cbrt(const double x) const
	{
		return pow(x, 1.0 / 3.0);
	}
	double inner_orbit(void) const
	{
		double m2 = M*M;
		double z1, z2;

		if (a / M > 1.0)
		{
			fprintf(stderr, "a/m > 1.0!!\n");
			return 0;
		}

		if (M > 1.0e-10)
		{
			z1 = 1 + cbrt(1 - a2 / m2)*(cbrt(1 + a / M) + cbrt(1 - a / M));
			z2 = sqrt(3 * a2 / m2 + z1*z1);
		}
		else
		{
			//z1 = 1;
			//z2 = sqrt(3  + z1*z1);
			z1 = 1 + cbrt(1)*(cbrt(1) + cbrt(1));
			z2 = sqrt(z1*z1);

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

		if (energy2 <= 0)
		{
			fprintf(stderr, "energy2=%.6f s1=%.6f delta=%.6f r0=%.6f\n",
				energy2, s1, delta, r0);
		}
		double energy = sqrt(fabs(energy2))+1.0e-10;

		/* Rescale */
		y0[4] = y0[4] / energy;
		y0[5] = y0[5] / energy;

		/* Angular Momentum with E = 1 */
		prm.L = ((sigma*delta*phidot0 - 2.0*a*r0*energy)*sin2 / s1) / energy;

		prm.kappa = y0[5] * y0[5] + a2*sin2 + prm.L*prm.L / sin2;

		/* Hack - make sure everything is normalized correctly by a call to geodesic */
		geodesic(prm, y0, ydot0);


		if (1) {
#if 10
			double y1[KERRBLACKHOLE_ODE_N], dydx1[KERRBLACKHOLE_ODE_N];

			memcpy(y1, y0, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx1, ydot0, KERRBLACKHOLE_ODE_N*sizeof(double));

			//接ベクトルを求める
			NextPosition(prm, y1, dydx1, 0, 0, 0.1, &tnv, 0);
#endif
		}
	}

	void initial2(KerrBlackHoleWrkParam& prm, double *y0, double *ydot0, const Vector3d& tnv)
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

		const double POLE_THRESHOLD = 1e-8;
		bool near_pole = (sin2 < POLE_THRESHOLD * POLE_THRESHOLD);
		if (near_pole) {
			// 極付近では phidot を L/sin2 の発散を避けるため再定義
			// 実際の角運動量 L → 0 として扱う（光線が極を通過する場合）
			// tnv（接線ベクトル）から直接 L と kappa を計算する

			// 局所正規直交基底での速度成分
			double sigma = r0 * r0 + a2; // theta=0 なので cos2=1, sin2=0
			double delta = r0 * r0 - 2.0 * M * r0 + a2;

			// 極では phidot の寄与は sin(theta) * (phidot) の形で現れるので
			// phidot * sin(theta) を有限量として扱う
			double vphi_sinθ = ydot0[2] * sintheta; // 有限値として保持

			double rdot0 = ydot0[0];
			double thetadot0 = ydot0[1];
			double r2 = r0 * r0;

			double energy2 = (sigma - 2.0 * r0) * (rdot0 * rdot0 / delta + thetadot0 * thetadot0)
				+ delta * vphi_sinθ * vphi_sinθ; // sin2*phidot2 → vphi_sinθ2

			double energy = sqrt(fabs(energy2)) + 1e-10;

			y0[4] = rdot0 * sigma / delta / energy;
			y0[5] = thetadot0 * sigma / energy;

			// 極では L = 0（または微小）が物理的に自然
			// L = sin2θ * (...) → 0 as sinθ → 0
			prm.L = 0.0;

			// kappa = pθ2 + a2sin2θ + L2/sin2θ
			// L = C*sin2θ のとき L2/sin2θ = C2*sin2θ → 0
			// 極では kappa = pθ2 + 0 + 0
			prm.kappa = y0[5] * y0[5] + a2 * sin2; // sin2?0 なので pθ2 のみ

			// それでも tnv から L を推定したい場合
			// L ? sigma * delta * vphi_sinθ * sintheta / ((sigma-2r0) * energy)
			// → sin2 が掛かるので自然に 0 に近づく
		}
		else {
			if (sin2 < 1e-16)
			{
				sincos(theta0, &sintheta, &costheta);
				sincos(theta0, &sintheta, &costheta);
				sincos(theta0, &sintheta, &costheta);
				fprintf(stderr, "@@@@@@@@@@@@@@@@@:%.16f (%.16f -> %.16f)\n", theta0, sintheta, sin2);
				sintheta = 1e-16;
				sin2 = 1e-32;
			}

			double rdot0 = ydot0[0];// cos(y)*cos(x);
			double thetadot0 = ydot0[1];// sin(y) / r0;

			double r2 = r0 * r0;
			double sigma = r2 + a2 * cos2;			//Σ = r^2 + a^2 cos^2(θ)
			double delta = r2 - 2.0 * M * r0 + a2;	//Δ = r^2 - 2Mr + a^2
			double s1 = sigma - 2.0 * r0;

			//if ( fabs(sigma) < 1.0e-8 ) fprintf(stderr, "sigma %.16f\n", sigma);
			//if ( fabs(delta) < 1.0e-8 ) fprintf(stderr, "delta %.16f\n", delta);
			//if ( fabs(s1) < 1.0e-8 ) fprintf(stderr, "s1 %.16f\n", s1);

			y0[4] = rdot0 * sigma / delta;
			y0[5] = thetadot0 * sigma;


			double phidot0 = ydot0[2];
			double energy2 = s1 * (rdot0 * rdot0 / delta + thetadot0 * thetadot0)
				+ delta * sin2 * phidot0 * phidot0;

			if (energy2 <= 0)
			{
				fprintf(stderr, "energy2=%.6f s1=%.6f delta=%.6f r0=%.6f\n",
					energy2, s1, delta, r0);
			}
			double energy = sqrt(fabs(energy2)) + 1.0e-10;

			if (fabs(energy) < 1.0e-8) fprintf(stderr, "energy %.16f\n", energy);

			/* Rescale */
			y0[4] = y0[4] / energy;
			y0[5] = y0[5] / energy;

			/* Angular Momentum with E = 1 */
			prm.L = ((sigma * delta * phidot0 - 2.0 * a * r0 * energy) * sin2 / s1) / energy;

			prm.kappa = y0[5] * y0[5] + a2 * sin2 + prm.L * prm.L / sin2;

			if (___debug)
			{
				printf("frame: L=%.6f kappa=%.6f energy=%.6f s1=%.6f delta=%.6f\n",
					prm.L, prm.kappa, energy, s1, delta);
				printf("rdot0=%.6f thetadot0=%.6f phidot0=%.6f\n",
					rdot0, thetadot0, phidot0);
				printf("dir=%.6f,%.6f,%.6f\n", tnv.x, tnv.y, tnv.z);
			}
		}
		//if ( fabs(prm.L) < 1.0e-8 ) fprintf(stderr, "prm.L %.16f\n", prm.L);
		//if ( fabs(prm.kappa) < 1.0e-8 ) fprintf(stderr, "prm.kappa %.16f\n", prm.kappa);
		/* Hack - make sure everything is normalized correctly by a call to geodesic */
		geodesic(prm, y0, ydot0);


		if (1) {
#if 0
			Verctor3d tnv1;
			double y1[KERRBLACKHOLE_ODE_N], dydx1[KERRBLACKHOLE_ODE_N];

			memcpy(y1, y0, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx1, ydot0, KERRBLACKHOLE_ODE_N*sizeof(double));

			//接ベクトルを求める
			NextPosition(prm, y1, dydx1, 0, 0, 0.1, &tnv1, 0);
#endif
		}
	}

	inline Vector3d WorldAxisSystem(double xx, double yy, double zz) const
	{
		return Vector3d(xx, yy, zz) + position;
	}

	inline Vector3d WorldAxisSystem(Spherical& pos) const
	{
		Cartesian q = pos.ToCartesian();
		return WorldAxisSystem(q.x, q.y, q.z);
	}
	inline Vector3d WorldAxisSystem(Cartesian& pos) const
	{
		return WorldAxisSystem(pos.x, pos.y, pos.z);
	}

	inline Vector3d Coordinate_transformation(double xx, double yy, double zz) const
	{
		return Coordinate_transformation(Vector3d(xx, yy, zz));
	}
	inline Vector3d Coordinate_transformation(Vector3d& p) const
	{
		return p - position;
	}
	inline Vector3d Coordinate_transformation(const Vector3d& p) const
	{
		return p - position;
	}

	//微分方程式の次の位置と接ベクトルを計算する
	double NextPosition(KerrBlackHoleWrkParam& prm, double* y, double* dydx, double* y_next, double* dydx_next, double htry, Vector3d* tnv, double *hdid)
	{
		double escal = 1e11, hnext = 0.0;

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
			//Cartesian c_org = org.ToCartesian();
			Cartesian c_org = org.ToBoyerLindquist(a);

			Spherical p(y1[0], y1[1], y1[2]);

			//Cartesian q = p.ToCartesian();
			Cartesian q = p.ToBoyerLindquist(a);

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
			memcpy(y_next, y1, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx_next, dydx1, KERRBLACKHOLE_ODE_N*sizeof(double));
		}

		return hnext;
	}

	inline int RayToOrdinaryDifferentialEquationInitialFnc(KerrBlackHoleWrkParam& prm, Vector3d& dir, double* y, double* dydx, double x_, double y_, double z_, double f[3])
	{
		//const double cost = a*cos(theta0 + x_);
		//const double r = y[0];
		//const double sigma = r*r + a*a*cost*cost;
		//const double R = sqrt(a2 + r*r);
		//const double v = -sin(z_)*cos(y_);
		//const double zdot = 1.0;

		//const double sint = sin(theta0 + x_);
		//const double rdot0 = zdot*(r*R*v*sint - R*R*cos(z_)*cost) / sigma;
		//const double thetadot0 = zdot*(r*sint*cos(z_) + R*v*cost) / sigma;
		//const double phidot0 = zdot*sin(z_)*sin(y_) / (R*sint);

		//dydx[0] = rdot0;
		//dydx[1] = thetadot0;
		//dydx[2] = phidot0;


		const double r = y[0];
		const double theta = y[1];
		const double phi = y[2];

#if 10
		Cartesian c(dir.x, dir.y, dir.z);
		const Spherical s = c.ToSpherical();
		//const Spherical s = c.ToBoyerLindquist(a);
		const double theta_obs = s.th;
		const double phi_obs = s.ph;
		const double Phi = phi - phi_obs;

		const double cost = cos(theta);
		const double sintobs = sin(theta_obs);

		const double rr = r*r;
		const double sigma = rr + (a*cost)*(a*cost);
		const double R = sqrt(a2 + rr);
		const double v = -sintobs*cos(Phi);
		const double zdot = 1.;

		const double inv_sigma = 1.0 / sigma;
		double sint = sin(theta);

		const double rdot0 = x_*zdot*(r*R*v*sint - R*R*cos(theta_obs)*cost) * inv_sigma;
		const double thetadot0 = y_*zdot*(r*sint*cos(theta_obs) + R*v*cost) * inv_sigma;
		
#if 0
		if (fabs(sint) < 1.0e-16)
		{
			fprintf(stderr, "sint %.16f\n", sint);
			if (sint < 0) sint = -1.0e-16;
			else  sint = 1.0e-16;
		}
		const double phidot0 = z_*zdot*sintobs*sin(Phi) / (R*sint);
#else
		double phidot0;
		const double POLE_THRESHOLD = 1e-8;
		if (fabs(sint) < POLE_THRESHOLD)
		{
			// 極では phi 方向の速度は不定だが、
			// 物理的な寄与 sin2θ * phidot2 → 0 なので 0 で近似
			phidot0 = 0.0;
		}
		else
		{
			phidot0 = z_ * zdot * sintobs * sin(Phi) / (R * sint);
		}
#endif

#else
		const double sin_th = sin(theta);
		const double cos_th = cos(theta);
		const double sin_ph = sin(phi);
		const double cos_ph = cos(phi);
		const Vector3d e_r(sin_th * cos_ph, sin_th * sin_ph, cos_th);
		const Vector3d e_th(cos_th * cos_ph, cos_th * sin_ph, -sin_th);
		const Vector3d e_ph(-sin_ph, cos_ph, 0.0);

		// dir をローカルフレームに投影
		const double dir_r = dot(dir, e_r);
		const double dir_th = dot(dir, e_th);
		const double dir_ph = dot(dir, e_ph);

		// cos/sin を直接計算（atan2の不連続を回避）
		const double R_ph = sqrt(dir_r * dir_r + dir_ph * dir_ph);
		const double cos_Phi = (R_ph > 1.0e-16) ? dir_r / R_ph : 1.0;
		const double sin_Phi = (R_ph > 1.0e-16) ? dir_ph / R_ph : 0.0;
		const double cos_theta_obs = dir_th;
		const double sin_theta_obs = R_ph;

		const double cost = cos_th;  // cos(theta) は既に計算済み

		const double rr = r * r;
		const double sigma = rr + (a * cos_th) * (a * cos_th);
		const double R = sqrt(a2 + rr);
		const double v = -sin_theta_obs * cos_Phi;
		const double zdot = 1.;

		const double inv_sigma = 1.0 / sigma;
		double sint = sin_th;  // sin(theta) は既に計算済み

		if (fabs(sint) < 1.0e-16)
		{
			fprintf(stderr, "sint %.16f\n", sint);
			if (sint < 0) sint = -1.0e-16;
			else  sint = 1.0e-16;
		}
		const double rdot0 = x_ * zdot * (r * R * v * sint - R * R * cos_theta_obs * cost) * inv_sigma;
		const double thetadot0 = y_ * zdot * (r * sint * cos_theta_obs + R * v * cost) * inv_sigma;
		const double phidot0 = z_ * zdot * sin_theta_obs * sin_Phi / (R * sint);
#endif

		dydx[0] = rdot0;
		dydx[1] = thetadot0;
		dydx[2] = phidot0;

		initial2(prm, y, dydx, traceDir*dir);

		double hdid = 0;
		Vector3d tnv;

		NextPosition(prm, y, dydx, 0, 0, 0.001, &tnv, &hdid);

		f[0] = (dir.x - tnv.x);
		f[1] = (dir.y - tnv.y);
		f[2] = (dir.z - tnv.z);

		f[0] *= f[0];
		f[1] *= f[1];
		f[2] *= f[2];
		return 0;
	}

#define MAX_THREADS_LC 1

	//Rayベクトルと一致する微分方程式初期条件を求める
	int RayToOrdinaryDifferentialEquationInitial(int depth, KerrBlackHoleWrkParam& prm, Vector3d& dir, double* y, double* dydx)
	{
		////Vector3d nrm = -1.0*normalize(Spherical(r0, theta0, phi0).ToVector());

		//if (dot(dir, nrm) < 0)
		//{
		//	return -1;
		//}

		//y[0] = r0;
		//y[1] = theta0;
		//y[2] = phi0;
		y[3] = 0;
		y[4] = 0;
		y[5] = 0;

#if 0
		Vector3d X, Y;
		X = normalize(cross(Vector3d(0, 1, 0), nrm));
		if (X.length() < 1.0e-8)
		{
			Y = normalize(cross(nrm, Vector3d(1, 0, 0)));
			X = normalize(cross(Y, nrm));
		}
		else
		{
			Y = normalize(cross(nrm, X));
		}

		double x00 = dot(X, dir)*-1.0;
		double y00 = dot(Y, dir)*-1.0;

		double rdot0 = cos(y00)*cos(x00);
		double thetadot0 = sin(y00) / r0;

		dydx[0] = rdot0;
		dydx[1] = thetadot0;
		dydx[2] = cos(y00)*sin(x00) / (r0*sin(theta0));

		initial2(prm, y, dydx, traceDir*dir);

		double y1[KERRBLACKHOLE_ODE_N];
		double dydx1[KERRBLACKHOLE_ODE_N];
		memcpy(y1, y, KERRBLACKHOLE_ODE_N*sizeof(double));
		memcpy(dydx1, dydx, KERRBLACKHOLE_ODE_N*sizeof(double));
		double hdid = 0;
		Vector3d tnv;

		//次の位置を求める
		//geodesic(prm, y1, dydx1);
		NextPosition(prm, y1, dydx1, y1, dydx1, 0.001, &tnv, &hdid);

		//Spherical ss(r0, theta0, phi0);
		//Spherical sss(y1[0], y1[1], y1[2]);
		//double inner = (sss.ToVector() - (ss.ToVector()+0.001*dir)).length();
		//double inner = (tnv - dir).length();
		//fprintf(stderr, "%f\n", inner);

		return 0;
#else
		if (initial_condition == 0)
		{
			double* y0 = y;
			double* ydot0 = dydx;
			double A = a;


			double r = y0[0];
			double theta = y0[1];
			double phi = y0[2];


			//Spherical ss(r, theta, phi);
			//double z = ss.ToCartesian().z;
			//double u = r*r -a2;
			//r =  sqrt((u+sqrt(u*u+(2.*A*z)*(2.*A*z)))/2.);   
			//theta = acos(Clamp(z / r, -1.0, 1.0));
			//phi = atan2(ss.ToCartesian().y / r, ss.ToCartesian().x / r);

#if 10
			Cartesian c(dir.x, dir.y, dir.z);
			//Spherical s = c.ToBoyerLindquist(a);
			Spherical s = c.ToSpherical();
			double theta_obs = s.th;
			double phi_obs = s.ph;
			double Phi = phi - phi_obs;

			//double d = (sqrt(r*r + a2) - ss.ToCartesian().z)*sin(theta_obs) - ss.ToCartesian().y*cos(theta_obs);
			//double xx = d*cos(theta_obs) - ss.ToCartesian().x*sin(theta_obs);
			//double yy = d*sin(theta_obs) + ss.ToCartesian().x*cos(theta_obs);
			//double zz = (r - ss.ToCartesian().z)*cos(theta_obs) + ss.ToCartesian().y*sin(theta_obs);
			//
			//y0[0] = sqrt((u + sqrt(u*u + (2.*A*z)*(2.*A*z))) / 2.);
			//y0[1] = Clamp(acos(z / r), -1.0, 1.0);
			//y0[2] = atan2(ss.ToCartesian().y, ss.ToCartesian().x);

			const double cost = cos(theta);
			const double sintobs = sin(theta_obs);

			const double rr = r*r;
			const double sigma = rr + (a*cost)*(a*cost);
			const double R = sqrt(a2 + rr);
			const double v = -sintobs*cos(Phi);
			const double zdot = 1.;

			const double inv_sigma = 1.0 / sigma;
			double sint = sin(theta);

			const double rdot0 = zdot*(r*R*v*sint - R*R*cos(theta_obs)*cost) * inv_sigma;
			const double thetadot0 = zdot*(r*sint*cos(theta_obs) + R*v*cost) * inv_sigma;
			
#if 0
			if (fabs(sint) < 1.0e-16)
			{
				fprintf(stderr, "sint %.16f\n", sint);
				if (sint < 0) sint = -1.0e-16;
				else  sint = 1.0e-16;
			}
			const double phidot0 = zdot*sintobs*sin(Phi) / (R*sint);
#else
			double phidot0;
			const double POLE_THRESHOLD = 1e-8;

			if (fabs(sint) < POLE_THRESHOLD)
			{
				// 極では phi 方向の速度は不定だが、
				// 物理的な寄与 sin2θ * phidot2 → 0 なので 0 で近似
				phidot0 = 0.0;
			}
			else
			{
				phidot0 = zdot * sintobs * sin(Phi) / (R * sint);
			}
#endif
#else
			const double sin_th = sin(theta);
			const double cos_th = cos(theta);
			const double sin_ph = sin(phi);
			const double cos_ph = cos(phi);
			const Vector3d e_r(sin_th * cos_ph, sin_th * sin_ph, cos_th);
			const Vector3d e_th(cos_th * cos_ph, cos_th * sin_ph, -sin_th);
			const Vector3d e_ph(-sin_ph, cos_ph, 0.0);

			// dir をローカルフレームに投影
			const double dir_r = dot(dir, e_r);
			const double dir_th = dot(dir, e_th);
			const double dir_ph = dot(dir, e_ph);

			// cos/sin を直接計算（atan2の不連続を回避）
			const double R_ph = sqrt(dir_r * dir_r + dir_ph * dir_ph);
			const double cos_Phi = (R_ph > 1.0e-16) ? dir_r / R_ph : 1.0;
			const double sin_Phi = (R_ph > 1.0e-16) ? dir_ph / R_ph : 0.0;
			const double cos_theta_obs = dir_th;
			const double sin_theta_obs = R_ph;

			const double cost = cos_th;

			const double rr = r * r;
			const double sigma = rr + (a * cos_th) * (a * cos_th);
			const double R = sqrt(a2 + rr);
			const double v = -sin_theta_obs * cos_Phi;
			const double zdot = 1.;

			const double inv_sigma = 1.0 / sigma;
			double sint = sin_th;
			if (fabs(sint) < 1.0e-16)
			{
				fprintf(stderr, "sint %.16f\n", sint);
				if (sint < 0) sint = -1.0e-16;
				else  sint = 1.0e-16;
			}

			const double rdot0 = zdot * (r * R * v * sint - R * R * cos_theta_obs * cost) * inv_sigma;
			const double thetadot0 = zdot * (r * sint * cos_theta_obs + R * v * cost) * inv_sigma;
			const double phidot0 = zdot * sin_theta_obs * sin_Phi / (R * sint); 
#endif

			dydx[0] = rdot0;
			dydx[1] = thetadot0;
			dydx[2] = phidot0;

			initial2(prm, y, dydx, traceDir*dir);

			double check_x = dydx[0] * sint * cost - dydx[1] * cost * cos(Phi) - dydx[2] * sin(Phi);
			double check_y = dydx[0] * sint * sin(Phi) - dydx[1] * cost * sin(Phi) + dydx[2] * cos(Phi);
			double check_z = dydx[0] * cost + dydx[1] * sint;
			Vector3d dydx_dir(check_x, check_y, check_z);
			double alignment = dot(normalize(dydx_dir), dir);
			//printf("alignment=%.6f (1=正常, -1=反転)\n", alignment);

			return 0;

		}

		if (initial_condition == 1)
		{
			//int threadNum = omp_get_max_threads();

			//omp_set_num_threads(MAX_THREADS_LC);

			double y1[MAX_THREADS_LC][KERRBLACKHOLE_ODE_N];
			double dydx1[MAX_THREADS_LC][KERRBLACKHOLE_ODE_N];

			int stat[MAX_THREADS_LC];
			int count[MAX_THREADS_LC];
			float err[MAX_THREADS_LC];
			double x__[MAX_THREADS_LC];
			double y__[MAX_THREADS_LC];
			double z__[MAX_THREADS_LC];

			MTRandom rnd[MAX_THREADS_LC];
			for (int i = 0; i < MAX_THREADS_LC; i++)
			{
				count[i] = 0;
				stat[i] = -1;
				err[i] = 99999999.0f;
				rnd[i].seed(i + 1);
			}

			bool errJ = false;

			const double ew = 0.01;	//0.01
			//#pragma omp parallel for
			for (int kk = 0; kk < 500000; kk++)
			{
				errJ = false;
				int id = 0;
				//int id = omp_get_thread_num();

				if (stat[id] == 0) continue;

				KerrBlackHoleWrkParam prm0;
				double x_ = 1.0 + (-1.0 + 2.0*rnd[id].next01()) * ew;
				double y_ = 1.0 + (-1.0 + 2.0*rnd[id].next01()) * ew;
				double z_ = 1.0 + (-1.0 + 2.0*rnd[id].next01()) * ew;
				double fnc0[3];
				double fnc1[3];
				double fnc2[3];
				double dh = 0.0001;
				double J[3][3];

				const double inv_dh = 1.0 / dh;
				float delta[3] = { 0.f, 0.f, 0.f };

				for (int ii = 0; ii < 50; ii++)	// 500
				{

					//bool flag = false;
					//for (int i = 0; i < MAX_THREADS_LC; i++)
					//{
					//	if (stat[i] == 0)
					//	{
					//		flag = true;
					//		break;
					//	}
					//}
					//if (flag) break;

					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_, y_, z_, fnc0);

					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_ + dh, y_, z_, fnc1);

					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_ - dh, y_, z_, fnc2);

					//∂f[0]/∂x ∂f[1]/∂x ∂f[2]/∂x
					J[0][0] = 0.5*(fnc1[0] - fnc2[0]) * inv_dh;
					J[1][0] = 0.5*(fnc1[1] - fnc2[1]) * inv_dh;
					J[2][0] = 0.5*(fnc1[2] - fnc2[2]) * inv_dh;
					//fprintf(stderr, "%f %f %f\n", J[0][0], J[1][0], J[2][0]);

					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_, y_ + dh, z_, fnc1);

					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_, y_ - dh, z_, fnc2);

					//∂f[0]/∂y ∂f[1]/∂y ∂f[2]/∂y
					J[0][1] = 0.5*(fnc1[0] - fnc2[0]) * inv_dh;
					J[1][1] = 0.5*(fnc1[1] - fnc2[1]) * inv_dh;
					J[2][1] = 0.5*(fnc1[2] - fnc2[2]) * inv_dh;
					//fprintf(stderr, "%f %f %f\n", J[0][1], J[1][1], J[2][1]);


					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_, y_, z_ + dh, fnc1);

					y1[id][0] = y[0]; y1[id][1] = y[1];	y1[id][2] = y[2];
					RayToOrdinaryDifferentialEquationInitialFnc(prm0, dir, y1[id], dydx1[id], x_, y_, z_ - dh, fnc2);

					//∂f[0]/∂z ∂f[1]/∂z ∂f[2]/∂z
					J[0][2] = 0.5*(fnc1[0] - fnc2[0]) * inv_dh;
					J[1][2] = 0.5*(fnc1[1] - fnc2[1]) * inv_dh;
					J[2][2] = 0.5*(fnc1[2] - fnc2[2]) * inv_dh;
					//fprintf(stderr, "%f %f %f\n", J[0][2], J[1][2], J[2][2]);

					/*
							  ∂f[0]/∂x  ∂f[0]/∂y ∂f[0]/∂z
							  J =[ ∂f[1]/∂x  ∂f[1]/∂y ∂f[1]/∂z ]
							  ∂f[2]/∂x  ∂f[2]/∂y ∂f[2]/∂z
							  */

					//double detJ;
					//detJ = J[0][0] * J[1][1] * J[2][2];
					//detJ += J[1][0] * J[2][1] * J[0][2];
					//detJ += J[2][0] * J[0][1] * J[1][2];
					//detJ -= J[2][0] * J[1][1] * J[0][2];
					//detJ -= J[1][0] * J[0][1] * J[2][2];
					//detJ -= J[0][0] * J[2][1] * J[1][2];

					//if (fabs(detJ) < 1.0e-16)
					//{
					//	fprintf(stderr, "detJ %.16f\n", detJ);
					//	break;
					//}
					double inv_J[3][3]; //ここに逆行列が入る

					{
						const int n = 3;  //配列の次数

						//単位行列を作る
						if (n == 3)
						{
							inv_J[0][0] = 1.0;
							inv_J[0][1] = 0.0;
							inv_J[0][2] = 0.0;
							inv_J[1][0] = 0.0;
							inv_J[1][1] = 1.0;
							inv_J[1][2] = 0.0;
							inv_J[2][0] = 0.0;
							inv_J[2][1] = 0.0;
							inv_J[2][2] = 1.0;
						}
						else
						{
							for (int i = 0; i < n; i++)
							{
								for (int j = 0; j < n; j++)
								{
									inv_J[i][j] = (i == j) ? 1.0 : 0.0;
								}
							}
						}

						//掃き出し法
						for (int i = 0; i < n; i++)
						{
							if (fabs(J[i][i]) < 1.0e-16)
							{
								errJ = true;
								break;
							}
							const double buf = 1.0 / J[i][i];

							if (n == 3)
							{
								J[i][0] *= buf;
								J[i][1] *= buf;
								J[i][2] *= buf;
								inv_J[i][0] *= buf;
								inv_J[i][1] *= buf;
								inv_J[i][2] *= buf;

								for (int j = 0; j < n; j++)
								{
									if (i != j)
									{
										const double buf = J[j][i];
										J[j][0] -= J[i][0] * buf;
										J[j][1] -= J[i][1] * buf;
										J[j][2] -= J[i][2] * buf;
										inv_J[j][0] -= inv_J[i][0] * buf;
										inv_J[j][1] -= inv_J[i][1] * buf;
										inv_J[j][2] -= inv_J[i][2] * buf;
									}
								}
							}
							else
							{
								for (int j = 0; j < n; j++)
								{
									J[i][j] *= buf;
									inv_J[i][j] *= buf;
								}
								for (int j = 0; j < n; j++)
								{
									if (i != j)
									{
										const double buf = J[j][i];
										for (int k = 0; k < n; k++)
										{
											J[j][k] -= J[i][k] * buf;
											inv_J[j][k] -= inv_J[i][k] * buf;
										}
									}
								}
							}
						}
					}
					if (errJ)
					{
						break;
					}

					delta[0] = inv_J[0][0] * (-fnc0[0]) + inv_J[0][1] * (-fnc0[1]) + inv_J[0][2] * (-fnc0[2]);
					delta[1] = inv_J[1][0] * (-fnc0[0]) + inv_J[1][1] * (-fnc0[1]) + inv_J[1][2] * (-fnc0[2]);
					delta[2] = inv_J[2][0] * (-fnc0[0]) + inv_J[2][1] * (-fnc0[1]) + inv_J[2][2] * (-fnc0[2]);

					x_ += delta[0] * 1;
					y_ += delta[1] * 1;
					z_ += delta[2] * 1;

					if (fabs(x_) > 1.0e2 || fabs(y_) > 1.0e2 || fabs(z_) > 1.0e2)
					{
						break;
					}

					//err[id] = sqrtf(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
					err[id] = sqrtf(fnc0[0] * fnc0[0] + fnc0[1] * fnc0[1] + fnc0[2] * fnc0[2]);
					if (err[id] < 1.0e-10)
					{
						if ( ii > 40 || kk > 200 )fprintf(stderr, "depth %d 収束[%d](%d %d)!!\n", depth, id, kk, ii);
						x__[id] = x_;
						y__[id] = y_;
						z__[id] = z_;
						stat[id] = 0;
						count[id] = ii;
						break;
					}
				}
			}

			for (int i = 0; i < MAX_THREADS_LC; i++)
			{
				if (stat[i] == 0)
				{
					double fnc0[3];
					RayToOrdinaryDifferentialEquationInitialFnc(prm, dir, y, dydx, x__[i], y__[i], z__[i], fnc0);
					//omp_set_num_threads(threadNum);
					//fprintf(stderr, "deptah %d 収束[%d]!!\n", depth, i);
					return 0;
				}
			}
			//omp_set_num_threads(threadNum);
			fprintf(stderr, "depth %d 収束しなかった\n", depth);
			return -1;
		}

		Spherical s(r0, theta0, phi0);
		int threadNum = omp_get_max_threads();

		omp_set_num_threads(MAX_THREADS_LC);

		double y0[MAX_THREADS_LC + 1][KERRBLACKHOLE_ODE_N];
		double dydx0[MAX_THREADS_LC + 1][KERRBLACKHOLE_ODE_N];

		//fprintf(stderr, "=>\n");
		double minval[MAX_THREADS_LC + 1];
		int n = 3;

		MTRandom rnd[MAX_THREADS_LC];
		for (int i = 0; i < MAX_THREADS_LC; i++)
		{
			rnd[i].seed(i + 1);
		}

		int cnt = 0;
		for (int i = 0; i < MAX_THREADS_LC + 1; i++) minval[i] = 999999999.0;

	A:;
		const double costheta0 = cos(theta0);
		const double sintheta0 = sin(theta0);
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			int id = omp_get_thread_num();

			double y1[KERRBLACKHOLE_ODE_N];
			double dydx1[KERRBLACKHOLE_ODE_N];

			double y2[KERRBLACKHOLE_ODE_N];
			double dydx2[KERRBLACKHOLE_ODE_N];

			KerrBlackHoleWrkParam prm;

			memcpy(y1, y, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx1, dydx, KERRBLACKHOLE_ODE_N*sizeof(double));

			//double x00 = (-1.0 + 2.0*rnd[id].next01()) * 3;
			double y00 = (-1.0 + 2.0*rnd[id].next01()) * 3;
			double z00 = (-1.0 + 2.0*rnd[id].next01()) * 3;

			double cosz00 = cos(z00);
			double sinz00 = sin(z00);
			double cosy00 = cos(y00);
			double siny00 = sin(y00);

			double r = y1[0];
			double sigma = r*r + (a* costheta0)*(a* costheta0);
			double R = sqrt(a2 + r*r);
			double v = -sinz00*cosy00;
			double zdot = 1.;

			double rdot0 = zdot*(r*R*v*sintheta0 - R*R*cosz00* costheta0) / sigma;
			double thetadot0 = zdot*(r* sintheta0 *cosz00 + R*v* costheta0) / sigma;
			double phidot0 = zdot*sinz00*siny00 / (R* sintheta0);

			dydx1[0] = rdot0;
			dydx1[1] = thetadot0;
			dydx1[2] = phidot0;


			memcpy(y2, y1, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx2, dydx1, KERRBLACKHOLE_ODE_N*sizeof(double));

			initial2(prm, y2, dydx2, traceDir*dir);

			double hdid = 0;
			Vector3d tnv;

			double y3[KERRBLACKHOLE_ODE_N];
			double dydx3[KERRBLACKHOLE_ODE_N];
			//次の位置を求める
			NextPosition(prm, y2, dydx2, y3, dydx3, 0.01, &tnv, &hdid);

			Spherical ss(y3[0], y3[1], y3[2]);
			//double inner = ( ss.ToVector() -(s.ToVector() + 0.5*dir)).length();
			double inner = fabs(1.0 - dot(tnv, dir));
			if (inner < minval[id])
			{
				minval[id] = inner;
				memcpy(y0[id], y1, KERRBLACKHOLE_ODE_N*sizeof(double));
				memcpy(dydx0[id], dydx1, KERRBLACKHOLE_ODE_N*sizeof(double));
			}
		}


		for (int i = 0; i < MAX_THREADS_LC; i++)
		{
			if (minval[i] < minval[MAX_THREADS_LC])
			{
				minval[MAX_THREADS_LC] = minval[i];
				memcpy(y0[MAX_THREADS_LC], y0[i], KERRBLACKHOLE_ODE_N*sizeof(double));
				memcpy(dydx0[MAX_THREADS_LC], dydx0[i], KERRBLACKHOLE_ODE_N*sizeof(double));
			}
		}

		cnt++;
		if (minval[MAX_THREADS_LC] > 0.0001 && cnt < 4000)
		{
			goto A;
		}

		omp_set_num_threads(threadNum);

		//if (minval[MAX_THREADS_LC] > 0.1)
		//{
		//	fprintf(stderr, "%d NG %f\n\n", cnt, minval[MAX_THREADS_LC]);
		//	return -1;
		//}

		//if ( cnt >= 400 )
		//{
		//	fprintf(stderr, "%d NG %f\n\n", cnt, minval[MAX_THREADS_LC]);
		//	return -1;
		//}
		//if ( cnt > 3000 ) fprintf(stderr, "%d OK %f\n\n", cnt, minval[MAX_THREADS_LC]);
		//else fprintf(stderr, "%d NG %f\n\n", cnt, minval[MAX_THREADS_LC]);
		memcpy(y, y0[MAX_THREADS_LC], KERRBLACKHOLE_ODE_N*sizeof(double));
		memcpy(dydx, dydx0[MAX_THREADS_LC], KERRBLACKHOLE_ODE_N*sizeof(double));
		initial2(prm, y, dydx, traceDir*dir);

		return 0;
#endif
	}

	/* Coupled differential equations describing motion of photon */
	inline void geodesic(KerrBlackHoleWrkParam& prm, double *y, double *dydx)
	{
		const double r = y[0];
		const double theta = y[1];
		const double pr = y[4];
		const double ptheta = y[5];

		const double r2 = r*r;
		const double twor = 2.0*M*r;

		double sintheta = sin(theta);
		double sin2 = sintheta*sintheta;

		const double costheta = cos(theta);
		const double cos2 = costheta*costheta;

		const double sigma = r2 + a2*cos2;
		const double delta = r2 - twor + a2;
		const double sd = sigma*delta;
		const double siginv = 1.0 / sigma;
		const double bot = 1.0 / sd;

		/* Prevent problems with the axis */
		//if ( fabs(sintheta) < 1e-8)
		//{
		//	fprintf(stderr, "sintheta %.16f\n", sintheta);
		//	sintheta = 1e-8;
		//	sin2 = 1e-16;
		//}

		if (fabs(sin2) < 1.0e-16)
		{
			if (sin2 < 0.0) sin2 = -1.0e-16;
			else sin2 = 1.0e-16;
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

	h *= kerrBlackHole.traceDir;

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


	const int loopMax = 40;

	int loopCnt = 0;
	while (loopCnt < loopMax)
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
	if (loopCnt == loopMax)
	{
		fprintf(stderr, "integration step h[%f] error[errmax %f]\n", h, errmax);
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
inline void binarysearch(KerrBlackHoleWrkParam& prm, KerrBlackHole& kerrBlackHole, double *y, double *dydx, double hbig, Vector3d* tnv)
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

			Spherical org(y[0], y[1], y[2]);
			Cartesian c_org = org.ToCartesian();

			Spherical p(yout[0], yout[1], yout[2]);
			Cartesian q = p.ToCartesian();

			tnv->x = q.x - c_org.x;
			tnv->y = q.y - c_org.y;
			tnv->z = q.z - c_org.z;

			*tnv = normalize(*tnv);

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

// ケプラー回転速度を計算（GM=1の単位系）
// v = sqrt(GM/r) = sqrt(1/r)
// 光速 c=1 の単位系なので β = v/c = v
inline double getKeplerianVelocity(double radius) {
	// 簡略化: v/c = sqrt(1/r)
	// 実際にはブラックホール近傍では一般相対論的補正が必要
	return std::sqrt(1.0 / radius);
}
// 降着円盤上の点の速度ベクトルを計算
inline Vector3d getDiskVelocity(const Vector3d& position, const Vector3d& disk_center, const Vector3d& rotation_axis, const double rotation_sign=1.0)
{
	// ブラックホール中心からの相対位置
	Vector3d r_vec = position - disk_center;

	// 円盤平面上での距離
	double r = r_vec.length();

	if (r < 1e-6) return Vector3d(0, 0, 0);

	// ケプラー速度の大きさ
	double v_magnitude = getKeplerianVelocity(r);

	// 速度方向 = rotation_axis × r_vec (外積)
	Vector3d velocity_direction(
		rotation_axis.y * r_vec.z - rotation_axis.z * r_vec.y,
		rotation_axis.z * r_vec.x - rotation_axis.x * r_vec.z,
		rotation_axis.x * r_vec.y - rotation_axis.y * r_vec.x
	);

	velocity_direction = normalize(velocity_direction);

	// 回転方向を適用
	return velocity_direction * (v_magnitude * rotation_sign);
}

// 観測者から見たドップラー因子を計算
inline double calculateDopplerFactorFromObserver(const Vector3d& disk_point,
	const Vector3d& disk_velocity,
	const Vector3d& observer_position) {
	// 観測方向（観測者から円盤の点へ）
	Vector3d to_disk = disk_point - observer_position;
	//Vector3d to_disk = observer_position - disk_point;
	Vector3d viewing_direction = normalize(to_disk);

	// 速度の大きさ（光速の割合）
	double beta = disk_velocity.length();

	if (beta < 1e-6) return 1.0;

	// 速度方向
	Vector3d velocity_direction = normalize(disk_velocity);

	// 視線方向と速度方向の内積
	// 注意: ここでは視線は「観測者→円盤」なので、
	// 円盤が観測者に近づく場合は cos_theta > 0
	double cos_theta = dot(viewing_direction,velocity_direction);

	// 相対論的ドップラー因子
	// 近づく場合: doppler > 1 (青方偏移)
	// 遠ざかる場合: doppler < 1 (赤方偏移)
	double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	double doppler = 1.0 / (gamma * (1.0 - beta * cos_theta));

	return doppler;
}

// 標準降着円盤の温度分布（幾何学単位系 G = c = 1）
inline double accretionDiskTemperatureGeom(
	double r,           // ブラックホールからの距離 [M単位]
	double M,           // ブラックホール質量（通常は1に規格化）
	double r_in         // 内縁半径（通常はISCO） [M単位]
) {

	// 降着率（エディントン降着率の10%）
	double Mdot_Edd = (4.0 * M_PI * 1.67262192369e-27) /(0.1 *  6.6524587158e-29);
	double Mdot = 0.1 * Mdot_Edd;

	// 幾何学単位系では: G = c = 1
	// 温度を実際の物理単位に変換する際、M_SI が必要

	if (r <= r_in) return 0.0;  // 内縁より内側

	// 無次元温度パラメータ（実際の温度を得るには M_SI が必要）
	// T^4 ∝ (M_dot / r^3) * (1 - sqrt(r_in/r))

	double factor = Mdot / (r * r * r);
	double innerBoundary = 1.0 - std::sqrt(r_in / r);

	if (innerBoundary < 0) innerBoundary = 0;

	// 無次元温度の4乗
	double T4_dimensionless = factor * innerBoundary;

	return std::pow(T4_dimensionless, 0.25);
}

// 黒体放射の近似的なRGB変換（温度[K] -> RGB）
inline Color temperatureToRGB(double temp) {
	// 簡略化した黒体放射の色近似
	// 実際のプランク関数を使う場合はより正確
	temp = Clamp(temp, 1000.0, 40000.0);
	temp /= 100.0;

	double r, g, b;

	// 赤チャンネル
	if (temp <= 66) {
		r = 255;
	}
	else {
		r = temp - 60;
		r = 329.698727446 * std::pow(r, -0.1332047592);
		r = Clamp(r, 0.0, 255.0);
	}

	// 緑チャンネル
	if (temp <= 66) {
		g = temp;
		g = 99.4708025861 * std::log(g) - 161.1195681661;
		g = Clamp(g, 0.0, 255.0);
	}
	else {
		g = temp - 60;
		g = 288.1221695283 * std::pow(g, -0.0755148492);
		g = Clamp(g, 0.0, 255.0);
	}

	// 青チャンネル
	if (temp >= 66) {
		b = 255;
	}
	else if (temp <= 19) {
		b = 0;
	}
	else {
		b = temp - 10;
		b = 138.5177312231 * std::log(b) - 305.0447927307;
		b = Clamp(b, 0.0, 255.0);
	}

	return Color(r / 255.0, g / 255.0, b / 255.0);
}
#endif
