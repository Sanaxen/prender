#ifndef _BRSSDF_H

#define _BRSSDF_H

//蜂須賀先生のアイデア
// dirpole, directional dipole by T. Hachisuka and J. R. Frisvad
// originally smallpt, a path tracer by Kevin Beason, 2008
/*
LICENSE

Directional Dipole Model modifications:
Copyright (c) 2014 Toshiya Hachisuka and Jeppe Revall Frisvad

smallpt:
Copyright (c) 2006-2008 Kevin Beason (kevin.beason@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "intersection.h"
#include "random.h"

#include "entity.h"
#include "scene_env.h"
#include "Spectrum.h"

#include "phase.h"
#include "sampler.h"

namespace prender {

	inline double min3(const double x, const double y, const double z) {
		const double r = x < y ? x : y; return r < z ? r : z;
	}

	inline double C1(const double n) {
		double r;
		if (n > 1.0) {
			r = -9.23372 + n * (22.2272 + n * (-20.9292 + n * (10.2291 + n * (-2.54396 + 0.254913 * n))));
		}
		else {
			r = 0.919317 + n * (-3.4793 + n * (6.75335 + n *  (-7.80989 + n *(4.98554 - 1.36881 * n))));
		}
		return r / 2.0;
	}
	inline double C2(const double n) {
		double r = -1641.1 + n * (1213.67 + n * (-568.556 + n * (164.798 + n * (-27.0181 + 1.91826 * n))));
		r += (((135.926 / n) - 656.175) / n + 1376.53) / n;
		return r / 3.0;
	}

	//σ_tr
	inline Color sigma_tr(const Material& m, const Vector3d& D)
	{
		return Sqrt(m.participatingMediaParam.absorbingColor/D);
	}

	class DirectionalDipole
	{
	public:
		Vector3d sigma_s;
		Vector3d sigma_a;
		Vector3d sigma_t;
		Vector3d sigma_sp;
		Vector3d sigma_tp;
		Vector3d albedo_p;
		Vector3d D;
		Vector3d sigma_tr;
		Vector3d de;
		double Cp_norm;
		double Cp;
		double Ce;
		double A;
		double min_sigma_tr;
		double eta;

		DirectionalDipole(const Material& m, const double eta_)
		{
			eta = eta_;
			sigma_a = m.participatingMediaParam.absorbingColor;
			sigma_s = m.participatingMediaParam.scatteringColor;
			sigma_t = sigma_s + sigma_a;
			sigma_sp = sigma_s * (Vector3d(1.0, 1.0, 1.0) - Vector3d(m.participatingMediaParam.phase_prm));
			sigma_tp = sigma_sp + sigma_a;
			albedo_p = sigma_sp / sigma_tp;
			D = Vector3d(1.0, 1.0, 1.0) / (3.0 * sigma_tp);		//式(34)の下にある説明
			sigma_tr = Sqrt(sigma_a / D);
			de = 2.131 * D / Sqrt(albedo_p);		//式(21)
			
			Cp_norm = 1.0 / (1.0 - 2.0 * C1(1.0 / eta));

			//式(30),31)
			Cp = (1.0 - 2.0 * C1(eta)) / 4.0;
			Ce = (1.0 - 3.0 * C2(eta)) / 2.0;

			//式(22)
			A = (1.0 - Ce) / (2.0 * Cp);
			min_sigma_tr = sigma_tr.Min();

			//sigma_s.print("sigma_s");
			//sigma_a.print("sigma_a");
			//sigma_t.print("sigma_t");
			//sigma_sp.print("sigma_sp");
			//albedo_p.print("albedo_p");
			//D.print("D");
			//sigma_tr.print("sigma_tr");
			//de.print("de");
			//fprintf(stderr, "Cp_norm %f\n", Cp_norm);
			//fprintf(stderr, "Cp %f\n", Cp);
			//fprintf(stderr, "Ce %f\n", Ce);
			//fprintf(stderr, "A %f\n", A);
			//fprintf(stderr, "min_sigma_tr %f\n", min_sigma_tr);

		}

	};
	inline double Sp_d(const DirectionalDipole& directional_dipole, const Vector3d& x, const Vector3d& ω, const double r, const Vector3d& n0, const int j)
	{
		// evaluate the profile
		const double s_tr_r = directional_dipole.sigma_tr[j] * r;	//σ_tr
		const double s_tr_r_one = 1.0 + s_tr_r;				//(1 + σ_tr*r)
		const double x_dot_w = dot(x, ω);					//x・ω12
		const double r_sqr = r * r;

		//式(20)の計算

		//式(20)の先頭の係数
		// (1.0/Cφ(1/η))*(1/(4π^2))*exp(-σ_tr*r)/r^3
		double exp_tr = exp(-s_tr_r);
		//if (exp_tr < 1.0e-6) exp_tr = 1.0e-6;
		//const double t0 = directional_dipole.Cp_norm * (1.0 / (4.0 * PS_PI * PS_PI)) * exp_tr / (r * r_sqr);
		const double t0 = 1.0 / (4.0*directional_dipole.Cp) * (PS_INV_FORPI2)* exp_tr / (r * r_sqr);
		//fprintf(stderr, "directional_dipole.sigma_tr %f\n", directional_dipole.sigma_tr[j]);
		//fprintf(stderr, "%f %f\n", s_tr_r,  exp(-s_tr_r));
		//fprintf(stderr, "directional_dipole.Cp_norm %f t0 %.16f\n", directional_dipole.Cp_norm, t0);

		//式(20)の1項目
		// r^2/D + 3(1 + σ_tr*r) x・ω12
		const double t1 = r_sqr / directional_dipole.D[j] + 3.0 * s_tr_r_one * x_dot_w;

		//式(20)の2項目
		// 3D (1 + σ_tr*r)*ω12・n0 
		const double t2 = 3.0 * directional_dipole.D[j] * s_tr_r_one * dot(ω, n0);

		//式(20)の3項目
		//  (( 1 + σ_tr*r) + (3D(3( 1 + σ_tr*r) + σ_tr^2)/r^2) x・ω12)x・n0
		const double t3 = (s_tr_r_one + 3.0 * directional_dipole.D[j] * (3.0 * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * dot(x, n0);

		//fprintf(stderr, "t0 %f t1 %f t2 %f t3%f\n", t0, t1, t2, t3);
		//fprintf(stderr, "Cp %f Ce %f ==> %f\n", directional_dipole.Cp, directional_dipole.Ce, t0 * (directional_dipole.Cp * t1 - directional_dipole.Ce * (t2 - t3)));

		return t0 * (directional_dipole.Cp * t1 - directional_dipole.Ce * (t2 - t3));
	}

	inline Vector3d modified_tangent_plane_normal(const Vector3d& xi, const Vector3d& ni, const Vector3d& xo, const Vector3d& no)
	{
		// distance
		const Vector3d xoxi = xo - xi;
		const double r = xoxi.length();

		// modified normal
		//式(23)
		Vector3d ni_s = cross(normalize(xoxi), normalize(cross(ni, xoxi)));
		if (r < PS_EPS8) ni_s = ni;	//追加2015.3.11
		
		return ni_s;
	}

	inline double bssrdf(const DirectionalDipole& directional_dipole,  const Vector3d& xi, const Vector3d& ni, const Vector3d& wi, const Vector3d& xo, const Vector3d& no, const Vector3d& wo, const int j)
	{
		// distance
		const Vector3d xoxi = xo - xi;
		const double r = xoxi.length();

		// modified normal
		//式(23)
		Vector3d ni_s = cross(normalize(xoxi), normalize(cross(ni, xoxi)));
		if (r < PS_EPS8) ni_s = ni;	//追加2015.3.11

		// directions of ray sources
		const double nnt = 1.0 / directional_dipole.eta, ddn = -dot(wi, ni);
		const Vector3d ω12 = normalize(wi * -nnt - ni * (ddn * nnt + sqrt(1.0 - nnt * nnt * (1.0 - ddn * ddn))));
		const Vector3d ωv = ω12 - ni_s * (2.0 * dot(ω12, ni_s));

		// distance to real sources
		//式(25)
		const double cos_beta = -sqrt((r * r - dot(xoxi, ω12) * dot(xoxi, ω12)) / (r * r + directional_dipole.de[j] * directional_dipole.de[j]));

		//式(24)の計算
		double dr;
		const double μ0 = -dot(no, ω12);
		if (μ0 > 0.0) 
		{
			double tmp = (directional_dipole.D[j] * μ0) * ((directional_dipole.D[j] * μ0) - directional_dipole.de[j] * cos_beta * 2.0) + r * r;
			if ( tmp < 0 ) tmp = 0.0;
			dr = sqrt(tmp);
		}
		else 
		{
			double tmp = 1.0 / (3.0 * directional_dipole.sigma_t[j] * 3.0 * directional_dipole.sigma_t[j]) + r * r;
			if ( tmp < 0 ) tmp = 0.0;
			dr = sqrt(tmp);
		}

		// distance to virtual source
		const Vector3d xoxv = xo - (xi + ni_s * (2.0 * directional_dipole.A * directional_dipole.de[j]));
		const double dv = xoxv.length();

		// BSSRDF 式(26)
		//Sd(xo - xi, ω12, dr) - Sd(xo - xv, ωv, dv)
		const double result = Sp_d(directional_dipole, xoxi, ω12, dr, no, j) - Sp_d(directional_dipole, xoxv, ωv, dv, no, j);

		//fprintf(stderr, "%f - %f = %f\n", Sp_d(directional_dipole, xoxi, ω12, dr, no, j), Sp_d(directional_dipole, xoxv, ωv, dv, no, j), result);

		// clamping to zero
		return (result < 0.0) ? 0.0 : result;
	}
	// --------------------------------

	//eta = 外部屈折率/内部屈折率    η = η2/η1
	inline double fresnel(const double cos_theta, const double eta) 
	{
		const double sin_theta_t_sqr = 1.0 / (eta * eta) * (1.0 - cos_theta * cos_theta);
		if (sin_theta_t_sqr >= 1.0) return 1.0;
		const double cos_theta_t = sqrt(1.0 - sin_theta_t_sqr);
		const double r_s = (cos_theta - eta * cos_theta_t) / (cos_theta + eta * cos_theta_t);
		const double r_p = (eta * cos_theta - cos_theta_t) / (eta * cos_theta + cos_theta_t);
		return (r_s * r_s + r_p * r_p) * 0.5;
	}
};

#endif