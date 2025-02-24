#ifndef _PHASE_H_
#define _PHASE_H_

#include "scene_env.h"

namespace prender {

	//Isotropic位相関数( 等方散乱)
	inline double PhaseFunction_Isotropic()
	{
		return PS_INV_FORPI;
	}
	//Henyey-Greenstein関数
	inline double PhaseFunction_HenyeyGreenstein(const SceneEnv* env, const double g, const double costheta)
	{
		const double w = 1.0 + g*g - 2.0*g*costheta;
		if (fabs(w) < PS_EPS12) return 0.0;

		return (PS_INV_FORPI)* (1.0 - g*g) / pow(w, 3.0 / 2.0);
	}

	//Schlick位相関数(Henyey-Greenstein位相関数を近似)
	inline double PhaseFunction_Schlick(const SceneEnv* env, const double g, const double costheta)
	{
		const double k = 1.55*g - 0.55*pow(g, 3.0);
		const double w = 1.0 - k*costheta;
		if (fabs(w) < PS_EPS12) return 0.0;

		return (PS_INV_FORPI)* (1.0 - g*g) / pow(w, 2.0);
	}

	//Rayleigh位相関数
	inline double PhaseFunction_Rayleigh(const SceneEnv* env, const double g, const double costheta, const double wavelength)
	{
		return (3.0 / (4.0*PS_PI))* (1.0 + costheta*costheta)/ pow(wavelength, 4.0);
	}
	//Lorenz-Mie位相関数
	inline double PhaseFunction_Lorenz_Mie_murky(const SceneEnv* env, const double g, const double costheta)
	{
		return (PS_INV_FORPI)* (1.0 / 2.0 + (33.0 / 2.0)*pow((1.0 + costheta) / 2.0, 32));
	}
	//Lorenz-Mie位相関数
	inline double PhaseFunction_Lorenz_Mie_hazy(const SceneEnv* env, const double g, const double costheta)
	{
		return (PS_INV_FORPI)* (1.0 / 2.0 + (9.0 / 2.0)*pow((1.0 + costheta) / 2.0, 8));
	}
};

#endif
