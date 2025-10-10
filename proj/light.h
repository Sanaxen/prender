#ifndef __LIGHT_H_
#define __LIGHT_H_

#include <vector>
#include "entity.h"
#include "random.h"

#include <algorithm>
#include <functional>

namespace prender {

class Light
{
public:

	double power;
	Vector3d dir;		//光源方向ベクトル

	int parallel_light;	//平行光源
	int infinity_light;
	Color infinity_light_emitssin;

	//スポットライト
	int spot_light;
	double spot_light_cos_angle[2];	//[0]:cosTotalWidth   [1]:cosFalloffStart
	double falloff;
	double attenuation[3];

	Entity* light;
	int light_id;	//LightList.listのインデックス

	Light()
	{
		dir = Vector3d(0, 0, 0);

		parallel_light = 0;
		infinity_light = 0;

		spot_light = 0;
		falloff = 2.0;
		spot_light_cos_angle[0] = cos(30.0*PS_RAD);
		spot_light_cos_angle[1] = cos(25.0*PS_RAD);	//25.0 = 30.0-5.0
		attenuation[0] = 1.0;
		attenuation[1] = 0.0;
		attenuation[2] = 0.0;
	}
	inline double SpotLightPower(const double Intensity) const
	{
		return Intensity * 2.0 * PS_PI *(1.0 - .5 * (spot_light_cos_angle[1] + spot_light_cos_angle[0]));
	}
	inline double SpotLightFalloff(const Vector3d &w) const
	{
		double cosTotalWidth = spot_light_cos_angle[0];
		double cosFalloffStart = spot_light_cos_angle[1];

		double costheta = dot(dir, w);
		if (costheta < cosTotalWidth)     return 0.0;
		if (costheta > cosFalloffStart)   return 1.0;
		// Compute falloff inside spotlight cone
		double delta = (costheta - cosTotalWidth) / (cosFalloffStart - cosTotalWidth);
		return pow(delta, falloff);
	}
	inline Color SpotLightEmission(const Color& Intensity, const Vector3d& lightPos, const Vector3d& p) const
	{
		Vector3d wi = normalize(lightPos - p);
		const double dist = (lightPos - p).length();

		double invA = attenuation[0] + dist * (attenuation[1] + attenuation[2] * dist);
		if (invA < 1.0e-3) invA = 1.0;

		return Intensity * SpotLightFalloff(wi*-1.0) / invA;
	}
};

inline bool operator<(const Light& left, const Light& right)
{
	return left.power < right.power;
}

inline bool operator>(const Light& left, const Light& right)
{
	return left.power > right.power;
}

class LightList
{
	std::vector<double> eachLightProbability;
	std::vector<double> accumulatedProbability;
	double totalPower;
public:
	std::vector<Light> list;


	void  random_light_setup()
	{
		std::sort(list.begin(), list.end(), std::greater<Light>());//降順ソート

		totalPower = 0.0;
		MTRandom r;
		for ( int i = 0; i < list.size(); i++ )
		{
			//randomLightで選ぶとき同じpowerがあると最初に見つかった光源しか選ばれないのでPowerにcを加算しておく
			double c = r.next01()*10.0;

			eachLightProbability.push_back(list[i].power+c);
			accumulatedProbability.push_back(totalPower+list[i].power+c);
			totalPower += list[i].power+c;
		}

		// 確率へ正規化
		for (size_t i = 0; i < list.size(); i++) 
		{
			eachLightProbability[i] /= totalPower;
			accumulatedProbability[i] /= totalPower;
		}
		for (size_t i = 0; i < list.size(); i++) 
		{
			if ( list[i].infinity_light )
			{
				printf("inifinity light:");
			}else
			if (list[i].light->material()->IBL())
			{
				printf("IBL:");
			}
			printf("light %d power %f cdf %f pdf %f\n", (int)i, list[i].power, accumulatedProbability[i], eachLightProbability[i]);
		}

	}

	inline Light& randomLight(Random *rnd, double& probability)
	{
#if 0
		const int sz = list.size();
		double next = rnd->_next01();
		for (size_t i = 0; i < sz; i++) 
		{
			if (next <= accumulatedProbability[i]) 
			{
				probability = eachLightProbability[i];
				return list[i];
			}
		}
#else
		probability = 1.0;
		double next = rnd->Next01();
		int index = (int)((next)*(double)(list.size()-1) + 0.5);
		if ( index >= list.size() ) index = list.size()-1;
		if ( index < 0 ) index = 0;
		return list[index];

#endif
		return list[0];
	}

	inline Spectrum InfinityLight(const Vector3d& dir, int& exist_infinity_light, const double wavelength)
	{
		exist_infinity_light = 0;
		for ( int i = 0; i < list.size(); i++ )
		{
			if ( list[i].infinity_light )
			{
				if ( fabs(dot(list[i].dir, dir*-1.0 ) -1.0) < 1.0e-3 )
				{
					exist_infinity_light++;
					return RGB2Spectrum(list[i].infinity_light_emitssin, wavelength);
				}
			}
		}
		return ZERO();
	}
};


};

#endif
