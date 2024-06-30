#ifndef _RADIANCE_H_
#define _RADIANCE_H_

#include <algorithm>
#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "intersection.h"
#include "random.h"

#include "entity.h"
#include "scene_env.h"
#include "Spectrum.h"

#include "directLight_radiance.h"

#include "ParticipatingMedia.h"
#include "Subsurface_Scattering.h"
#include "brdf.h"

namespace prender {


// ray方向からの放射輝度を求める
Spectrum radiance(int* hit_direct_light, const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
{
	Intersection intersection;

	// シーンと交差判定
	if (!intersect_scene(ray, &intersection, depth))
	{
		if (intersection.object_id == EVENT_HORIZON)
		{
			//return RGB2Spectrum(Color(0.5, 0.5, 0.5), wavelength);
			return RGB2Spectrum(Color(0, 0, 0), wavelength);
		}
		if (intersection.object_id == BACK_GLOUND || intersection.object_id == BACK_GLOUND2)
		{
			return RGB2Spectrum(intersection.hitpoint.material.color, wavelength);
		}

		return RGB2Spectrum(env->BackgroundColor, 0.380);
	}

#if !IBL_LIGHT
	if (intersection.hitpoint.material.IBL() )
	{
		return RGB2Spectrum(intersection.hitpoint.material.emission, wavelength);
	}
#endif

	if ( depth == 0 && intersection.hitpoint.material.background )
	{
		return RGB2Spectrum(intersection.hitpoint.material.color, 0.380);
	}

	const int DepthLimit = env->DepthLimit;
	const int Depth      = env->Depth;

	const Material& hitpos_material = intersection.hitpoint.material;

	const Entity *now_object = EntList->List[intersection.object_id];
	const IntersectionPos &hitpoint = intersection.hitpoint;


	//if ( now_object->back == 0)
	//{
	//	if ( dot(hitpoint.normal , ray.dir) >= 0.0 )
	//	{
	//		return RGB2Spectrum(Color(0.0));
	//	}
	//}

	// 交差位置の法線（物体からのレイの入出を考慮）
	const Vector3d orienting_normal = dot(hitpoint.normal , ray.dir) < 0.0 ? hitpoint.normal: (-1.0 * hitpoint.normal); 

	// 色の反射率最大のものを得る。ロシアンルーレットで使う。
	// ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
	Color reflection_w = hitpos_material.emission;
	switch (hitpoint.material.reflection_type)
	{
		case REFLECTION_TYPE_DIFFUSE:		reflection_w = reflection_w + hitpoint.material.color;break;
		case REFLECTION_TYPE_SPECULAR:		reflection_w = reflection_w + hitpoint.material.color;break;
		case REFLECTION_TYPE_REFRACTION:	reflection_w = reflection_w + hitpoint.material.color;break;
		case REFLECTION_TYPE_WARD_BRDF:		reflection_w = reflection_w + hitpoint.material.color+hitpoint.material.specular;break;
		case REFLECTION_TYPE_PHONG_BRDF:	reflection_w = reflection_w + hitpoint.material.specular; break;
		default:			reflection_w = reflection_w + hitpoint.material.color; break;
	}
	double russian_roulette_probability = RGB_MAX(RGB2Spectrum(reflection_w, wavelength));


	// 反射回数が一定以上になったらロシアンルーレットの確率を急上昇させる。（スタックオーバーフロー対策）
	if (depth > DepthLimit)
	{
		russian_roulette_probability *= pow(0.5, depth - DepthLimit);
	}

	Color emission = ZERO();					//自己発光
	Spectrum direct_light = ZERO();				//直接光(反射)
	Spectrum direct_light_refraction = ZERO();	//直接光(屈折)

	emission = hitpos_material.emission;		// Init color to emissive value

	// ロシアンルーレットを実行し追跡を打ち切るかどうかを判断する。
	// ただしDepth回の追跡は保障する。
	if (depth > Depth)
	{
		if (rnd->next01() >= russian_roulette_probability)
		{
			Spectrum transmittance_ratio = ONE();

			//関与媒質を考慮する場合は光の減衰を考慮する
			if ( participatingMedia ) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance,wavelength);

			//return RGB2Spectrum( hitpos_material.emission, wavelength);
			return RGB2Spectrum(transmittance_ratio * hitpos_material.emission, wavelength);
			//if ( nextEventEstimation && depth == 0 || !nextEventEstimation )
			//{
			//	return RGB2Spectrum(MULTIPLY(transmittance_ratio, hitpos_material.emission), wavelength);
			//}
			//return ZERO();
		}
	} else
	{
		russian_roulette_probability = 1.0; // ロシアンルーレット実行しなかった
	}


	Spectrum incoming_radiance = 0.0;
	Spectrum weight = 1.0;

	// no absorption inside glass
	bool from_outside = now_object->isInsideDir2(ray, hitpoint);


	// レイと物体の交差点からの放射輝度を計算するか、途中の点における周囲からの影響を計算するかを選択するロシアンルーレット
	double scattering_probability = 0.0;
	
	//物体内部では関与媒質との干渉は無い
	if (participatingMedia && from_outside)
	{
		scattering_probability = 1.0- exp(-env->participatingMediaParam.scattering*hitpoint.distance);
		if ( scattering_probability < 0.01 ) scattering_probability = 0.01;
		else if ( scattering_probability > 0.99 ) scattering_probability = 0.99;

		if (rnd->next01() < scattering_probability)
		{
			return participatingMedia_radiance(env, rnd, ray, hitpoint, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}
	}

	switch (hitpos_material.reflection_type) 
	{
	case REFLECTION_TYPE_SSS_REFRACTION:
	case REFLECTION_TYPE_SSS_DIFFUSE:
	case REFLECTION_TYPE_SSS_WARD_BRDF:
	{
		emission = hitpos_material.emission;		// Init color to emissive value

		if (1 * hitpos_material.participatingMediaParam.use_brssdf && (hitpos_material.participatingMediaParam.level == -2 || hitpos_material.participatingMediaParam.level == -1))
		{
			if (hitpos_material.reflection_type == REFLECTION_TYPE_SSS_REFRACTION)
			{
				const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか
				// Snellの法則
				const double nc = 1.0;								// 真空の屈折率
				double nt = hitpos_material.refractive_index;		// オブジェクトの屈折率(kIor)

				nt = hitpos_material.refractive(wavelength);

				const double nnt = into ? nc / nt : nt / nc;
				const double ddn = dot(ray.dir, orienting_normal);
				const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

				// SchlickによるFresnelの反射係数の近似を使う
				const double a = nt - nc, b = nt + nc;
				const double R0 = (a * a) / (b * b);

				// 屈折の方向
				Ray tdir = Ray(hitpoint.position,
					normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(fabs(cos2t)))));
				const double c = 1.0 - (into ? -ddn : dot(tdir.dir, -1.0 * orienting_normal));
				const double Re = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
				const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
				double Tr = (1.0 - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

				if (cos2t < 0)
				{
					//IntersectionPos hitp = hitpoint;
					//hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION_FRESNEL;
					////hitp.material.color = participating_media.albedoColor*hitp.material.color;
					//return BRDF_Refraction(env, rnd, Ray(hitpoint.position, ray.dir), hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, 1 * nextEventEstimation, participatingMedia);
					
					tdir = Ray(hitpoint.position, normalize(ray.dir - orienting_normal * 2.0 * dot(orienting_normal, ray.dir)));
					Tr = Re;
				}

				const ParticipatingMedia& participating_media = const_cast<Entity*>(now_object)->material()->participatingMediaParam;
				Spectrum tr = ONE();
				if (!into)
				{
					//出て行く方向だから直前からの距離分だけ減衰している
					tr = participating_media.transmittanceRatio(hitpoint.distance, wavelength);
				}

				const double probability = 1;
				if (rnd->next01() < probability)
				{
					IntersectionPos hitp = hitpoint;
					//hitp.material.color =  participating_media.albedoColor;
					return tr*BRDF_SSS2(env, rnd, tdir, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, 1 * nextEventEstimation, participatingMedia, true) / probability;
				}

				IntersectionPos hitp = hitpoint;
				hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION;
				hitp.material.color = hitp.material.color;
				return tr*BRDF_Refraction(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia) / (1.0 - probability);
			}
		}
		else
		{
			Spectrum Li = ZERO();
			Spectrum weight = ONE();


			const ParticipatingMedia& participating_media = const_cast<Entity*>(now_object)->material()->participatingMediaParam;

			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか
			Spectrum tr = ONE();
			if (!into)
			{
				const ParticipatingMedia& participating_media = const_cast<Entity*>(now_object)->material()->participatingMediaParam;
				//出て行く方向だから直前からの距離分だけ減衰している
				tr = participating_media.transmittanceRatio(hitpoint.distance, wavelength);
			}


			// Snellの法則
			const double nc = 1.0;								// 真空の屈折率
			double nt = hitpos_material.refractive_index;		// オブジェクトの屈折率(kIor)

			nt = hitpos_material.refractive(wavelength);

			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = dot(ray.dir, orienting_normal);
			const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

			// SchlickによるFresnelの反射係数の近似を使う
			const double a = nt - nc, b = nt + nc;
			const double R0 = (a * a) / (b * b);

			// 屈折の方向
			Ray tdir = Ray(hitpoint.position,
				normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(fabs(cos2t)))));
			const double c = 1.0 - (into ? -ddn : dot(tdir.dir, -1.0 * orienting_normal));
			double Re = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
			const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
			double Tr = (1.0 - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

			double probability = 0.5;//std::min(0.5, std::max(1.0,participating_media.scattering/participating_media.transmittance)); //0.5;
			int sc_max = -1;
			sc_max = participating_media.level;
			if ( sc_max == 0 )
			{
				probability = 0;
				Re = 1.0;
			}

			if (!into )
			{
				//内側から来るのは散乱過程がキャンセルされた場合なのでそのまま屈折して出ていく
				IntersectionPos hitp = hitpoint;
				hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION_FRESNEL;
				//hitp.material.color = participating_media.albedoColor*hitp.material.color;
				return tr*BRDF_Refraction(env, rnd, Ray(hitpoint.position, ray.dir), hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, 1 * nextEventEstimation, participatingMedia);
			}

			if (cos2t < 0)
			{
				probability = 0;
				//Re = 1.0;
				//IntersectionPos hitp = hitpoint;
				//hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION_FRESNEL;
				//hitp.material.color = participating_media.albedoColor*hitp.material.color;
				//return tr*BRDF_Refraction(env, rnd, Ray(hitpoint.position, ray.dir), hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, 1 * nextEventEstimation, participatingMedia);
			}

			if (rnd->next01() < probability)
			{
				//散乱過程に入るので物体の内側に少しだけ入っておく
				if (cos2t < 0)
				{
					tdir = Ray(hitpoint.position - 0.0001*hitpoint.normal, ray.dir);
				}
				else
				{
					tdir.org = hitpoint.position - 0.0001*hitpoint.normal;
				}
				const double T21 = 1.0;// std::max(1.0 - fresnel(dot(tdir.dir*-1.0, hitpoint.normal), hitpoint.material.refractive_index), 0.0);

				double totalDist = 0.0;
				Tr = 1.0;
				Li = Tr*radiance_through_media(env, hitpoint, tdir, rnd, depth + 1, now_object, 0, sc_max, russian_roulette_probability, scattering_probability, wavelength, totalDist, 1 * nextEventEstimation, participatingMedia);

				//Li = RGB2Spectrum(participating_media.albedoColor, wavelength)*Li;

				weight = tr*RGB2Spectrum(hitpoint.material.color, wavelength)
					/ russian_roulette_probability / (1.0 - scattering_probability) / probability;

				//ここでnextEventEstimationをやってはダメ（出口でnextEventEstimationをしている）
				if (0)
				{
					IntersectionPos hitp = hitpoint;
					hitp.material.color = participating_media.albedoColor;
					if (nextEventEstimation && env->light_list.list.size())
					{
						// 直接光の評価
						//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
						double light_probability = 1.0;
						for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
						{
							const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);

							//weightが係って返ってくる
							Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitp, orienting_normal, ray.dir, depth, wavelength) / light_probability;

							direct_light = direct_light + dl;
						}
						direct_light = direct_light / SHADOW_RAY_SAMPLING / russian_roulette_probability / (1.0 - scattering_probability) / probability;
						;
					}
					// 直接光の評価
					if (!nextEventEstimation && env->light_list.list.size())
					{
						// 直接光の評価
						//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
						double light_probability = 1.0;
						for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
						{
							const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
							if (light.parallel_light || light.spot_light || light.infinity_light)
							{

								//weightが係って返ってくる
								Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitp, orienting_normal, ray.dir, depth, wavelength) / light_probability;

								direct_light = direct_light + dl;
							}
						}
						direct_light = direct_light / SHADOW_RAY_SAMPLING / russian_roulette_probability / (1.0 - scattering_probability) / probability;
						;
					}
				}
			}
			else
			{
				if (hitpoint.material.roughness)
				{
					double diffuse_probability = hitpoint.material.roughness;

					double η = 1.495;//polypropylene;
					double m_invη2 = 1.0 / (η*η);

					//●●●
					//Ks + Kd <= 1.0
					// weight = Ks * weight_s + Kd * weight_d = 
					//        = (hitpoint.material.roughness-1) * weight_s + hitpoint.material.roughness * weight_d
					// ロシアンルーレットでスペキュラ―とディフュースを分けて計算する
					// ディフュースのpdf は　diffuse_probabilityなので　ディフュース=hitpoint.material.roughness * weight_d/diffuse_probability
					// diffuse_probability == hitpoint.material.roughness なのでpdfと相殺して　ディフュース=weight_d
					// 同様にスペキュラ―=(hitpoint.material.roughness-1) * weight_s
					//スペキュラ―のpdf=(hitpoint.material.roughness-1)なのでpdfと相殺してスペキュラ―=weight_s
					if (rnd->next01() < diffuse_probability)
					{
						IntersectionPos hitp = hitpoint;

						hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
						hitp.material.color = hitp.material.color*participating_media.albedoColor;

						weight = tr*m_invη2*(1.0 - fresnelDiffuseReflectance(η)) / (1.0 - probability);

						Li = BRDF_Diffuse(0, env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
					}
					else
					{
						IntersectionPos hitp = hitpoint;
						hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION;
						hitp.material.refractive_index = η;
						hitp.material.color = hitp.material.color;

						weight = tr / (1.0 - probability);

						Li = BRDF_Refraction(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
					}
				}
				else
				{
					IntersectionPos hitp = hitpoint;
					hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION;

					weight = tr / (1.0 - probability);
					Li = BRDF_Refraction(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
				}
			}

			incoming_radiance = Li;

			//関与媒質を考慮する場合は光の減衰を考慮する
			if (participatingMedia)
			{
				// no absorption inside glass
				bool inside_dir = const_cast<Entity*>(now_object)->isInsideDir2(ray, hitpoint);
				
				//関与媒質を考慮する場合は光の減衰を考慮する
				Spectrum transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

				//環境内（物体の外側)
				//現在の交点ではレイは物体の中に入ろうとする方向なら物体の外側から
				if (inside_dir)
				{
					direct_light = transmittance_ratio*direct_light;
					emission = transmittance_ratio * emission;
					incoming_radiance = transmittance_ratio * incoming_radiance;
				}
			}

			return RGB2Spectrum(env->Environment_Light, wavelength) + direct_light + incoming_radiance*weight;
		}
	}
	break;
	case REFLECTION_TYPE_PHONG_BRDF:
	{
#if 0
		const double rd = (hitpos_material.color.x + hitpos_material.color.y + hitpos_material.color.z) / 3.0;
		const double rs = (hitpos_material.specular.x + hitpos_material.specular.y + hitpos_material.specular.z) / 3.0;
		const double diffuse_probability = rd / (rs + rd);
		const double specular_probability = 1.0 - diffuse_probability;
		const double u = rnd->next01();

		double costh = dot(orienting_normal, ray.dir);
		if (costh < 0) costh = -1.0*costh;

		const double temp1 = 1.0 - costh;
		const double R0 = hitpos_material.roughness;
		double Re = R0 + (1.0 - R0) * temp1 * temp1 * temp1 * temp1 * temp1;
		double Tr = 1.0 - Re;
		double P = (Re + 0.5) / 2.0;

		if (diffuse_probability == 0.0 || hitpos_material.roughness == 0.0)
		{
			Re = 1.0;
			Tr = 0.0;
		}
		if (u < P)
		{
			IntersectionPos hitp = hitpoint;

			hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
			hitp.material.color = (hitp.material.color) * Tr / P;

			return BRDF_Diffuse(hit_direct_light, env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}
		IntersectionPos hitp = hitpoint;
		hitp.material.specular = (hitp.material.specular) * Re / ( 1.0 - P);
		return BRDF_Phong(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
#else
		//double costh = dot(orienting_normal, ray.dir);
		//if (costh < 0) costh = -1.0*costh;

		//const double temp1 = 1.0 - costh;
		//const double R0 = hitpos_material.roughness;
		//double Re = R0 + (1.0 - R0) * temp1 * temp1 * temp1 * temp1 * temp1;
		//double Tr = 1.0 - Re;
		//double P = (Re + 0.5) / 2.0;

		double diffuse_probability = hitpos_material.roughness;
		if (rnd->next01() < diffuse_probability)
		{
			IntersectionPos hitp = hitpoint;
			hitp.material.color = hitp.material.color;
			hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
			return BRDF_Diffuse(hit_direct_light, env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}
		else
		{
			IntersectionPos hitp = hitpoint;
			hitp.material.color = hitp.material.color;
			return BRDF_Phong(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}
#endif
	} break;
	// Ward BRDF
	case REFLECTION_TYPE_WARD_BRDF: 
	{
		const double rd = (hitpos_material.color.x + hitpos_material.color.y + hitpos_material.color.z)/3.0;
		const double rs = (hitpos_material.specular.x + hitpos_material.specular.y + hitpos_material.specular.z) / 3.0;
		const double diffuse_probability = rd / (rs + rd);
		const double specular_probability = 1.0 - diffuse_probability;
		const double u = rnd->next01();
		
		if (u < diffuse_probability)
		{
			IntersectionPos hitp = hitpoint;

			hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
			hitp.material.color = hitp.material.color;
			return BRDF_Diffuse(hit_direct_light, env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}

		IntersectionPos hitp = hitpoint;
		hitp.material.specular = hitp.material.specular;
		return BRDF_Ward(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
	} break;
	case REFLECTION_TYPE_DIFFUSE: 
	{
		const double roughness = hitpoint.material.roughness;
		if (rnd->next01() < roughness)
		{
			IntersectionPos hitp = hitpoint;

			hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
			hitp.material.color = hitp.material.color;
			return BRDF_Diffuse(hit_direct_light, env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}
		IntersectionPos hitp = hitpoint;
		hitp.material.reflection_type = REFLECTION_TYPE_SPECULAR;
		hitp.material.color = hitp.material.color;
		return BRDF_Specular(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
	}break;
	// 完全鏡面
	case REFLECTION_TYPE_SPECULAR:
	{
		const double specular = 1.0 - hitpoint.material.roughness;
		if (rnd->next01() < specular)
		{
			IntersectionPos hitp = hitpoint;

			hitp.material.specular = hitp.material.specular;
			return BRDF_Specular(env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
		}
		IntersectionPos hitp = hitpoint;
		hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
		hitp.material.color = hitp.material.color;
		return BRDF_Diffuse(hit_direct_light, env, rnd, ray, hitp, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
	} break;

	// 屈折
	case REFLECTION_TYPE_REFRACTION:
	case REFLECTION_TYPE_REFRACTION_FRESNEL:{
		return BRDF_Refraction(env, rnd, ray, hitpoint, now_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
	} break;

	}

	return ZERO();
}

};

#endif
