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

namespace prender {
Spectrum radiance(const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);

Color radiance_through_media(const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const Entity* in_object, const int sc_depth, const double wavelength, const int nextEventEstimation=0, const int participatingMedia=0) 
{
	const double eps = 0.0001;

	/*
	const double scattering = 0.999;
	const double absorbing = 0.001;
	const double transmittance = scattering + absorbing;
	*/

//		const Vec scattering(0.9, 0.999, 0.9);
	//const Vector3d scattering(0.7, 0.7, 0.7);
	//const Vector3d absorbing(0.9, 0.001, 0.9);
	const ParticipatingMedia& participating_media = const_cast<Entity*>(in_object)->material()->participatingMediaParam;
	const double phase = 1.0 / (4.0 * PS_PI);


	double russian_roulette_probability = std::min(0.9, participating_media.scattering_albedo);
	if (depth > env->DepthLimit*4) 
	{
		if (rnd->next01() >= russian_roulette_probability)
		{
			//fprintf(stderr, "DepthLimit!!\n");
			return const_cast<Entity*>(in_object)->material()->emission;
		}
	} else
		russian_roulette_probability = 1.0; // ロシアンルーレット実行しなかった
	
	//出口の交点	rayは表面上の点から中に入っていく方向のレイ
	Intersection out;
	if (!intersect_scene(ray, &out)) 
	{
		fprintf(stderr, "@Impossible!!\n");
		// 絶対ここには入らないはず
		return env->BackgroundColor;
	}
	Entity* now_object = env->EntList.List[out.object_id];

	//中に入らなかった（おそらく接線方向でかすめて抜けて行って他の物体に当たった?)
	if ( now_object != in_object )
	{
		return ZERO();	//この寄与は無いものとする

		return radiance(env, ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)
			/ russian_roulette_probability;
	}

	//出口の交点
	const Vector3d hitpoint = out.hitpoint.position;
	//const Vector3d hitpoint = ray.org + out.hitpoint.distance * ray.dir;

	//{
	//	Vector3d dd = hitpoint-((Sphere*)in_object)->position;

	//	if ( fabs(sqrt(dot(dd,dd)) - ((Sphere*)in_object)->radius) > 0.0001 )
	//	{
	//		fprintf(stderr, "@@@@@@@@@@@@ hitpoint が表面に無い!!\n");
	//	}
	//}

	bool insidedir = now_object->isInsideDir(ray, out.hitpoint);


	//出口の交点では、、、
	if ( 0 )
	{
		if ( now_object == in_object )
		{
			if ( insidedir )
			{
				//表からレイが当たっている
				fprintf(stderr, "*中に入って行く\n");
			}else
			{
				//裏からレイが当たっている（内部だが他の表面だと間違い...)
				//fprintf(stderr, "*外に出て行く=>OK\n");
			}
		}else
		{
			if ( insidedir )
			{
				fprintf(stderr, "*他の中に入って行く\n");
			}else
			{
				fprintf(stderr, "*他の外に出て行く\n");
			}
		}
	}

	//ほとんど表面なのでそのまま抜ける
	if (out.hitpoint.distance  <= eps)
	{
		const Spectrum Tr = participating_media.transmittanceRatio(out.hitpoint.distance,wavelength);
		return Tr * radiance(env, ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia) / russian_roulette_probability;
	}


	Color Li, L;
	
	//double probability = std::min(0.1, sc_average);
	double probability = std::max(0.001, std::min(0.99, participating_media.scattering));
	//if (probability < 0.01) probability = 0.01;
	//const double probability = 1.0 - exp(-sc_average*out.hitpoint.distance);


	//表面位置から中に入って行く(ray方向)
	if (rnd->next01() < probability )
	{
		//散乱位置への移動
		double d = 0.0;
		do {
			d = -log(rnd->next01()) / participating_media.transmittance;
		} while (d >= out.hitpoint.distance*(1.0-eps));

		// 等方散乱
		const double r1 = 2 * PS_PI * rnd->next01();
		const double r2 = 1.0 - 2.0 * rnd->next01();

		//あたらなレイを飛ばす位置（散乱位置）
		Vector3d new_org = ray.org + d * ray.dir;
		Ray next_ray(new_org, normalize(Vector3d(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2)));
	
		//物体が球ならnew_orgはその球の中にあるはず
		//{
		//	
		//	if (((Sphere*)in_object)->isIn(new_org) != 1 )
		//	{
		//		fprintf(stderr, "外になった\n");
		//	}
		//}

		if (1)
		{
			Intersection intersectionWrk;
			if (!intersect_scene(next_ray, &intersectionWrk))
			{
				fprintf(stderr, "-Impossible!!\n");
				// 絶対ここには入らないはず
				return env->BackgroundColor;
			}
			Entity* now_object = env->EntList.List[intersectionWrk.object_id];
			bool insidedir = now_object->isInsideDir(next_ray, intersectionWrk.hitpoint);

			insidedir = now_object->isInsideDir(next_ray, intersectionWrk.hitpoint);

			if (now_object == in_object)
			{
				if (insidedir)
				{

					//散乱の結果、レイは表面にあたる（これはあり得ないはず)
					fprintf(stderr, "-Impossible( ->) ->()\n");
					return ZERO();

					return radiance(env, next_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)*0.5
						/ russian_roulette_probability;

				}
				else
				{
					//散乱の結果、レイは裏面にあたる
					//fprintf(stderr, "-外に出て行く\n");
				}
			}
			else
			{
				//中からレイを飛ばしているのに他の物体に当たった（ありえないはず、、、)
				//fprintf(stderr, "-Impossible( ->) ->{} or ->{->}%f (d=%f)\n", out.hitpoint.distance, d);
				return ZERO();

				///fprintf(stderr, "%f %d (%f,%f,%f)\n", dot(next_ray.dir, intersectionWrk.hitpoint.normal), insidedir ? 1 : 0, intersectionWrk.hitpoint.normal.x, intersectionWrk.hitpoint.normal.y, intersectionWrk.hitpoint.normal.z);

				const Spectrum Tr = participating_media.transmittanceRatio(d,wavelength);
				return Tr * radiance(env, next_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)
					/ russian_roulette_probability;
			}
		}


		//fprintf(stderr, "/内部で減衰して散乱\n");
		//内部で減衰して散乱(in-scatting)
		const Spectrum Tr = participating_media.transmittanceRatio(d,wavelength);

		Li = Tr * participating_media.scatteringColor * radiance_through_media(env, next_ray, rnd, depth+1, in_object, sc_depth+1, wavelength, nextEventEstimation, participatingMedia)
			* phase
			/ exp(-participating_media.transmittance * d)
			/ (1.0 / (4.0 * PS_PI)) 
			/ probability
			/ russian_roulette_probability;

		const Spectrum Trr = participating_media.transmittanceRatio(out.hitpoint.distance,wavelength);
		L = Trr * radiance(env, next_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)
			/ (1.0 - probability)
			/ russian_roulette_probability;

		Li = L + Li;
		//if (participating_media.absorbingColor.length() > 0 && participating_media.scatteringColor.length())
		//{
		//	//const Spectrum a = participating_media.absorbingRatio(d, wavelength);

		//	//Li = MULTIPLY((ONE() - a)*participating_media.absorbingColor, ONE() + Li);
		//	Li = MULTIPLY(participating_media.scatteringColor, Li);
		//}
		return Li;
	} else {
		//散乱しなかった場合

		//出口から飛ばす
		Ray next_ray(hitpoint + ray.dir*eps, ray.dir);


		if(1)
		{
			Intersection intersectionWrk;
			if (!intersect_scene(next_ray, &intersectionWrk)) 
			{
				fprintf(stderr, "*Impossible!!\n");
				// 絶対ここには入らないはず
				return env->BackgroundColor;
			}
			Entity* now_object = env->EntList.List[intersectionWrk.object_id];
			bool insidedir = now_object->isInsideDir(next_ray, intersectionWrk.hitpoint);

			//出口から出ていくレイなので再度自分に入って行くなんてありえない
			if ( now_object == in_object )
			{
				fprintf(stderr, "@Impossible ( )->  ->( ) or ( ->)%d\n", insidedir?1:0);
			}else
			{
				if ( insidedir )
				{
					//そのまま進んだ結果、レイは他の表面にあたる
					//fprintf(stderr, "@他の中に入って行く\n");
				}else
				{

					//そのまま進んだ結果、レイは他の裏面にあたる（これはあり得ないはず)
					//fprintf(stderr, "@Impossible ( )->  { ->}\n");

					//よくわからんが以下をスルーしないと駄目(2014.09.25)
					//return ZERO();

					//fprintf(stderr, "%f %d (%f,%f,%f)\n", dot(next_ray.dir, intersectionWrk.hitpoint.normal), insidedir ? 1 : 0, intersectionWrk.hitpoint.normal.x, intersectionWrk.hitpoint.normal.y, intersectionWrk.hitpoint.normal.z);

					//const double Tr = exp(-participating_media.transmittance * out.hitpoint.distance);
					//
					//return Tr * radiance(env, next_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)
					//	/ russian_roulette_probability;
				}			
			}
		}

		//内部で減衰してそのまま抜ける
		//fprintf(stderr, "内部で減衰してそのまま抜ける\n");
		const Spectrum Tr = participating_media.transmittanceRatio(out.hitpoint.distance,wavelength);
		L = Tr * radiance(env, next_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)
			/ (1.0 - probability)
			/ russian_roulette_probability;


		//if (participating_media.absorbingColor.length() > 0 )
		//{
		//	L = MULTIPLY(participating_media.absorbingColor, L);
		//	//L = MULTIPLY((ONE() - Tr)*participating_media.absorbingColor, ONE() + L);
		//}
		return L;
	}
}


// ray方向からの放射輝度を求める
Spectrum radiance(const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const double wavelength, const int nextEventEstimation=0, const int participatingMedia=0) 
{
	Intersection intersection;

	// シーンと交差判定
	if (!intersect_scene(ray, &intersection))
	{
		return RGB2Spectrum(env->BackgroundColor, 0.380);
	}
	//if ( intersection.hitpoint.material.IBL() )
	//{
	//	return intersection.hitpoint.material.color;
	//}

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
	switch( const_cast<Entity*>(now_object)->material()->reflection_type )
	{
		case REFLECTION_TYPE_DIFFUSE:		reflection_w = reflection_w+hitpoint.material.color;break;
		case REFLECTION_TYPE_SPECULAR:		reflection_w = reflection_w+hitpoint.material.color;break;
		case REFLECTION_TYPE_REFRACTION:	reflection_w = reflection_w+hitpoint.material.color;break;
		case REFLECTION_TYPE_WARD_BRD:		reflection_w = reflection_w+hitpoint.material.color+hitpoint.material.ward_brdf.specular;break;
		default:			reflection_w = reflection_w+hitpoint.material.color;break;
	}
	double russian_roulette_probability = RGB_MAX(RGB2Spectrum(reflection_w, wavelength));

	// 反射回数が一定以上になったらロシアンルーレットの確率を急上昇させる。（スタックオーバーフロー対策）
	if (depth > DepthLimit)
		russian_roulette_probability *= pow(0.5, depth - DepthLimit);

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
	bool inside_dir = now_object->isInsideDir(ray, hitpoint);


	// レイと物体の交差点からの放射輝度を計算するか、途中の点における周囲からの影響を計算するかを選択するロシアンルーレット
	double scattering_probability = 0.0;
	
	if ( participatingMedia )
	{
		//scattering_probability = std::min(0.9, env->participatingMediaParam.scattering);
		//scattering_probability = 1.0 - exp(-env->participatingMediaParam.scattering*hitpoint.distance);
		//if (scattering_probability < 1.0e-8)
		//{
		//	scattering_probability = 1.0e-8;
		//}
		//test
		
		//scattering_probability = 0.5;
		scattering_probability = std::min(0.99, exp(-env->participatingMediaParam.scattering_albedo));		

		if (rnd->next01() < scattering_probability) 
		{   
#if 10
			// 周囲からの影響を考慮（つまりまわりから来た光が散乱してray方向に来た）
			// 減衰に応じた重点サンプリング
#if 10
			const double u = rnd->next01();
			const double d = -log(1.0 + u * (exp(-env->participatingMediaParam.transmittance * hitpoint.distance) - 1.0)) / env->participatingMediaParam.transmittance;
			const double pdf = exp(-env->participatingMediaParam.transmittance * d) * (-env->participatingMediaParam.transmittance / (exp(-env->participatingMediaParam.transmittance * hitpoint.distance) - 1.0));
#else		
			const double eps = 0.0001;
			double d = 0.0;
			do {
				d = -log(rnd->next01()) / env->participatingMediaParam.transmittance;
			} while (d >= hitpoint.distance*(1.0-eps));
			const double pdf = exp(-env->participatingMediaParam.transmittance * d);
#endif
			// 等方散乱（始点から距離dの位置から方向が変わる)
			const double r1 = 2 * PS_PI * rnd->next01();
			const double r2 = 1.0 - 2.0 * rnd->next01() ;
			Vector3d next_dir = normalize(Vector3d(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2));

			//レイ視点と衝突点の途中位置からnext_dirへのレイ
			const Ray next_ray = Ray(ray.org + d * ray.dir, next_dir);
			const Spectrum transmittance_ratio = env->participatingMediaParam.transmittanceRatio(d,wavelength);
		
			IntersectionPos hitpoint_media = hitpoint;
			hitpoint_media.distance = d;
			hitpoint_media.position = next_ray.org;

			if ( !now_object->light )
			{
				double light_probability = 1.0;
				if ( env->light_list.list.size() )
				{
					for ( int i = 0; i < SHADOW_RAY_SAMPLING; i++ )
					{	
						const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
						direct_light = direct_light + direct_radiance_sample_media(env, rnd, light, hitpoint_media, depth, wavelength) /  light_probability;
					}
					direct_light = direct_light / SHADOW_RAY_SAMPLING;
				}
			}

			if(0)
			{
				Intersection intersectionWrk;
				if (!intersect_scene(next_ray, &intersectionWrk))
				{
					//fprintf(stderr, "-Impossible!!\n");
					// 絶対ここには入らないはず
					return env->BackgroundColor;
				}
				Entity* object = env->EntList.List[intersectionWrk.object_id];
				bool insidedir = object->isInsideDir(next_ray, intersectionWrk.hitpoint);

				insidedir = object->isInsideDir(next_ray, intersectionWrk.hitpoint);
				if ( !insidedir )
				{
					return ZERO();
				}
			}

			if (pdf == 0.0) {
				return ZERO();
			} else {

				//const double phase = PhaseFunction_HenyeyGreenstein(env, dot(next_dir, ray.dir));
				const double phase = 1.0 / (4.0 * PS_PI); // 等方散乱

				Spectrum incoming_radiance = ZERO();
				
				if ( now_object->light )
				{
					if ( depth == 0 )
					{
						incoming_radiance = const_cast<Entity*>(now_object)->material()->emission;
						return transmittance_ratio*incoming_radiance;
					}else
					{
						return ZERO();
					}
				}
				Spectrum L = radiance(env, next_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
				
				L =  RGB2Spectrum(
					(transmittance_ratio * env->participatingMediaParam.scattering *( direct_light+L) )
					* phase/ pdf
					/ (1.0 / (4.0 * PS_PI))
					/ scattering_probability
					/ russian_roulette_probability,
					wavelength);

				//const Spectrum Tr = env->participatingMediaParam.transmittanceRatio(hitpoint.distance,wavelength);
				//L = L + Tr * radiance(env, next_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)
				//	/ scattering_probability
				//	/ russian_roulette_probability;

				//if (env->participatingMediaParam.absorbingColor.length() > 0 && env->participatingMediaParam.scatteringColor.length())
				//{
				//	const double a = env->participatingMediaParam.absorbingRatio(d);
				//	L = MULTIPLY((1.0 - a)*env->participatingMediaParam.scatteringColor, ONE() + L);
				//}
				return L;
			}
#endif
			return ZERO();
		}
	}

	Spectrum transmittance_ratio = ONE();
	//関与媒質を考慮する場合は光の減衰を考慮する
	if ( participatingMedia ) transmittance_ratio  = env->participatingMediaParam.transmittanceRatio(hitpoint.distance,wavelength);


	switch (hitpos_material.reflection_type) 
	{
	case REFLECTION_TYPE_SSS_REFRACTION:
	case REFLECTION_TYPE_SSS_DIFFUSE:
	case REFLECTION_TYPE_SSS_WARD_BRD:
		{
		 int direct_light_meshod = 1;
		
		 Vector3d tdir = ray.dir;

		// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
		Vector3d w, u, v;
		OrthonormalBasis(orienting_normal, w, u, v);
		

		// コサイン項を使った重点的サンプリング
		const double r1 = 2 * PS_PI * rnd->next01();
		const double r2 = rnd->next01(), r2s = sqrt(r2);
		Vector3d dir = normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));


		if (hitpos_material.reflection_type==REFLECTION_TYPE_SSS_WARD_BRD)
		{
			const double alpha_x = hitpos_material.ward_brdf.alp_x;
			const double alpha_y = hitpos_material.ward_brdf.alp_y;
			const Vector3d in = -1.0 * ray.dir;
			Vector3d halfv;
			do {
				const double u1 = rnd->next01();
				const double u2 = rnd->next01();
				double phi = atan(alpha_x / alpha_y * tan(PS_TWOPI * u2));
				if (0.25 <= u2 && u2 <= 0.75)
					phi += PS_PI;
				else if (0.75 < u2)
					phi += PS_TWOPI;
				const double theta = atan(sqrt(-log(u1) / (pow(cos(phi), 2) / pow(alpha_x, 2) + pow(sin(phi), 2) / pow(alpha_y, 2))));

				halfv = normalize((
					u * cos(phi) * sin(theta) + 
					v * sin(phi) * sin(theta) + 
					w * cos(theta)));

				dir = 2.0 * dot(in, halfv) * halfv - in;
			} while (dot(orienting_normal, dir) < 0.0);
		}

		double Re = 0.5;
		double Tr = 0.5;
		double cos2t = 0.0;
		if (hitpos_material.reflection_type == REFLECTION_TYPE_SSS_REFRACTION)
		{
			direct_light_meshod = 0;
			dir = normalize(ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir));
			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか
			// Snellの法則
			const double nc = 1.0;								// 真空の屈折率
			double nt = hitpos_material.refractive_index;		// オブジェクトの屈折率(kIor)

			nt = hitpos_material.refractive(wavelength);

			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = dot(ray.dir, orienting_normal);
			double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

			if ( cos2t >= 0.0 )
			{
				// 屈折の方向
				tdir = 	normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

				// SchlickによるFresnelの反射係数の近似を使う
				const double a = nt - nc, b = nt + nc;
				const double R0 = (a * a) / (b * b);

				const double c = 1.0 - (into ? -ddn : dot(tdir, -1.0 * orienting_normal));
				Re = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
				const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
				Tr = (1.0 - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

				// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
				// ロシアンルーレットで決定する。
				const double ref_probability = 0.25 + 0.5 * Re;
				Re = Re / ref_probability;
				Tr = Tr / (1.0 - ref_probability);
			}
			//tdir = ray.dir;
		}


		if ( direct_light_meshod && nextEventEstimation && env->light_list.list.size())
		{
			// 直接光の評価
			//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
			double light_probability = 1.0;
			for ( int i = 0; i < SHADOW_RAY_SAMPLING; i++ )
			{
				const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
				direct_light = direct_light + direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;
					
			}
			direct_light = direct_light / SHADOW_RAY_SAMPLING;
			//if ( shadow_ray_hit_count ) fprintf(stderr, "shadow_ray_hit_count %d\n", shadow_ray_hit_count);

			//関与媒質による減衰を考慮
			if (participatingMedia)
			{
				direct_light = env->participatingMediaParam.transmittanceRatio( hitpoint.distance,wavelength)* direct_light;
			}
		}
		//direct_light はすでに反射率が乗算済みなので、weightを掛ける必要はない
		//なので後でweight*incoming_radianceをするがincoming_radianceには加算しない
		//つまり　direct_light + weight*incoming_radiance

		if ( !direct_light_meshod && nextEventEstimation )
		{
			// 直接光の評価
			Intersection lid;
			if ( intersect_scene(Ray(hitpoint.position, dir), &lid) )
			{
				const Entity *now_object2 = EntList->List[lid.object_id];
				if ( now_object2->light )
				{
					direct_light =  RGB2Spectrum(env->light_list.list[now_object2->light_id].light->material()->emission, wavelength);
						
					//関与媒質による減衰を考慮
					if ( participatingMedia )
					{
						direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance,wavelength) * direct_light;
					}
				}
			}
		}

#if 10
		if (hitpos_material.reflection_type == REFLECTION_TYPE_SSS_REFRACTION)
		{
			if (cos2t < 0.0) 
			{
				const double probability = 0.5;
				if (rnd->next01() < probability )
				{
					// 全反射
					incoming_radiance = radiance(env, Ray(hitpoint.position,dir), rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
					weight = RGB2Spectrum(hitpos_material.color, wavelength) 
						/ probability / russian_roulette_probability/(1.0 - scattering_probability);

					//関与媒質を考慮する場合は光の減衰を考慮する
					if ( participatingMedia )
					{
						incoming_radiance = transmittance_ratio * incoming_radiance;
					}


					// emission は反射率が乗算されていないので後でweightを掛ける必要があるためincoming_radianceに加算しておく
					return (MULTIPLY(weight, incoming_radiance + direct_light));
				}else
				{
					Spectrum through_media;
					//表面下散乱(今レイが飛んできた方向のまま->中に入る)
					through_media = radiance_through_media(env, Ray(hitpoint.position, tdir), rnd, 0, now_object, 0, wavelength, nextEventEstimation, participatingMedia)*Tr;
							
					weight = const_cast<Entity*>(now_object)->material()->color / (1.0 - probability) / russian_roulette_probability / (1.0 - scattering_probability);


					//関与媒質を考慮する場合は光の減衰を考慮する
					if (participatingMedia)
					{
						through_media = transmittance_ratio * through_media;
					}
					return ( through_media+MULTIPLY(weight,direct_light));
				}
			}
		}
#endif

#if 10
		if (depth < 2) 
		{
			Spectrum through_media;
			//反射成分(普通に反射方向へ)
			incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength,nextEventEstimation, participatingMedia)*Re;
			weight = const_cast<Entity*>(now_object)->material()->color 
				/ russian_roulette_probability / (1.0 - scattering_probability);

			//表面下散乱(今レイが飛んできた方向のまま->中に入る)
			through_media = radiance_through_media(env, Ray(hitpoint.position, tdir), rnd, depth + 1, now_object, 0, wavelength, nextEventEstimation, participatingMedia)*Tr
				/ russian_roulette_probability
				/ (1.0 - scattering_probability);


			//関与媒質を考慮する場合は光の減衰を考慮する
			if ( participatingMedia )
			{
				incoming_radiance = transmittance_ratio * incoming_radiance;
				through_media = transmittance_ratio * through_media;
			}

			if ( direct_light_meshod )
			{
				direct_light = direct_light / russian_roulette_probability/(1.0 - scattering_probability);
				return (direct_light + MULTIPLY(weight, incoming_radiance) + through_media);
			}
			return (MULTIPLY(weight, incoming_radiance + direct_light) + through_media);
			
			//return 
			//	direct_light + (MULTIPLY( const_cast<Entity*>(now_object)->material()->color, radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength,nextEventEstimation, participatingMedia)) * 0.5 +
			//	radiance_through_media(env, Ray(hitpoint.position, ray.dir), rnd, 0, now_object, wavelength, nextEventEstimation, participatingMedia) * 0.5) / russian_roulette_probability/(1.0 - scattering_probability);
			//return 
			//	(MULTIPLY( const_cast<Entity*>(now_object)->material()->color, radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength,0,0)) * 0.5 +
			//	radiance_through_media(env, Ray(hitpoint.position, ray.dir), rnd, 0, now_object, wavelength, 0, 0) * 0.5) / russian_roulette_probability;
		} else
#endif
		{
			//fprintf(stderr, "%f %f\n", hitpos_material.participatingMediaParam.scattering_albedo,1.0-std::max(0.001, std::min(0.9, hitpos_material.participatingMediaParam.scattering_albedo)));

			//反射か表面下か確率的（５０％）
			const double probability = 0.5;
			if (rnd->next01() < probability )
			{
				//fprintf(stderr, "Re %f\n", Re);
				//反射(普通に反射)
				incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)*Re;
				weight = const_cast<Entity*>(now_object)->material()->color 
					/ probability / russian_roulette_probability / (1.0 - scattering_probability);
	

				//関与媒質を考慮する場合は光の減衰を考慮する
				if ( participatingMedia )
				{
					incoming_radiance = transmittance_ratio * incoming_radiance;
				}
				if ( direct_light_meshod )
				{
					direct_light = direct_light / probability / russian_roulette_probability / (1.0 - scattering_probability);
					return (direct_light + MULTIPLY(weight, incoming_radiance))*0.5;
				}
				return ( MULTIPLY(weight, incoming_radiance+direct_light))*0.5;

				//direct_light = direct_light / probability;
				//direct_light = direct_light * 0.5;
				//return direct_light + MULTIPLY(const_cast<Entity*>(now_object)->material()->color, radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength,0,0)) * 0.5
				//		/ probability
				//		/ russian_roulette_probability;
				//return MULTIPLY(const_cast<Entity*>(now_object)->material()->color, radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength,0,0)) * 0.5
				//		/ probability
				//		/ russian_roulette_probability;
			} else 
			{
				//fprintf(stderr, "Tr %f\n", Tr);
				Spectrum through_media;
				//表面下散乱(今レイが飛んできた方向のまま->中に入る)
				through_media = radiance_through_media(env, Ray(hitpoint.position, tdir), rnd, depth+1, now_object, 0, wavelength, nextEventEstimation, participatingMedia)*Tr
						/ (1.0 - probability)
						/ russian_roulette_probability / (1.0 - scattering_probability);
				weight = const_cast<Entity*>(now_object)->material()->color / (1.0 - probability) / russian_roulette_probability / (1.0 - scattering_probability);


				//関与媒質を考慮する場合は光の減衰を考慮する
				if (participatingMedia)
				{
					through_media = transmittance_ratio * through_media;
				}

				if ( direct_light_meshod )
				{
					//direct_light = direct_light / (1.0 - probability) / russian_roulette_probability / (1.0 - scattering_probability);
					//direct_light = ZERO();
					return ( through_media + MULTIPLY(weight, direct_light))*0.5;
					//return (direct_light + through_media)*0.5;
				}
				
				//direct_light = direct_light / (1.0 - probability) / russian_roulette_probability / (1.0 - scattering_probability);
				return ( through_media+MULTIPLY(weight,direct_light))*0.5;

				//direct_light = direct_light / (1.0 - probability);
				//
				//direct_light = direct_light * 0.5;
				//return radiance_through_media(env, Ray(hitpoint.position, ray.dir), rnd, 0, now_object, wavelength, nextEventEstimation, participatingMedia) * 0.5
				//		/ (1.0 - probability)
				//		/ russian_roulette_probability;
				//return radiance_through_media(env, Ray(hitpoint.position, ray.dir), rnd, 0, now_object, wavelength, 0, 0) * 0.5
				//		/ (1.0 - probability)
				//		/ russian_roulette_probability;
			}
		}
	}
	break;
	// Ward BRDF
	case REFLECTION_TYPE_WARD_BRD: {

		//fprintf(stderr, "%f %f %f,%f,%f\n", hitpos_material.ward_brdf.alp_x, hitpos_material.ward_brdf.alp_y, hitpos_material.ward_brdf.specular.x, hitpos_material.ward_brdf.specular.y, hitpos_material.ward_brdf.specular.z);

		const double alpha_x = hitpos_material.ward_brdf.alp_x;
		const double alpha_y = hitpos_material.ward_brdf.alp_y;

		const Vector3d in = -1.0 * ray.dir;
		Vector3d halfv;
		Vector3d dir;
		
		// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
		Vector3d w, u, v;
		OrthonormalBasis(orienting_normal, w, u, v);


#if 10	//重点サンプリング（Notes on the Ward BRDF)
		//Bruce Walter
		//Technical Report PCG-05-06
		//Cornell Program of Computer Graphics
		//April 29, 2005）

		// 重点サンプリングする。ただし、生成されるハーフベクトルしだいでは
		// 反射方向が半球外に出てしまうのでそういう場合は棄却する。
		do {
			const double u1 = rnd->next01();
			const double u2 = rnd->next01();
			double phi = atan(alpha_y / alpha_x * tan(PS_TWOPI * u2));
			if (0.25 <= u2 && u2 <= 0.75)
				phi += PS_PI;
			else if (0.75 < u2)
				phi += PS_TWOPI;
			const double theta = atan(sqrt(-log(u1) / (pow(cos(phi), 2) / pow(alpha_x, 2) + pow(sin(phi), 2) / pow(alpha_y, 2))));

			halfv = normalize((
				u * cos(phi) * sin(theta) + 
				v * sin(phi) * sin(theta) + 
				w * cos(theta)));

			dir = normalize(2.0 * dot(in, halfv) * halfv - in);
		} while (dot(orienting_normal, dir) < 0.0);

		// 重み。brdf * cosθ / pdf(ω)をまとめたものになっている。
		const double ww = dot(halfv, in) * pow(dot(halfv, orienting_normal), 3) * 
			sqrt(dot(dir, orienting_normal) / dot(in, orienting_normal));
		
		incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
		weight = RGB2Spectrum(hitpos_material.ward_brdf.specular, wavelength) * ww/russian_roulette_probability/(1.0 - scattering_probability);

		if ( nextEventEstimation && env->light_list.list.size() )
		{
			double light_probability = 1.0;
			//直接光サンプリングする
			//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
			for ( int i = 0; i < SHADOW_RAY_SAMPLING; i++ )
			{	
				const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
				direct_light = direct_light + direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) /  light_probability;
			}
			direct_light = direct_light / SHADOW_RAY_SAMPLING;
		}

		// direct_light はすでに反射率が乗算済みなので、weightを掛ける必要はない
		//なので後でweight*incoming_radianceをするがincoming_radianceには加算しない
		//つまり　direct_light + weight*incoming_radiance
		

#else	//一様サンプリング（拡散項付加)
		while( 1 )
		{
			// 生成されるハーフベクトル反射方向が半球外に出てしまう場合はやり直し。
			do {
				//一様サンプリング
				const double t = rnd->next01();
				const double s = PS_TWOPI * rnd->next01();
				const double tt = 1 - t;
				const double tt2 = sqrt(1.0 - tt*tt);

				halfv = normalize((
						u * cos(s) * tt2 + 
						v * sin(s) * tt2 + 
						w * tt));

				dir = normalize( 2.0 * dot(in, halfv) * halfv - in);
			} while (dot(orienting_normal, dir) < 0.0);

			double ww0 = pow(dot(halfv, orienting_normal),2.0);
			//if ( fabs(ww0) < 1.0e-14 )
			//{
			//	continue;
			//}
			double ww1 = pow(dot(halfv, u)/alpha_x,2.0) + pow(dot(halfv, v)/alpha_y,2.0)/ww0;
			double ww2 = 4.0*PS_PI*alpha_x*alpha_y*sqrt(dot(dir, orienting_normal)* dot(in, orienting_normal));
		
			//if ( fabs(ww2) < 1.0e-14 )
			//{
			//	continue;
			//}
			const double ww = exp(-ww1)/ww2;
			incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
			weight = RGB2Spectrum(hitpos_material.ward_brdf.specular, wavelength) * ww * dot(w, dir)*PS_TWOPI/russian_roulette_probability/(1.0 - scattering_probability);

			if ( nextEventEstimation && env->light_list.list.size() )
			{
				double light_probability = 1.0;
				//直接光サンプリングする
				for ( int i = 0; i < SHADOW_RAY_SAMPLING; i++ )
				{	
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
					direct_light = direct_light + direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) /  light_probability;
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
			}
			// direct_light はすでに反射率が乗算済みなので、weightを掛ける必要はない
			//なので後でweight*incoming_radianceをするがincoming_radianceには加算しない
			//つまり　direct_light + weight*incoming_radiance
	
			Spectrum ward =  MULTIPLY(weight,incoming_radiance);


			// レンダリング方程式に対するモンテカルロ積分を考えると、outgoing_radiance = weight * incoming_radiance。
			// ここで、weight = (ρ/π) * cosθ / pdf(ω) / R になる。 pdf(ω)=1/(2π)
			// ρ/πは完全拡散面のBRDFでρは反射率、cosθはレンダリング方程式におけるコサイン項
			weight = RGB2Spectrum(hitpos_material.color/PS_PI , wavelength)*dot(w, dir)*PS_TWOPI/russian_roulette_probability;

			Spectrum diffuse_term = MULTIPLY(weight,incoming_radiance);
			
			incoming_radiance = (diffuse_term + ward);
		}		
#endif
	}break;

	// 完全拡散面
	case REFLECTION_TYPE_DIFFUSE: {

		//平面光源の裏側は光らないようにする処理
		if ( now_object->light && now_object->type == ENTITY_TYPE_UVPLANE )
		{
			if ( dot(hitpoint.normal, ray.dir) >= 0.0 )
			{
				return ZERO();
			}
		}

		// 直接光を計算しているのでeye から直接 light に hit した場合のみ、emission を返する
		if ( nextEventEstimation )
		{
			if ( now_object->light )
			{
				if ( depth == 0 )
				{
					//return RGB2Spectrum( hitpos_material.emission, wavelength);
					return RGB2Spectrum(transmittance_ratio * hitpos_material.emission, wavelength);
				}else
				{
					return ZERO();
				}
			}
		}



		// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
		Vector3d w, u, v;
		OrthonormalBasis(orienting_normal, w, u, v);

#if 10
		// コサイン項を使った重点的サンプリング
		const double r1 = PS_TWOPI * rnd->next01();
		const double r2 = rnd->next01(), r2s = sqrt(r2);
		Vector3d dir = normalize((
			u * cos(r1) * r2s +
			v * sin(r1) * r2s +
			w * sqrt(1.0 - r2)));

		incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
		// レンダリング方程式に対するモンテカルロ積分を考えると、outgoing_radiance = weight * incoming_radiance。
		// ここで、weight = (ρ/π) * cosθ / pdf(ω) / R になる。
		// ρ/πは完全拡散面のBRDFでρは反射率、cosθはレンダリング方程式におけるコサイン項、pdf(ω)はサンプリング方向についての確率密度関数。
		// Rはロシアンルーレットの確率。
		// 今、コサイン項に比例した確率密度関数によるサンプリングを行っているため、pdf(ω) = cosθ/π
		// よって、weight = ρ/ R。
		weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability/(1.0 - scattering_probability);
#else
		//一様サンプリング
		Vector3d dir;
		//do {
			const double t = rnd->next01();
			const double s = PS_TWOPI * rnd->next01();
			const double tt = 1 - t;
			const double tt2 = sqrt(1.0 - tt*tt);
			dir = normalize((
					u * cos(s) * tt2 + 
					v * sin(s) * tt2 + 
					w * tt));
		//} while (dot(w, dir) < 0.0);

		incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
		// レンダリング方程式に対するモンテカルロ積分を考えると、outgoing_radiance = weight * incoming_radiance。
		// ここで、weight = (ρ/π) * cosθ / pdf(ω) / R になる。 pdf(ω)=1/(2π)
		// ρ/πは完全拡散面のBRDFでρは反射率、cosθはレンダリング方程式におけるコサイン項
		weight = RGB2Spectrum(hitpos_material.color/PS_PI , wavelength)*dot(w, dir)*PS_TWOPI/russian_roulette_probability;
#endif

		if ( nextEventEstimation && env->light_list.list.size() )
		{
			// 直接光の評価
			//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
			double light_probability = 1.0;
			for ( int i = 0; i < SHADOW_RAY_SAMPLING; i++ )
			{
				const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
				Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;
					
				direct_light = direct_light + dl;
			}
			direct_light = direct_light / SHADOW_RAY_SAMPLING;
			direct_light = direct_light / russian_roulette_probability/(1.0 - scattering_probability);
			//if ( shadow_ray_hit_count ) fprintf(stderr, "shadow_ray_hit_count %d\n", shadow_ray_hit_count);
		}
		// direct_light はすでに反射率が乗算済みなので、weightを掛ける必要はない
		//なので後でweight*incoming_radianceをするがincoming_radianceには加算しない
		//つまり　direct_light + weight*incoming_radiance


#if 0
		//光源が比較的隠れてているような状況で光源の面積が非常に大きいと光源の位置サンプリングで
		//光源からのレイがヒットしないようなケースが殆どになってしまう
		//その場合は今の衝突点からのレイを飛ばして光源にヒットする可能性高い場合は以下のように考慮する。
		//※このような現象はありえそうだがnextEventEstimationの計算が間違っていてそうなるのかも知れない

		//↑はやはり間違いだった(2014.9.10) 平面光源の座標変換で元平面は変換していなかったためヒットしなかった。
		//nextEventEstimationを指定しないと光源にヒットしてからもトレースをするためどこかの経路でヒットする可能性が大きかったが
		//nextEventEstimationを指定していると光源からのトレースは無いためヒットする可能性は殆ど無かった
		if ( nextEventEstimation && shadow_ray_hit_count == 0 )
		{
			const bool flg = true;

			// 直接光の評価
			Intersection lid;
			if ( intersect_scene(Ray(hitpoint.position, dir), &lid) )
			{
				const Entity *now_object2 = EntList->List[lid.object_id];
				if ( now_object2->light )
				{
					if ( flg )
					{
						Color emit = env->light_list.list[now_object2->light_id].light->material.emission;
						//反射率が乗算されていないので後で反射率を掛ける
						//incoming_radiance =  incoming_radiance + emit / (SHADOW_RAY_SAMPLING+1);

						direct_light = emit / (SHADOW_RAY_SAMPLING+1) / russian_roulette_probability/(1.0 - scattering_probability);
					}else
					{
						const Vector3d light_pos = lid.hitpoint.position;
						const Light& light = env->light_list.list[now_object2->light_id];

						direct_light = eval_direct_radiance(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, light_pos, depth, wavelength);
						direct_light = direct_light / (SHADOW_RAY_SAMPLING+1) / russian_roulette_probability/(1.0 - scattering_probability);
					}
				}
			}
		}
#endif


	} break;

	// 完全鏡面
	case REFLECTION_TYPE_SPECULAR:
	case REFLECTION_TYPE_SPECULAR_DIFFUSE: {

		double r_probability = 0.2;
		if (hitpos_material.reflection_type == REFLECTION_TYPE_SPECULAR) r_probability = 1.0;

		if (rnd->next01() < r_probability)
		{
			const Vector3d dir = ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir);
			// 完全鏡面なのでレイの反射方向は決定的。
			// ロシアンルーレットの確率で除算するのは上と同じ。
			incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);

			weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability / (1.0 - scattering_probability) / r_probability;


			if (nextEventEstimation)
			{
				//const bool flg = (!participatingMedia)?true:false;
				const bool flg = true;


				// 直接光の評価
				Intersection lid;
				if (intersect_scene(Ray(hitpoint.position, dir), &lid))
				{
					const Entity *now_object2 = EntList->List[lid.object_id];
					if (now_object2->light)
					{
						if (flg)
						{
							direct_light = RGB2Spectrum(lid.hitpoint.material.emission, wavelength);

							//関与媒質による減衰を考慮
							if (participatingMedia)
							{
								direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance, wavelength) * direct_light;
							}
						}
						else
						{
							const Vector3d light_pos = lid.hitpoint.position;
							const Light& light = env->light_list.list[now_object2->light_id];

							direct_light = eval_direct_radiance(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, light_pos, depth, wavelength);
							direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
						}
					}
				}
				if (flg)
				{
					// emission は反射率が乗算されていないので後でweightを掛ける必要があるためincoming_radianceに加算しておく
					incoming_radiance = incoming_radiance + direct_light;
					direct_light = ZERO();
				}
			}
			else
			{
				//平面光源の裏側は光らないようにする処理
				if (now_object->light && now_object->type == ENTITY_TYPE_UVPLANE)
				{
					if (dot(hitpoint.normal, ray.dir) >= 0.0)
					{
						return ZERO();
					}
				}

				// 直接光を計算しているのでeye から直接 light に hit した場合のみ、emission を返する
				if (nextEventEstimation)
				{
					if (now_object->light)
					{
						if (depth == 0)
						{
							//return RGB2Spectrum( hitpos_material.emission, wavelength);
							return RGB2Spectrum(transmittance_ratio * hitpos_material.emission, wavelength) / (1.0 - r_probability);
						}
						else
						{
							return ZERO();
						}
					}
				}



				// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
				Vector3d w, u, v;
				OrthonormalBasis(orienting_normal, w, u, v);

				// コサイン項を使った重点的サンプリング
				const double r1 = PS_TWOPI * rnd->next01();
				const double r2 = rnd->next01(), r2s = sqrt(r2);
				Vector3d dir = normalize((
					u * cos(r1) * r2s +
					v * sin(r1) * r2s +
					w * sqrt(1.0 - r2)));

				incoming_radiance = radiance(env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
				// レンダリング方程式に対するモンテカルロ積分を考えると、outgoing_radiance = weight * incoming_radiance。
				// ここで、weight = (ρ/π) * cosθ / pdf(ω) / R になる。
				// ρ/πは完全拡散面のBRDFでρは反射率、cosθはレンダリング方程式におけるコサイン項、pdf(ω)はサンプリング方向についての確率密度関数。
				// Rはロシアンルーレットの確率。
				// 今、コサイン項に比例した確率密度関数によるサンプリングを行っているため、pdf(ω) = cosθ/π
				// よって、weight = ρ/ R。
				weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability / (1.0 - scattering_probability) / (1.0 - r_probability);

				if (nextEventEstimation && env->light_list.list.size())
				{
					// 直接光の評価
					//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
					double light_probability = 1.0;
					for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
					{
						const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
						Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

						direct_light = direct_light + dl;
					}
					direct_light = direct_light / SHADOW_RAY_SAMPLING;
					direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability) / (1.0 - r_probability);
					//if ( shadow_ray_hit_count ) fprintf(stderr, "shadow_ray_hit_count %d\n", shadow_ray_hit_count);
				}
				// direct_light はすでに反射率が乗算済みなので、weightを掛ける必要はない
				//なので後でweight*incoming_radianceをするがincoming_radianceには加算しない
				//つまり　direct_light + weight*incoming_radiance

			}
		}

	} break;

	// 屈折
	case REFLECTION_TYPE_REFRACTION: {

		const Ray reflection_ray = Ray(hitpoint.position, normalize(ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir)));
		const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

		//const bool flg = (!participatingMedia)?true:false;
		const bool flg = true;

		// 反射方向の直接光の評価
		if ( nextEventEstimation )
		{
			Intersection lid;
			if ( intersect_scene(reflection_ray, &lid) )
			{
				const Entity *now_object2 = EntList->List[lid.object_id];
				if ( now_object2->light )
				{
					if ( flg )
					{
						direct_light =  RGB2Spectrum(lid.hitpoint.material.emission, wavelength);
						//関与媒質による減衰を考慮
						if (participatingMedia )
						{
							direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance,wavelength) * direct_light;
						}
					}else
					{
						const Vector3d light_pos = lid.hitpoint.position;
						const Light& light = env->light_list.list[now_object2->light_id];
						direct_light = eval_direct_radiance(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, reflection_ray.dir, light_pos, depth, wavelength);
					}
				}
			}
		}

		// Snellの法則
		const double nc = 1.0;								// 真空の屈折率
		double nt = hitpos_material.refractive_index;		// オブジェクトの屈折率(kIor)
		
		nt = hitpos_material.refractive(wavelength);

		const double nnt = into ? nc / nt : nt / nc;
		const double ddn = dot(ray.dir, orienting_normal);
		const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
		
		//printf("now_object->refractive_index %f\n", now_object->refractive_index);
		if (cos2t < 0.0) 
		{
			// 全反射
			incoming_radiance = radiance(env, reflection_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia);
			weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability/(1.0 - scattering_probability);

			if ( flg )
			{
				// emission は反射率が乗算されていないので後でweightを掛ける必要があるためincoming_radianceに加算しておく
				incoming_radiance = incoming_radiance + direct_light;
				direct_light = ZERO();
			}else
			{
				direct_light = direct_light / russian_roulette_probability/(1.0 - scattering_probability);
			}
			break;
		}

		// 屈折の方向
		const Ray tdir = Ray(hitpoint.position,
			normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(fabs(cos2t)))));
		
		// SchlickによるFresnelの反射係数の近似を使う
		const double a = nt - nc, b = nt + nc;
		const double R0 = (a * a) / (b * b);

		const double c = 1.0 - (into ? -ddn : dot(tdir.dir, -1.0 * orienting_normal));
		const double Re = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
		const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
		const double Tr = (1.0 - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合
		
		// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
		// ロシアンルーレットで決定する。
		const double probability  = 0.25 + 0.5 * Re;

		// 屈折方向の直接光の評価
		if ( nextEventEstimation )
		{
			Intersection lid;
			if ( intersect_scene(tdir, &lid) )
			{
				const Entity *now_object3 = EntList->List[lid.object_id];
				if ( now_object3->light )
				{
					if ( flg )
					{
						direct_light_refraction =  RGB2Spectrum(lid.hitpoint.material.emission, wavelength);
						
						//関与媒質による減衰を考慮
						if (participatingMedia)
						{
							direct_light_refraction = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance,wavelength) * direct_light_refraction;
						}
					}else
					{
						const Vector3d light_pos = lid.hitpoint.position;
						const Light& light = env->light_list.list[now_object3->light_id];
						direct_light_refraction = eval_direct_radiance(env, ray.dir, rnd, light, now_object, hitpoint, -1.0 * orienting_normal, tdir.dir, light_pos, depth, wavelength);
					}
				}
			}
		}

		if (depth > 2) {
			if (rnd->next01() < probability) { // 反射
				incoming_radiance = (radiance(env, reflection_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)) * Re;
				weight = RGB2Spectrum(hitpos_material.color, wavelength)/ (probability * russian_roulette_probability)/(1.0 - scattering_probability);

				if ( flg )
				{
					// emission は反射率が乗算されていないので後でweightを掛ける必要があるためincoming_radianceに加算しておく
					incoming_radiance = incoming_radiance + direct_light;
					direct_light = ZERO();
					direct_light_refraction = ZERO();
				}else
				{
					direct_light = direct_light / russian_roulette_probability/(1.0 - scattering_probability);
				}
			} else { // 屈折
				incoming_radiance = (radiance(env, tdir, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)) * Tr;
				weight = RGB2Spectrum(hitpos_material.color, wavelength) / ((1.0 - probability) * russian_roulette_probability)/(1.0 - scattering_probability);
				
				if ( flg )
				{
					// emission は反射率が乗算されていないので後でweightを掛ける必要があるためincoming_radianceに加算しておく
					incoming_radiance = incoming_radiance + direct_light_refraction;
					direct_light = ZERO();
					direct_light_refraction = ZERO();
				}else
				{
					direct_light_refraction = direct_light_refraction / russian_roulette_probability/(1.0 - scattering_probability);
				}
			}
		} else { // 屈折と反射の両方を追跡
			incoming_radiance = 
				(radiance(env, reflection_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)) * Re +
				(radiance(env, tdir, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)) * Tr;
			weight = RGB2Spectrum(hitpos_material.color, wavelength) / (russian_roulette_probability)/(1.0 - scattering_probability);
		
			if ( flg )
			{
				// emission は反射率が乗算されていないので後でweightを掛ける必要があるためincoming_radianceに加算しておく
				incoming_radiance = incoming_radiance + direct_light + direct_light_refraction;
				direct_light = ZERO();
				direct_light_refraction = ZERO();
			}else
			{
				direct_light = direct_light / russian_roulette_probability/(1.0 - scattering_probability);
				direct_light_refraction = direct_light_refraction / russian_roulette_probability/(1.0 - scattering_probability);
			}
		}
	} break;

	}

	//関与媒質を考慮する場合は光の減衰を考慮する
	if ( participatingMedia )
	{
		//環境内（物体の外側)
		//現在の交点ではレイは物体の中に入ろうとする方向なら物体の外側から
		//if ( /*inside_dir &&*/ hitpos_material.reflection_type != REFLECTION_TYPE_REFRACTION && hitpos_material.reflection_type != REFLECTION_TYPE_SPECULAR )
		{
			emission = transmittance_ratio * emission;
			incoming_radiance = transmittance_ratio * incoming_radiance;
		}
	}
	return RGB2Spectrum(emission,wavelength)+direct_light+direct_light_refraction + MULTIPLY(weight,incoming_radiance);
}

};

#endif
