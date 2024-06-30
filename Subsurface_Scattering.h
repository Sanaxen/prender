#ifndef _SUBSURFACE_SCATTERING_H
#define _SUBSURFACE_SCATTERING_H

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

#include "bssrdf.h"

#define OVER_RADIANCE_CUT0	100
#define OVER_RADIANCE_CUT1	200
#define OVER_RADIANCE_CUT2	100
#define OVER_RADIANCE_CUT3	200

namespace prender 
{
Spectrum radiance(int* direct_light_hit, const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);
inline Spectrum BRDF_Refraction(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);
inline Spectrum BRDF_Diffuse(int* hit_direct_light, const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);
inline double fresnelDiffuseReflectance(double eta);

inline Spectrum radiance_through_media(const SceneEnv* env, const IntersectionPos &hitpoint, const Ray &ray, Random *rnd, const int depth, const Entity* in_object, const int sc_depth,  const int sc_max, const double russian_roulette_probability, const double scattering_probability, const double wavelength, double& totalDist, const int nextEventEstimation = 0, const int participatingMedia = 0)
{
	const double eps = 0.00001;
	const ParticipatingMedia& participating_media = const_cast<Entity*>(in_object)->material()->participatingMediaParam;
	
	double dirPdf = 1.0;
	double phase = 1.0;


	double russian_roulette_probability2 = 1.0;
	if (depth > env->DepthLimit * 30)
	{
		russian_roulette_probability2 *= std::min(pow(0.5, depth - env->DepthLimit), pow(0.5, sc_depth));

		//fprintf(stderr, "depth %d scattering %d\n", depth, sc_depth);
		if (rnd->next01() >= russian_roulette_probability2)
		{
			//fprintf(stderr, "DepthLimit(%dx30)!! depth %d scattering %d\n", env->DepthLimit, depth, sc_depth);
			//fprintf(stderr, "\n");
			return RGB2Spectrum(const_cast<Entity*>(in_object)->material()->emission, wavelength);
		}
	}
	else
		russian_roulette_probability2 = 1.0; // ロシアンルーレット実行しなかった


	//if ( sc_depth > 10 )
	//{
	//	fprintf(stderr, "depth %d scattering %d\n", depth, sc_depth);
	//}

	//かならずオブジェクトの中である必要があるためチェック。
	bool inpos = false;
	Intersection out;
	Ray tmp_ray = ray;

	if (!intersect_scene(tmp_ray, &out, depth))
	{
		return ZERO();
	}

	Spectrum Li, L;

	//const int select = participating_media.sample(rnd->next01(),0);

	double σt = participating_media.transmittanceColor[participating_media.sample(rnd->next01(),0)];

	double probability = 1.0 - exp(-σt*out.hitpoint.distance);
	if (probability < 0) probability = 0;
	else if (probability > 1.0) probability = 1.0;

	//double absorbing_p = 1.0 - exp(-participating_media.absorbingColor.Av() * out.hitpoint.distance);
	//if (absorbing_p > probability && !(sc_depth == 0 || sc_depth == sc_max))
	//{
	//	return ZERO();
	//}

	//const double σt_av = participating_media.transmittanceColor.Av();
	//if (/*σt_av > 2.0 &&*/ probability > 0.99)
	//{
	//	probability = 0.99;
	//}

	//participating_media.scattering / participating_media.transmittance;

	//probability = 0;
	//fprintf(stderr, "%f\n", participating_media.scattering_albedo);


	if (sc_depth == 0) probability = 1.0;
	if (sc_depth == sc_max || out.hitpoint.distance < 1.0e-4) probability = 0;


	const double T21 = std::max(1.0 - fresnel(dot(ray.dir*-1.0, hitpoint.normal), hitpoint.material.refractive_index), 0.0);

	int stat = -1;
	//表面位置から中に入って行く(ray方向)
	bool sss = (rnd->next01() < probability);
	if ( sss )
	{
		//散乱位置への移動
		double d = 0.0;
		double pdf;
		double p_pdf=1;

		//あたらなレイを飛ばす位置（散乱位置）
		Vector3d new_org =ray.org;
		//散乱方向
		Ray next_ray = ray;

		sss = false;
		int trynum = 145;// 50;
		//while(1)
		for (int kk = 0; kk < trynum; kk++)
		{
			//while(1)//現在の位置と散乱位置が表面を跨がないように選択するのでこのループは本当は不要
			//{
			//	//現在の位置と散乱位置が表面を跨がないように選択
			//	const double u = rnd->next01();
			//	d = -log(1.0 + u * (exp(-participating_media.transmittance * out.hitpoint.distance) - 1.0)) / participating_media.transmittance;
			//	pdf = exp(-participating_media.transmittance * d) * (-participating_media.transmittance / (exp(-participating_media.transmittance * out.hitpoint.distance) - 1.0));
			//
			//	if ( d < out.hitpoint.distance) break;
			//}

			//const double u = rnd->next01();
			//d = -log(1.0 - u*(1.0 - exp(-out.hitpoint.distance*participating_media.transmittanceColor[select]))) / participating_media.transmittanceColor[select];
			//pdf = participating_media.transmittanceColor[select] / (exp(d*participating_media.transmittanceColor[select]) - exp((d - out.hitpoint.distance)*participating_media.transmittanceColor[select]));


			double pdf1, d1;
			double pdf2, d2;
			bool ss = false;
			
			//double σt = participating_media.transmittanceColor.Min();
			int ii = participating_media.sample(rnd->next01(), SAMPLE_ABSORBING, 1);
			double σt = participating_media.transmittanceColor[ii];
			//double σt = participating_media.transmittanceColor[std::min(2, (int)(3.0*rnd->next01()))];


			if (0)
			{
				for (int j = 0; j < 50; j++)
				{
					double t = -log(rnd->next01()) / σt;
					if (t < out.hitpoint.distance)
					{
						d1 = t;
						pdf1 = σt*exp(-σt*t);
						//fprintf(stderr, "%d\n", j);
						ss = true;
						break;
					}
				}
				//if (!ss ) fprintf(stderr, "--------------\n");
			}
			if ( 1 )
			{
				//double σt = participating_media.transmittanceColor.Max();
				int ii = participating_media.sample(rnd->next01(),SAMPLE_ABSORBING,1);
				double σt = participating_media.transmittanceColor[ii];
				//double σt = participating_media.transmittanceColor[std::min(2, (int)(3.0*rnd->next01()))];

				bool ss = false;
				//for (int j = 0; j < 1000; j++)
				//{
				//	double t = -log(rnd->next01()) / σt;
				//	if (t < out.hitpoint.distance)
				//	{
				//		d2 = t;
				//		pdf2 = σt*exp(-σt*t);
				//		//fprintf(stderr, "%d\n", j);
				//		ss = true;
				//		break;
				//	}
				//}

				if (!ss )
				{
					double b = out.hitpoint.distance;

					double t = -log(1.0 + rnd->next01() * (exp(-σt * b) - 1.0)) / σt;
					pdf2 = exp(-σt * t) * (-σt / (exp(-σt * b) - 1.0));
					d2 = t;
				}
			}
			if ( 0 )
			{
				//return ZERO();
				//double σt = participating_media.transmittanceColor.Max();
				int ii = participating_media.sample(rnd->next01(), SAMPLE_TRANSMITTANCE, 1);
				double σt = participating_media.transmittanceColor[ii];
				//double σt = participating_media.transmittanceColor[std::min(2, (int)(3.0*rnd->next01()))];
				
				double a = 0.0;
				double b = out.hitpoint.distance;

				double t = a - (1.0 / σt)*log(1.0 - rnd->next01()*(1.0 - exp(-(b - a)*σt)));
				pdf2 = σt / (exp((t - a)*σt) - exp((t - b)*σt));
				d2 = t;
			}
			//pdf = pdf2;
			//d = d2;
			//if (!ss)
			//{
			//	fprintf(stderr, "-------------\n");
			//	return ZERO();
			//}
			//pdf = pdf1;
			//d = d1;

#if 10
			if (ss)
			{
				double p = 0.0;
				//if ( participating_media.absorbingColor.Max() < 0.14 ) p = 0.7;

				//const double k = 200000.0;
				//const double n = 3.65;
				//const double x = 20000.0;
				//double c = 0.95;
				//double absorbing = participating_media.absorbingColor.Min();

				//double p = std::min(c,1.05-exp(-x*k*pow(absorbing,n)));

				if (rnd->next01() < p)
				{
					pdf = pdf1;
					d = d1;
				}else
				{
					pdf = pdf2;
					d = d2;
				}
			}
			else
			{
				pdf = pdf2;
				d = d2;
			}
#endif

			//あたらなレイを飛ばす位置（散乱位置）
			new_org = ray.org + d * ray.dir;
			
			//散乱方向
			next_ray = ray;

			//現在の位置と散乱位置が表面を跨いでいる？(つまり表面を飛び出した）
			//現在の位置と散乱位置が表面を跨がないように選択しているので殆ど不要だがSingleScatterring(BSSRDF)で
			//使う時に物体の中に物体があるようなケースでは判定しておく必要がある
			if (sc_depth < 1)
			{
				Intersection intersectionWrk;
				Ray temp_ray(ray.org, ray.dir);

				if (intersect_scene(temp_ray, &intersectionWrk, depth))
				{
					if (intersectionWrk.hitpoint.distance < (new_org - ray.org).length() && env->EntList.List[intersectionWrk.object_id] == in_object)
					{
						//sss = false;
						//probability = 0;
						//stat = 0;
						//break;
						//跨いでいるのはNG!!
						continue;
					}
				}
			}

			const int select = participating_media.sample(rnd->next01(),0);
			if (participating_media.phase_prm[select] == 0.0 || fabs(participating_media.phase_prm[select]) > 1.0)
			{
				// 等方散乱
				//const double r1 = 2 * PS_PI * rnd->next01();
				//const double r2 = 1.0 - 2.0 * rnd->next01();

				//
				//new_org = ray.org + d * ray.dir;
				//next_ray = Ray(new_org, normalize(Vector3d(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2)));

				next_ray.org = new_org;
				sampleSphere(rnd, ray.dir, next_ray.dir);
				phase = PhaseFunction_Isotropic();
				dirPdf = phase;
			}
			else
			{
				next_ray.org = new_org;
				sampleHG(rnd, ray.dir, participating_media.phase_prm[select], next_ray.dir);
				phase = PhaseFunction_HenyeyGreenstein(env, participating_media.phase_prm[select], dot(ray.dir, next_ray.dir));
				dirPdf = phase;
			}

			//現在の位置と散乱位置が表面を跨がないように選択しているので殆ど不要だがSingleScatterring(BSSRDF)で
			//使う時に物体の中に物体があるようなケースでは判定しておく必要がある
			if (/*sc_depth < */1)
			{
				Intersection intersectionWrk;
				if (!intersect_scene(next_ray, &intersectionWrk, depth))
				{
					continue;
					fprintf(stderr, "-Impossible!!\n");
					// 絶対ここには入らないはず
					return RGB2Spectrum(env->BackgroundColor, wavelength);
				}
				Entity* now_object = env->EntList.List[intersectionWrk.object_id];
				bool insidedir = now_object->isInsideDir(next_ray, intersectionWrk.hitpoint);

				insidedir = now_object->isInsideDir(next_ray, intersectionWrk.hitpoint);

				if (now_object == in_object)
				{
					if (insidedir)
					{
						//ここに入るのは散乱位置が悪すぎてオブジェクトの外側で方向はオブジェクトに向かう場合
						//散乱の結果、レイは表面にあたる（これはあり得ないはず)
						//fprintf(stderr, "-Impossible( ->) ->()\n");
						continue;
						return ZERO();

					}
					else
					{
						sss = true;
						stat = 0;
						break;
					}
				}
				else
				{
					//ここに入るのは散乱位置が悪すぎて他のオブジェクトに当たった
					//これはオブジェクトの中に他のオブジェクトがある場合に起きる場合がある。
					if (insidedir)
					{
						//他のオブジェクトの外側の場合はＯＫ
						sss = true;
						stat = 0;
						break;
					}
					else
					{
						//他のオブジェクトの内側から。これはあり得ない
						//中からレイを飛ばしているのに他の物体に当たった（ありえないはず、、、)
						//fprintf(stderr, "-Impossible( ->) ->{} or ->{->}%f (d=%f)\n", out.hitpoint.distance, d);
						continue;
						//return ZERO();
					}
				}
			}
			else
			{
				sss = true;
				stat = 0;
				break;
			}
		}
		if (!sss)
		{
			fprintf(stderr, "***Impossible\n");
			return ZERO();
		}
		if (stat == -1)
		{
			fprintf(stderr, "***Impossible------------\n");
			return ZERO();
		}
		if ( sss )
		{
			//fprintf(stderr, "/内部で減衰して散乱\n");
			//内部で減衰して散乱(in-scatting)
			const Spectrum Tr = participating_media.transmittanceRatio(d, wavelength);
			const Spectrum Ld = radiance_through_media(env, hitpoint, next_ray, rnd, depth + 1, in_object, sc_depth + 1, sc_max, russian_roulette_probability, scattering_probability, wavelength, totalDist, nextEventEstimation, participatingMedia);

			totalDist += d;
			Li = Tr * RGB2Spectrum(participating_media.scatteringColor, wavelength)
				* Ld
				* phase
				/ pdf
				/ dirPdf
				/ probability
				/ russian_roulette_probability2/p_pdf;

#if OVER_RADIANCE_CUT0
#ifndef SPECTRUM_USE
			if (Li.length() > OVER_RADIANCE_CUT0)
			{
				//fprintf(stderr, "pdf %f dirPdf %f probability %f russian_roulette_probability2 %f Ld %f\n",
				//	pdf, dirPdf, probability, russian_roulette_probability2, Ld.length());
				return ZERO();
			}
#else
			if ( Spectrum2RGB(Li).length() > OVER_RADIANCE_CUT0 ) return ZERO();
#endif
#endif
			return Li;
		}
	}

	if ( !sss )
	{
		totalDist += out.hitpoint.distance;
		//散乱しなかった場合

		//出口から飛ばす
		//Ray next_ray(out.hitpoint.position - 0*out.hitpoint.normal*0.0001, ray.dir);
		Spectrum L = ZERO();

		const Spectrum Tr = participating_media.transmittanceRatio(out.hitpoint.distance, wavelength);

		if (env->EntList.List[out.object_id] != in_object)
		{
			Ray next_ray(out.hitpoint.position + 0 * out.hitpoint.normal*0.0001, ray.dir);
			return Tr*radiance(0, env, next_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)
				/ (1.0 - probability)
				/ russian_roulette_probability2;
		}
		else
		{
			const Vector3d out_orienting_normal = dot(out.hitpoint.normal , ray.dir) < 0.0 ? out.hitpoint.normal: (-1.0 * out.hitpoint.normal); 
			//const Ray reflection_ray = Ray(out.hitpoint.position + 1 * out.hitpoint.normal*0.0001, normalize(ray.dir - out.hitpoint.normal*(-1.0) * 2.0 * dot(out.hitpoint.normal*(-1.0), ray.dir)));
			
			const Ray reflection_ray = Ray(out.hitpoint.position, normalize(ray.dir - out.hitpoint.normal * 2.0 * dot(out.hitpoint.normal, ray.dir)));

			Ray next_ray(out.hitpoint.position + 0*out.hitpoint.normal*0.0001, ray.dir);
			double Re2 = 1.0;
			double Tr2 = 1.0;
			{
				// Snellの法則
				const double nc = 1.0;								// 真空の屈折率
				double nt = out.hitpoint.material.refractive_index;		// オブジェクトの屈折率(kIor)

				nt = out.hitpoint.material.refractive(wavelength);
				bool into = dot(out.hitpoint.normal, ray.dir) < 0.0;

				const double nnt = into ? nc / nt : nt / nc;
				const double ddn = dot(ray.dir, out_orienting_normal);
				const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

				// SchlickによるFresnelの反射係数の近似を使う
				const double a = nt - nc, b = nt + nc;
				const double R0 = (a * a) / (b * b);

				// 屈折の方向
				if (cos2t > 0)
				{
					next_ray.dir = normalize(ray.dir * nnt - out.hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(fabs(cos2t))));
				}
				const double c = 1.0 - (into ? -ddn : dot(next_ray.dir, out_orienting_normal*-1.0));
				Re2 = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
				const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
				Tr2 = (1.0 - Re2) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

				if (cos2t < 0)
				{
					Tr2 = 0.0;
					Re2 = 1.0;

					Tr2 = Re2;
					next_ray = reflection_ray;
				}
			}
			//内部で減衰してそのまま抜ける
			//fprintf(stderr, "内部で減衰してそのまま抜ける\n");

			double cos_wi_n = dot(out.hitpoint.normal, next_ray.dir);

			if (cos_wi_n < 0)
			{
				//return ZERO();
			}
			//out.hitpoint.material.color = out.hitpoint.material.color*participating_media.albedoColor;
			const double T12 = std::max(1.0 - fresnel(cos_wi_n, hitpoint.material.refractive_index), 0.0);

			Spectrum weight = Tr*RGB2Spectrum(out.hitpoint.material.color, wavelength)
				/ (1.0 - probability) / russian_roulette_probability2;
			
			Li = radiance(0, env, next_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);

			Spectrum direct_light = ZERO();
			if( 0)
			{
				if (nextEventEstimation)
				{
					//out.hitpoint.material.color = participating_media.albedoColor;
					direct_light = direct_light_perfect_refrection(env, out.hitpoint, next_ray, depth, wavelength, participatingMedia);
				}
			}
			if (1)
			{
				const Vector3d orienting_normal = out.hitpoint.normal;// dot(out.hitpoint.normal, ray.dir) < 0.0 ? out.hitpoint.normal : (-1.0 * out.hitpoint.normal);

				if (nextEventEstimation && env->light_list.list.size())
				{
					// 直接光の評価
					//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
					double light_probability = 1.0;
					for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
					{
						const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);

						//weightが係って返ってくる
						Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, in_object, out.hitpoint, orienting_normal, next_ray.dir, depth, wavelength) / light_probability;

						direct_light = direct_light + dl;
					}
					direct_light = direct_light / SHADOW_RAY_SAMPLING;
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
							Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, in_object, out.hitpoint, orienting_normal, next_ray.dir, depth, wavelength) / light_probability;

							direct_light = direct_light + dl;
						}
					}
					direct_light = direct_light / SHADOW_RAY_SAMPLING;
				}
			}
			return (direct_light + Li)*weight;
		}
	}
}

inline Spectrum BRDF_SSS2(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* in_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia, const bool bump);


inline Spectrum BRDF_Refraction(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);

//蜂須賀先生のdirpole.cpp と比較
//#define DEMO_TEST
inline Spectrum BRDF_SSS2_(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* in_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia, const bool bump);
inline Spectrum BRDF_SSS2(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* in_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0, const bool bump=false)
{
	Spectrum c = ZERO();

	//const double eta = hitpoint.material.r_refractive(wavelength);
	//const double T21 = 1.0 - fresnel(dot(-1.0*ray.dir, orienting_normal), eta);
	const int num_samps = 1;
	for (int i = 0; i < num_samps; i++)
	{
		c = c + BRDF_SSS2_(env, rnd, ray, hitpoint, in_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia, bump);
	}
	c = c / (double)num_samps;
	return c;
}

inline Spectrum BRDF_Specular(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);

inline Vector3d uniformSampleDisk(double u1, double u2) {
	double r, theta;
	// transform [0, 1] sample to [-1, 1]
	double x = 2.0f * u1 - 1.0f;
	double y = 2.0f * u2 - 1.0f;

	if (x + y > 0) {
		if (x > y) {
			// right quarter of square
			r = x;
			theta = 0.25f * PS_PI * (y / x);
		}
		else {
			// up quarter of square
			r = y;
			// (-pi / 4) * (x / y) + pi / 2
			theta = 0.25f * PS_PI * (2.0f - x / y);
		}
	}
	else {
		if (x < y) {
			// left quarter of square
			r = -x;
			// (pi / 4) * (y / x) + pi
			theta = 0.25f * PS_PI * (4.0f + y / x);
		}
		else {
			// down quarter of square
			r = -y;
			// x and y may both be 0
			if (y != 0.0f) {
				// (-pi / 4) * (x / y) + (3pi / 2)
				theta = 0.25f * PS_PI * (6.0f - x / y);
			}
			else {
				theta = 0.0f;
			}
		}
	}

	return Vector3d(r * cos(theta), r * sin(theta), 0.0);
}
inline  Vector3d gaussianSample2D(double u1, double u2, double falloff, double Rmax) 
{
    double r = sqrt(log(1.0 - u1 * (1.0 - exp(-falloff * Rmax * Rmax))) /  -falloff);
	double theta = PS_TWOPI * u2;
    return Vector3d(r * cos(theta), r * sin(theta), 0);
}
inline double gaussianSample2DPdf(const Vector3d& pCenter, const Vector3d& pSample, const Vector3d& N, double falloff) 
{
        Vector3d d = pSample - pCenter;
        Vector3d projected = d - N * dot(d, N);
		return falloff * exp(-falloff * (projected).sqr())/PS_PI;
}
inline double gaussianSample2DPdf(const Vector3d& pCenter, const Vector3d& pSample, const Vector3d& N, double falloff, double Rmax) 
{
	return gaussianSample2DPdf(pCenter, pSample, N, falloff) / (1.0f - exp(-falloff * Rmax * Rmax));
}


inline Spectrum BRDF_SSS2_(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* in_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0, const bool bump = true)
{
	const Material& hitpos_material = hitpoint.material;
	Color emission = ZERO();					//自己発光
	Spectrum direct_light = ZERO();				//直接光(反射)
	Spectrum direct_light_refraction = ZERO();	//直接光(屈折)
	Spectrum incoming_radiance = ZERO();
	Spectrum weight = ONE();

	// no absorption inside glass
	bool inside_dir = const_cast<Entity*>(in_object)->isInsideDir2(ray, hitpoint);

	emission = hitpos_material.emission;		// Init color to emissive value

	Spectrum transmittance_ratio = ONE();
	//関与媒質を考慮する場合は光の減衰を考慮する
	if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

	// 相対屈折率 
	const double eta = hitpos_material.r_refractive(wavelength);

	DirectionalDipole directional_dipole(hitpos_material, eta);
	const ParticipatingMedia& participating_media = const_cast<Entity*>(in_object)->material()->participatingMediaParam;

	//const int sigma_sample = participating_media.sample_transmittance(rnd->next01(), 1);
	//const double σt = directional_dipole.A*directional_dipole.sigma_t[sigma_sample];
	double accept_prob = -1.0;


	Vector3d orienting_normal2;

	Ray tmp_ray = ray;
	Vector3d dir;
	Intersection out;

	double sampleDir = 1.0;
	double samplePDF = 1.0;



	//0:全域一様サンプリング 
	//1:重点サンプリング(Disk分布)
	//2:重点サンプリング(Gauss分布)
	int samplingType = 2;	

	const int sampleN = 40;
	int err = 1;

	out.hitpoint = hitpoint;
	out.object_id = in_object->id;


	double prob = 0.01;
	double rmax = -log(prob) / directional_dipole.sigma_t.Min();
	double sigmaTr = directional_dipole.sigma_tp.Min();
	//double sigmaTr = directional_dipole.sigma_tr.Min();

	//http://www.cs.cornell.edu/~ags/documents/BSSRDF_Notes.pdf
	//I use n n/σt where n ranges from 1 to 4 or 5.
	//rmax = (1.0 + rnd->next01()*4.0) / directional_dipole.sigma_t.Max();

	if ( samplingType == 1 )
	{
		sigmaTr = directional_dipole.sigma_tp.Min();
		//sigmaTr = directional_dipole.sigma_tp.Max();
		//sigmaTr = directional_dipole.sigma_tr.Min();
		//sigmaTr = directional_dipole.sigma_t.Min();
		double skipRatio = 1.0e-3;
		rmax = 1.5*sqrt(log(skipRatio) / -sigmaTr);
	}

	const double sigmaTrCoef = 0.001;
	if ( samplingType == 2 )
	{
		//sigmaTr = directional_dipole.sigma_tr.Min();
		sigmaTr = directional_dipole.sigma_tp.Min();
		//sigmaTr = directional_dipole.sigma_t.Max();
		double skipRatio = 1.0e-4;
		double x = participating_media.scatteringColor.Min();
		double y = 2.0-1.0/(1.0+exp(-4.0*(x-1.6)));
		if (participating_media.absorbing < 0.001)
		{
			y *= 1.3 - exp(participating_media.absorbingColor.Min()*1000.0);
		}
		rmax = 1.25*sqrt(log(skipRatio) / -sigmaTr)*y*participating_media.r_max;
	}
	//if ( rnd->next01() < 0.05 )
	//{
	//	samplingType = 0;
	//}

	const double ln = (in_object->boundingBox.max() - in_object->boundingBox.min()).length();
	//fprintf(stderr, "r_max %f ln %f\n", r_max, ln);
	if (rmax > ln) rmax = ln*0.3;

	const double areaPdf = 1.0 / const_cast<Entity*>(in_object)->area;
	const double disk_areaPdf = 1.0 / (rmax*rmax*PS_PI);

	const double step = 0.5*rmax / sampleN;
	for (int i = sampleN; i >= 1; i--)
	{
		err = 1;
		Vector3d ww, u, v;

		if ( samplingType < 1 || samplingType > 2)
		{
			Vector3d samplePoint = const_cast<Entity*>(in_object)->randomPoint(&orienting_normal2);
			out.hitpoint.position = samplePoint;
			out.hitpoint.normal = orienting_normal2;

			samplePDF = areaPdf;
		}else
		{
			OrthonormalBasis(orienting_normal, ww, u, v);
			if (samplingType == 1)
			{
				Vector3d pSample = rmax* uniformSampleDisk(rnd->next01(), rnd->next01());

				Vector3d p = hitpoint.position + pSample.x*u + pSample.y*v;

				double halfProbeLength = i *step;
				Ray tmp_ray(p - orienting_normal*halfProbeLength, orienting_normal);
				if (!intersect_scene(tmp_ray, &out, depth, in_object->id))
				{
					continue;
				}
				orienting_normal2 = out.hitpoint.normal;
				samplePDF = disk_areaPdf;
			}
			if (samplingType == 2)
			{
				Vector3d pSample = gaussianSample2D(rnd->next01(), rnd->next01(), sigmaTr*sigmaTrCoef, rmax);

				Vector3d p = hitpoint.position + pSample.x*u + pSample.y*v;

				double halfProbeLength = i*step;
				Ray tmp_ray(p - orienting_normal*halfProbeLength, orienting_normal);
				if (!intersect_scene(tmp_ray, &out, depth, in_object->id))
				{
					continue;
				}
				orienting_normal2 = out.hitpoint.normal;
				samplePDF = gaussianSample2DPdf(Vector3d(0, 0, 0), pSample, Vector3d(0, 0, 1), sigmaTr*sigmaTrCoef, rmax);
			}
		}
		//fprintf(stderr, "area %f\n", const_cast<Entity*>(in_object)->area);
		if (dot(hitpoint.normal, orienting_normal) < 0)
		{
			orienting_normal2 = -1.0*orienting_normal2;
		}

		if ( samplingType < 1 || samplingType > 2)
		{
			double dot12 = dot(orienting_normal2, orienting_normal);
			if (dot12 < 1.0e-6)
			{
				continue;
			}
		}
		//if ((hitpoint.position - out.hitpoint.position).length() > rmax)
		//{
		//	continue;
		//}

		OrthonormalBasis(orienting_normal2, ww, u, v);
		//半球重点サンプリング
		const double r1 = PS_TWOPI * rnd->next01();
		const double r2 = rnd->next01(), r2s = sqrt(r2);
		dir = normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + ww * sqrt(1.0 - r2)));

		sampleDir = dot(orienting_normal2, dir) * PS_INV_PI;
		if (sampleDir > 0)
		{
			err = 0;
			break;
		}
	}

	if (err)
	{
		fprintf(stderr, "==========================\n");
		return ZERO();
	}


	// 交差位置の法線（物体からのレイの入出を考慮）
	Color result;
	Vector3d sp, sn, sw;
	Vector3d n;
	Vector3d x;
	Vector3d w;

	//入口
	w = -1.0*ray.dir;
	x = hitpoint.position;
	n = orienting_normal;

	//出口
	sw = dir;
	sp = out.hitpoint.position;
	sn = orienting_normal2;


	const double cos_wi_n = dot(sn, sw);

	if (cos_wi_n < 0.0)
	{
		fprintf(stderr, "==========================\n");
		//return  radiance(0, env, ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
		return ZERO();
	}

	//fprintf(stderr, "η eta %f\n", eta);
	// compute Fresnel transmittance at point of emergence
	const double T21 = std::max(1.0 - fresnel(dot(w, n), eta), 0.0);
	const double T12 = std::max(1.0 - fresnel(cos_wi_n, eta), 0.0);
	//fprintf(stderr, "η eta %f T21 %f T12 %f\n", eta, T21, T12);

	// Russian roulette (note that accept_prob can be any value in (0, 1))
	//accept_prob = exp(-(sp - x).length() * directional_dipole.min_sigma_tr);
	//if (accept_prob < 0.0) accept_prob = 0.0;
	//else if (accept_prob > 1.0) accept_prob = 1.0;
	
	//accept_prob = 1.0 - exp(-(sp - x).length() * directional_dipole.min_sigma_tr);
	//if (accept_prob < 0.0) accept_prob = 0.0;
	//else if (accept_prob > 1.0) accept_prob = 1.0;

	//accept_prob = 1.0 - exp(-participating_media.transmittance*(sp - x).length());
	//if (accept_prob < 0.0) accept_prob = 0.0;
	//else if (accept_prob > 1.0) accept_prob = 1.0;

	//double σt = directional_dipole.sigma_t.Max();

	//accept_prob = exp(-σt*(sp - x).length());
	//if (accept_prob < 0) accept_prob = 0;
	//else if (accept_prob > 1.0) accept_prob = 1.0;

	//accept_prob = exp(-(sp - x).length() * directional_dipole.sigma_a.Max());
	//if (accept_prob < 0.0) accept_prob = 0.0;
	//else if (accept_prob > 1.0) accept_prob = 1.0;


	//fprintf(stderr, "accept_prob %.16f\n", accept_prob);
	//accept_prob = 1.0;
	//accept_prob = 0.0;
	//accept_prob = 1.0 - (directional_dipole.sigma_a / directional_dipole.sigma_t).Max();
	accept_prob = 0.5;

	const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

	if (!into ) accept_prob = 0;
	//if ( participating_media.level == -2 )
	//{
	//	accept_prob = 1;
	//}

	int sc = 0;
	if (rnd->next01() < accept_prob )
	{
		sc = 1;
		//式(1)の計算
		//result.x = bssrdf(directional_dipole, sp, sn, sw, x, n, w, 0) + bssrdf(directional_dipole, x, n, w, sp, sn, sw, 0);
		//result.y = bssrdf(directional_dipole, sp, sn, sw, x, n, w, 1) + bssrdf(directional_dipole, x, n, w, sp, sn, sw, 1);
		//result.z = bssrdf(directional_dipole, sp, sn, sw, x, n, w, 2) + bssrdf(directional_dipole, x, n, w, sp, sn, sw, 2);
		//result = 0.5*result;

		//result.x = bssrdf(directional_dipole, x, n, w, sp, sn, sw, 0);
		//result.y = bssrdf(directional_dipole, x, n, w, sp, sn, sw, 1);
		//result.z = bssrdf(directional_dipole, x, n, w, sp, sn, sw, 2);

		//result.x = bssrdf(directional_dipole, sp, sn, sw, x, n, w, 0);
		//result.y = bssrdf(directional_dipole, sp, sn, sw, x, n, w, 1);
		//result.z = bssrdf(directional_dipole, sp, sn, sw, x, n, w, 2);

		int multiscattering = 1;
		double probability_multiscattering = 0.5;
		//probability_multiscattering = 1.0;

		if ( participating_media.level == -2 )
		{
			probability_multiscattering = 1;
			//fprintf(stderr, "------------------------------------\n");
		}
		//if (!into )
		//{
		//	probability_multiscattering = 1;
		//}
		if (rnd->next01() < probability_multiscattering)
		{
			multiscattering = 1;
			double res[3];
			for (int j = 0; j < 3; j++)
			{
				//res[j] = bssrdf(directional_dipole, x, n, w, sp, sn, sw, j);
				res[j] = bssrdf(directional_dipole, sp, sn, sw, x, n, w, j);
				//res[j] = 0.5*(bssrdf(directional_dipole, sp, sn, sw, x, n, w, j) + bssrdf(directional_dipole, x, n, w, sp, sn, sw, j));
			}
			result.x = res[0];
			result.y = res[1];
			result.z = res[2];

			weight = RGB2Spectrum(hitpoint.material.color, wavelength) * T12 * RGB2Spectrum(result, wavelength) * cos_wi_n*T21
				/ samplePDF
				/ sampleDir
				/ accept_prob
				/ russian_roulette_probability / (1.0 - scattering_probability)
				/ probability_multiscattering
				;

			incoming_radiance = radiance(0, env, Ray(sp, sw), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);

			//fprintf(stderr, "-----------------------------%f %f %f\n", incoming_radiance.x, incoming_radiance.y, incoming_radiance.z);
			//fprintf(stderr, "                             %f %f %f\n", weight.x, weight.y, weight.z);
		}
		else
		{
			multiscattering = 0;
			//sigle scattering
			int sc_max = 1;
			weight = RGB2Spectrum(hitpoint.material.color, wavelength)
				/ accept_prob
				/ russian_roulette_probability / (1.0 - scattering_probability) 
				/ (1.0 - probability_multiscattering)
				;

			double totalDist = 0.0;

			incoming_radiance = radiance_through_media(env, hitpoint, Ray(hitpoint.position - hitpoint.normal*0.0001, ray.dir), rnd, depth + 1, in_object, 0, sc_max, russian_roulette_probability, scattering_probability, wavelength, totalDist, 1 * nextEventEstimation, participatingMedia);
		}

		if (multiscattering)
		{
			//out.hitpoint.material.color = participating_media.albedoColor;
			int direct_light_sample = 0;
			if (nextEventEstimation && env->light_list.list.size())
			{
				// 直接光の評価
				//光源と光源をサンプリングしてそのサンプリング位置から現在位置を結ぶ経路間で放射輝度を求める
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);

					//weightが係って返ってくる
					Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, in_object, out.hitpoint, orienting_normal2, dir, depth, wavelength) / light_probability;

					direct_light = direct_light + dl;
					if (SIZE(dl) > 1.0e-16)
					{
						direct_light_sample++;
					}
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light	*T12 * RGB2Spectrum(result, wavelength) * cos_wi_n*T21
					/ samplePDF
					/ sampleDir
					/ accept_prob
					/ russian_roulette_probability / (1.0 - scattering_probability)	
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
						Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, in_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

						direct_light = direct_light + dl;
					}
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light	*T12 * RGB2Spectrum(result, wavelength) * cos_wi_n*T21
					/ samplePDF
					/ sampleDir
					/ accept_prob
					/ russian_roulette_probability / (1.0 - scattering_probability)				
					;
			}


			//light intensity weight
			double lightIntensity = 1.0;
			if (hitpoint.bump && bump)
			{
				lightIntensity = std::max(0.1, dot(dir, out.hitpoint.bump_new_normal)) / dot(dir, out.hitpoint.normal);
			}
			weight = weight*lightIntensity;


			if (0 * nextEventEstimation && !direct_light_sample)
			{
				Intersection lid;
				if (intersect_scene(Ray(sp, sw), &lid, depth))
				{
					const Entity *now_object2 = EntList->List[lid.object_id];

					if (now_object2->light)
					{
						direct_light = RGB2Spectrum(lid.hitpoint.material.emission, wavelength);

						//平行光源からの寄与があるか
						const Light& light = env->light_list.list[now_object2->light_id];
						if (light.parallel_light)
						{
							direct_light = parallel_light_direct(hitpoint, lid, dir, light, wavelength);
						}
						if (light.spot_light)
						{
							direct_light = spot_light_direct(hitpoint, lid, dir, light, wavelength);
						}

						//関与媒質による減衰を考慮
						if (participatingMedia)
						{
							direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance, wavelength) * direct_light;
						}
						//後でweightを掛ける必要があるためincoming_radianceに加算しておく
						incoming_radiance = incoming_radiance + direct_light;
						direct_light = ZERO();
					}
				}
			}
		}
#if OVER_RADIANCE_CUT2
#ifndef SPECTRUM_USE
		if ((incoming_radiance*weight).length() > OVER_RADIANCE_CUT2)
		{
			return ZERO();
		}
#else
		if (Spectrum2RGB(incoming_radiance*weight).length() > OVER_RADIANCE_CUT1) return ZERO();
#endif
#endif

	}
	else
	if(1){
		if (hitpos_material.reflection_type == REFLECTION_TYPE_SSS_REFRACTION)
		{
			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

			Spectrum Tr = ONE();
			if (into)
			{
			}
			else
			{
				//出て行く方向だから直前からの距離分だけ減衰している
				Tr = participating_media.transmittanceRatio(hitpoint.distance, wavelength);
			}

			Ray tmp_ray = ray;
			tmp_ray.org = hitpoint.position;

			//fprintf(stderr, "hitpoint.material.roughness %f\n", hitpoint.material.roughness);
			if (hitpoint.material.roughness)
			{
				//abort();
				double diffuse_probability = hitpoint.material.roughness;

				double η = 1.495;//polypropylene;
				double m_invη2 = 1.0 / (η*η);

				if (rnd->next01() < diffuse_probability)
				{
					IntersectionPos hitp = hitpoint;

					hitp.material.reflection_type = REFLECTION_TYPE_DIFFUSE;
					hitp.material.color = hitp.material.color*participating_media.albedoColor;

					//pdf(=diffuse_probability)で割らないのは反射率にdiffuse_probabilityが掛かって相殺するため radiance.hの●●●に詳細説明
					weight = Tr*m_invη2*(1.0 - fresnelDiffuseReflectance(η))/(1.0 - accept_prob);

					incoming_radiance = BRDF_Diffuse(0, env, rnd, tmp_ray, hitp, in_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
#if OVER_RADIANCE_CUT1
#ifndef SPECTRUM_USE
					if ((incoming_radiance*weight).length() > OVER_RADIANCE_CUT1)
					{
						return ZERO();
					}
#else
					if (Spectrum2RGB(incoming_radiance*weight).length() > OVER_RADIANCE_CUT1) return ZERO();
#endif
#endif
				}
				else
				{
					//abort();
					IntersectionPos hitp = hitpoint;
					hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION;
					hitp.material.refractive_index = η;
					hitp.material.color = hitp.material.color;

					//pdf(=diffuse_probability-1)で割らないのは反射率に(diffuse_probability-1)が掛かって相殺するため radiance.hの●●●に詳細説明
					weight = Tr / (1.0 - accept_prob);

					incoming_radiance = BRDF_Refraction(env, rnd, tmp_ray, hitp, in_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);
				}
			}
			else
			{

				IntersectionPos hitp = hitpoint;
				hitp.material.reflection_type = REFLECTION_TYPE_REFRACTION;
				hitp.material.color = hitpoint.material.color;

				weight = Tr/(1.0 - accept_prob);

				incoming_radiance = BRDF_Refraction(env, rnd, tmp_ray, hitp, in_object, orienting_normal, russian_roulette_probability, scattering_probability, depth, wavelength, nextEventEstimation, participatingMedia);

#if OVER_RADIANCE_CUT3
#ifndef SPECTRUM_USE
				if ((incoming_radiance*weight).length() > OVER_RADIANCE_CUT3)
				{
					return ZERO();
				}
#else
				if (Spectrum2RGB(incoming_radiance*weight).length() > OVER_RADIANCE_CUT1) return ZERO();
#endif
#endif
			}
			direct_light = ZERO();
		}
	}

rtn:;


	//関与媒質を考慮する場合は光の減衰を考慮する
	if (participatingMedia)
	{
		//環境内（物体の外側)
		//現在の交点ではレイは物体の中に入ろうとする方向なら物体の外側から
		if (inside_dir)
		{
			direct_light = transmittance_ratio*direct_light;
			emission = transmittance_ratio * emission;
			incoming_radiance = transmittance_ratio * incoming_radiance;
		}
	}

	return RGB2Spectrum(env->Environment_Light, wavelength) + RGB2Spectrum(emission, wavelength) + direct_light + MULTIPLY(weight, incoming_radiance);
}
};
#endif