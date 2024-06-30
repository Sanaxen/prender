#ifndef __BRDF_H_
#define __BRDF_H_

#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "intersection.h"
#include "random.h"

#include "entity.h"
#include "scene_env.h"
#include "Spectrum.h"

namespace prender {

	inline Spectrum BRDF_Phong(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		const Material& hitpos_material = hitpoint.material;
		Color emission = ZERO();					//���Ȕ���
		Spectrum direct_light = ZERO();				//���ڌ�(����)
		Spectrum direct_light_refraction = ZERO();	//���ڌ�(����)
		Spectrum incoming_radiance = 0.0;
		Spectrum weight = 1.0;

		// no absorption inside glass
		bool inside_dir = const_cast<Entity*>(now_object)->isInsideDir2(ray, hitpoint);

		emission = hitpos_material.emission;		// Init color to emissive value

		Spectrum transmittance_ratio = ONE();
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

		{
			Vector3d dir = orienting_normal;

			// �g�U���˂ł�orienting_normal�̕�������Ƃ������K�������(w, u, v)����邪�A
			// Phong�ł͊��S���˕�������Ƃ���
			// ���̊��ɑ΂��锼�����Ŏ��̃��C���΂��B
			Vector3d w, u, v;
			OrthonormalBasis(normalize(ray.dir - orienting_normal*2.0*dot(orienting_normal, ray.dir)), w, u, v);
			//do
			{
				//���S���˕�������ɃT���v�����O����
				const double ph = PS_TWOPI*rnd->next01();
#if 10
				double cost = pow(rnd->next01(), 1.0 / (hitpos_material.brdfParameter.phong_brdf.specular_exponent + 1.0));
				if (cost > 1.0) cost = 1.0;
				if (cost < -1.0) cost = -1.0;
				const double sint = sqrt(1.0 - cost*cost);
#else
				double cost = pow(rnd->next01(), 1.0/(hitpos_material.specular_exponent+1.0));
				if (cost > 1.0) cost = 1.0;
				if (cost < -1.0) cost = -1.0;
				const double th = acos(cost);
				const double sint = sin(th);
#endif

				// �d�_�I�T���v�����O
				dir = normalize((
					u * cos(ph) * sint +
					v * sin(ph) * sint +
					w * cost));
			}
			//while( dot(orienting_normal, dir) < 0.0 );

			if (dot(orienting_normal, dir) < 0.0)
			{
				return ZERO();
			}
			// �d�݁Bbrdf * cos�� / pdf(��)�B
			/*��=���˕����Ɗ��S���˕����Ƃ̊p�x
			pdf = ((n + 1) / 2��) cos(��)^n
			brdf = ((n + 2) / 2��) cos(��)^n  ---> 		brdf = ((n + 2) / 2��) cos(��)^n / cos(��) ?


			brdf / pdf = ((n + 2) / (n + 1))
			*/
			const double ww = (hitpos_material.brdfParameter.phong_brdf.specular_exponent + 2) / (hitpos_material.brdfParameter.phong_brdf.specular_exponent + 1);

			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
			weight = RGB2Spectrum(hitpos_material.specular, wavelength) * ww /**dot(orienting_normal, dir)*/ / russian_roulette_probability / (1.0 - scattering_probability);

#if 0
			if (nextEventEstimation && env->light_list.list.size())
			{
				double light_probability = 1.0;
				//���ڌ��T���v�����O����
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
					direct_light = direct_light + direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				// direct_light �͂��łɔ��˗�����Z�ς݂Ȃ̂ŁAweight���|����K�v�͂Ȃ�
				//�Ȃ̂Ō��weight*incoming_radiance�����邪incoming_radiance�ɂ͉��Z���Ȃ�
				//�܂�@direct_light + weight*incoming_radiance
			}
#else

			//light intensity weight
			double lightIntensity = 1.0;
			if (hitpoint.bump)
			{
				lightIntensity = std::max(0.1, dot(dir, hitpoint.bump_new_normal)) / dot(dir, hitpoint.normal);;
			}
			weight = weight*lightIntensity;

#if 10
			if (nextEventEstimation && env->light_list.list.size())
			{
				// ���ڌ��̕]��
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);

					//weight���W���ĕԂ��Ă���
					Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

					direct_light = direct_light + dl;
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
			}
			if (!nextEventEstimation && env->light_list.list.size())
			{
				// ���ڌ��̕]��
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
					if (light.parallel_light || light.spot_light || light.infinity_light)
					{

						//weight���W���ĕԂ��Ă���
						Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

						direct_light = direct_light + dl;
					}
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
			}
#else
			if (nextEventEstimation)
			{
				// ���ڌ��̕]��
				Intersection lid;
				if (intersect_scene(Ray(hitpoint.position, dir), &lid, depth))
				{
					const Entity *now_object2 = EntList->List[lid.object_id];
					if (now_object2->light)
					{
						direct_light = RGB2Spectrum(lid.hitpoint.material.emission, wavelength);

						//���s��������̊�^�����邩
						const Light& light = env->light_list.list[now_object2->light_id];
						if (light.parallel_light)
						{
							direct_light = parallel_light_direct(hitpoint, lid, dir, light, wavelength);
						}
						if (light.spot_light)
						{
							direct_light = spot_light_direct(hitpoint, lid, dir, light, wavelength);
						}
						//�֗^�}���ɂ�錸�����l��
						if (participatingMedia)
						{
							direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance, wavelength) * direct_light;
						}
					}
				}
				// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
				incoming_radiance = incoming_radiance + direct_light;
				direct_light = ZERO();
			}
#endif
#endif

		}

		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia)
		{
			//�����i���̂̊O��)
			//���݂̌�_�ł̓��C�͕��̂̒��ɓ��낤�Ƃ�������Ȃ畨�̂̊O������
			if (inside_dir)
			{
				direct_light = transmittance_ratio*direct_light;
				emission = transmittance_ratio * emission;
				incoming_radiance = transmittance_ratio * incoming_radiance;
			}
		}
		return RGB2Spectrum(env->Environment_Light, wavelength) + RGB2Spectrum(emission, wavelength) + direct_light + MULTIPLY(weight, incoming_radiance);
	}

	inline Spectrum BRDF_Ward(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		const Material& hitpos_material = hitpoint.material;
		Color emission = ZERO();					//���Ȕ���
		Spectrum direct_light = ZERO();				//���ڌ�(����)
		Spectrum direct_light_refraction = ZERO();	//���ڌ�(����)
		Spectrum incoming_radiance = 0.0;
		Spectrum weight = 1.0;

		// no absorption inside glass
		bool inside_dir = const_cast<Entity*>(now_object)->isInsideDir2(ray, hitpoint);

		emission = hitpos_material.emission;		// Init color to emissive value

		Spectrum transmittance_ratio = ONE();
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

		{

			//fprintf(stderr, "%f %f %f,%f,%f\n", hitpos_material.ward_brdf.alp_x, hitpos_material.ward_brdf.alp_y, hitpos_material.ward_brdf.specular.x, hitpos_material.ward_brdf.specular.y, hitpos_material.ward_brdf.specular.z);

			const double alpha_x = hitpos_material.brdfParameter.ward_brdf.alp_x;
			const double alpha_y = hitpos_material.brdfParameter.ward_brdf.alp_y;

			const Vector3d in = -1.0 * ray.dir;
			Vector3d halfv;
			Vector3d dir;

			// orienting_normal�̕�������Ƃ������K�������(w, u, v)�����B���̊��ɑ΂��锼�����Ŏ��̃��C���΂��B
			Vector3d w, u, v;
			OrthonormalBasis(orienting_normal, w, u, v);


#if 10	//�d�_�T���v�����O�iNotes on the Ward BRDF)
			//Bruce Walter
			//Technical Report PCG-05-06
			//Cornell Program of Computer Graphics
			//April 29, 2005�j

			// �d�_�T���v�����O����B�������A���������n�[�t�x�N�g���������ł�
			// ���˕����������O�ɏo�Ă��܂��̂ł��������ꍇ�͊��p����B
			//do 
			{
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
			}
			//while (dot(orienting_normal, dir) < 0.0);
			if (dot(orienting_normal, dir) < 0.0)
			{
				return ZERO();
			}

			// �d�݁Bbrdf * cos�� / pdf(��)���܂Ƃ߂����̂ɂȂ��Ă���B
			const double ww = dot(halfv, in) * pow(dot(halfv, orienting_normal), 3) *
				sqrt(dot(dir, orienting_normal) / dot(in, orienting_normal));

			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
			weight = RGB2Spectrum(hitpos_material.specular, wavelength) * ww / russian_roulette_probability / (1.0 - scattering_probability);
#else	//��l�T���v�����O
			while (1)
			{
				// ���������n�[�t�x�N�g�����˕����������O�ɏo�Ă��܂��ꍇ�͂�蒼���B
				do {
					//��l�T���v�����O
					const double t = rnd->next01();
					const double s = PS_TWOPI * rnd->next01();
					const double tt = 1 - t;
					const double tt2 = sqrt(1.0 - tt*tt);

					halfv = normalize((
						u * cos(s) * tt2 +
						v * sin(s) * tt2 +
						w * tt));

					dir = normalize(2.0 * dot(in, halfv) * halfv - in);
				} while (dot(orienting_normal, dir) < 0.0);

				double ww0 = pow(dot(halfv, orienting_normal), 2.0);
				//if ( fabs(ww0) < 1.0e-14 )
				//{
				//	continue;
				//}
				double ww1 = pow(dot(halfv, u) / alpha_x, 2.0) + pow(dot(halfv, v) / alpha_y, 2.0) / ww0;
				double ww2 = 4.0*PS_PI*alpha_x*alpha_y*sqrt(dot(dir, orienting_normal)* dot(in, orienting_normal));

				//if ( fabs(ww2) < 1.0e-14 )
				//{
				//	continue;
				//}
				const double ww = exp(-ww1) / ww2;
				incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
				weight = RGB2Spectrum(hitpos_material.specular, wavelength) * ww * dot(w, dir)*PS_TWOPI / russian_roulette_probability / (1.0 - scattering_probability);
			}
#endif

			//light intensity weight
			double lightIntensity = 1.0;
			if (hitpoint.bump)
			{
				lightIntensity = std::max(0.1, dot(dir, hitpoint.bump_new_normal)) / dot(dir, hitpoint.normal);;
			}
			weight = weight*lightIntensity;

#if 10
			if (nextEventEstimation && env->light_list.list.size())
			{
				// ���ڌ��̕]��
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);

					//weight���W���ĕԂ��Ă���
					Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

					direct_light = direct_light + dl;
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
			}
			if (!nextEventEstimation && env->light_list.list.size())
			{
				// ���ڌ��̕]��
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
					if (light.parallel_light || light.spot_light || light.infinity_light)
					{

						//weight���W���ĕԂ��Ă���
						Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

						direct_light = direct_light + dl;
					}
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
			}
#else
			if (0*nextEventEstimation)
			{
				// ���ڌ��̕]��
				Intersection lid;
				if (intersect_scene(Ray(hitpoint.position, dir), &lid, depth))
				{
					const Entity *now_object2 = EntList->List[lid.object_id];
					if (now_object2->light)
					{
						direct_light = RGB2Spectrum(lid.hitpoint.material.emission, wavelength);

						//���s��������̊�^�����邩
						const Light& light = env->light_list.list[now_object2->light_id];
						if (light.parallel_light)
						{
							direct_light = parallel_light_direct(hitpoint, lid, dir, light, wavelength);
						}
						if (light.spot_light)
						{
							direct_light = spot_light_direct(hitpoint, lid, dir, light, wavelength);
						}
						//�֗^�}���ɂ�錸�����l��
						if (participatingMedia)
						{
							direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance, wavelength) * direct_light;
						}
					}
				}
				// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
				incoming_radiance = incoming_radiance + direct_light;
				direct_light = ZERO();
			}
#endif
		}

		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia)
		{
			//�����i���̂̊O��)
			//���݂̌�_�ł̓��C�͕��̂̒��ɓ��낤�Ƃ�������Ȃ畨�̂̊O������
			if (inside_dir)
			{
				direct_light = transmittance_ratio*direct_light;
				emission = transmittance_ratio * emission;
				incoming_radiance = transmittance_ratio * incoming_radiance;
			}
		}
		return RGB2Spectrum(env->Environment_Light, wavelength) + RGB2Spectrum(emission, wavelength) + direct_light + MULTIPLY(weight, incoming_radiance);
	}


	inline Spectrum BRDF_Diffuse(int* hit_direct_light, const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		const Material& hitpos_material = hitpoint.material;
		Color emission = ZERO();					//���Ȕ���
		Spectrum direct_light = ZERO();				//���ڌ�(����)
		Spectrum direct_light_refraction = ZERO();	//���ڌ�(����)
		Spectrum incoming_radiance = 0.0;
		Spectrum weight = 1.0;

		// no absorption inside glass
		bool inside_dir = const_cast<Entity*>(now_object)->isInsideDir2(ray, hitpoint);

		emission = hitpos_material.emission;		// Init color to emissive value

		if (hit_direct_light != NULL)
		{
			*hit_direct_light = 0;
		}

		Spectrum transmittance_ratio = ONE();
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

		{

			//���ʌ����̗����͌���Ȃ��悤�ɂ��鏈��
			if (now_object->light && now_object->type == ENTITY_TYPE_UVPLANE)
			{
				if (dot(hitpoint.normal, ray.dir) >= 0.0)
				{
					return ZERO();
				}
			}

			// ���ڌ����v�Z���Ă���̂�eye ���璼�� light �� hit �����ꍇ�̂݁Aemission ��Ԃ���
			if (nextEventEstimation)
			{
				if (now_object->light)
				{
					if (depth == 0)
					{
						if (hit_direct_light != NULL)
						{
							*hit_direct_light = 1;
						}

						//return RGB2Spectrum( hitpos_material.emission, wavelength);
						return RGB2Spectrum(transmittance_ratio * hitpos_material.emission, wavelength) / russian_roulette_probability / (1.0 - scattering_probability);
					}
					else
					{
						return ZERO();
					}
				}
			}
			else
			{
				if (now_object->light)
				{
					//const Light& light = env->light_list.list[now_object->light_id];
					//if (light.parallel_light || light.spot_light)
					//{
					//	return ZERO();
					//}

					return RGB2Spectrum(transmittance_ratio * hitpos_material.emission, wavelength) / russian_roulette_probability / (1.0 - scattering_probability);

				}
			}



			// orienting_normal�̕�������Ƃ������K�������(w, u, v)�����B���̊��ɑ΂��锼�����Ŏ��̃��C���΂��B
			Vector3d w, u, v;
			OrthonormalBasis(orienting_normal, w, u, v);

			//if (!nextEventEstimation)
			//{
			//	if (now_object->light)
			//	{
			//		const Light& light = env->light_list.list[now_object->light_id];
			//		//���s�����������ꍇ
			//		if (light.parallel_light)
			//		{
			//			//���˕����͈��
			//			Vector3d dir = light.dir*-1.0;
			//			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
			//			weight = RGB2Spectrum(hitpos_material.color, wavelength)*(dot(w, dir) / PS_PI) / russian_roulette_probability / (1.0 - scattering_probability);
			//			goto FIN;
			//		}
			//		if (light.spot_light)
			//		{
			//			//���˕����͈��
			//			Vector3d dir = normalize(hitpoint.position - ((Sphere*)light.light)->position)*-1.0;
			//			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
			//			weight = RGB2Spectrum(hitpos_material.color, wavelength)*(dot(w, dir) / PS_PI) / russian_roulette_probability / (1.0 - scattering_probability);
			//			goto FIN;
			//		}
			//	}
			//}

#if 10
			// �R�T�C�������g�����d�_�I�T���v�����O
			const double r1 = PS_TWOPI * rnd->next01();
			const double r2 = rnd->next01(), r2s = sqrt(r2);
			Vector3d dir = normalize((
				u * cos(r1) * r2s +
				v * sin(r1) * r2s +
				w * sqrt(1.0 - r2)));

			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
			// �����_�����O�������ɑ΂��郂���e�J�����ϕ����l����ƁAoutgoing_radiance = weight * incoming_radiance�B
			// �����ŁAweight = (��/��) * cos�� / pdf(��) / R �ɂȂ�B
			// ��/�΂͊��S�g�U�ʂ�BRDF�Ńς͔��˗��Acos�Ƃ̓����_�����O�������ɂ�����R�T�C�����Apdf(��)�̓T���v�����O�����ɂ��Ă̊m�����x�֐��B
			// R�̓��V�A�����[���b�g�̊m���B
			// ���A�R�T�C�����ɔ�Ⴕ���m�����x�֐��ɂ��T���v�����O���s���Ă��邽�߁Apdf(��) = cos��/��
			// ����āAweight = ��/ R�B
			weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability / (1.0 - scattering_probability);

			//if ( depth+1 >= env->Depth && rnd->next01() < russian_roulette_probability)
			//{
			//	incoming_radiance = ONE()*dot(w, dir)*PS_TWOPI;
			//	//direct_light = Color(0.01, 0.01, 0.01);
			//}
#else
			//��l�T���v�����O
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

			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
			// �����_�����O�������ɑ΂��郂���e�J�����ϕ����l����ƁAoutgoing_radiance = weight * incoming_radiance�B
			// �����ŁAweight = (��/��) * cos�� / pdf(��) / R �ɂȂ�B pdf(��)=1/(2��)
			// ��/�΂͊��S�g�U�ʂ�BRDF�Ńς͔��˗��Acos�Ƃ̓����_�����O�������ɂ�����R�T�C����
			weight = RGB2Spectrum(hitpos_material.color / PS_PI, wavelength)*dot(w, dir)*PS_TWOPI / russian_roulette_probability;
#endif
			//light intensity weight
			double lightIntensity = 1.0;
			if (hitpoint.bump)
			{
				lightIntensity = std::max(0.1, dot(dir, hitpoint.bump_new_normal)) / dot(dir, hitpoint.normal);
			}
			weight = weight*lightIntensity;

			//shadow�����̕��ʂɂԂ����Ă���
			if (now_object->type == ENTITY_TYPE_PLANE || now_object->type == ENTITY_TYPE_UVPLANE)
			{
				//Angurer IBL�ŉe�����
				bool shadow_plane = false;
				if (now_object->type == ENTITY_TYPE_PLANE && ((Plane*)now_object)->shadow)
				{
					shadow_plane = true;
				}
				if (now_object->type == ENTITY_TYPE_UVPLANE && ((UVPlane*)now_object)->shadow)
				{
					shadow_plane = true;
				}
				if (shadow_plane)
				{
					// ���ڌ��̕]��
					Intersection lid;

					//���˃��C�͉����Ɍ��������H
					if (intersect_scene(Ray(hitpoint.position + PS_EPS8*dir, dir), &lid, depth))
					{
						Entity *now_object2 = EntList->List[lid.object_id];

						//IBL�ɒ��ڂԂ�����
						if (now_object2->material()->IBL())
						{
							//���̂܂܊ђʂ���B
							return radiance(0, env, Ray(hitpoint.position + PS_EPS8*ray.dir, ray.dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);;
						}
					}
				}
			}

			if (nextEventEstimation && env->light_list.list.size())
			{
				// ���ڌ��̕]��
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);

					//weight���W���ĕԂ��Ă���
					Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

					direct_light = direct_light + dl;
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
				//if ( shadow_ray_hit_count ) fprintf(stderr, "shadow_ray_hit_count %d\n", shadow_ray_hit_count);

			}

			if (!nextEventEstimation && env->light_list.list.size())
			{
				// ���ڌ��̕]��
				//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
				double light_probability = 1.0;
				for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
				{
					const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
					if (light.parallel_light || light.spot_light || light.infinity_light)
					{

						//weight���W���ĕԂ��Ă���
						Spectrum dl = direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

						direct_light = direct_light + dl;
					}
				}
				direct_light = direct_light / SHADOW_RAY_SAMPLING;
				direct_light = direct_light / russian_roulette_probability / (1.0 - scattering_probability);
				//if ( shadow_ray_hit_count ) fprintf(stderr, "shadow_ray_hit_count %d\n", shadow_ray_hit_count);

			}


			// direct_light �͂��łɔ��˗�����Z�ς݂Ȃ̂ŁAweight���|����K�v�͂Ȃ�
			//�Ȃ̂Ō��weight*incoming_radiance�����邪incoming_radiance�ɂ͉��Z���Ȃ�
			//�܂�@direct_light + weight*incoming_radiance

#if 0
			//��������r�I�B��ĂĂ���悤�ȏ󋵂Ō����̖ʐς����ɑ傫���ƌ����̈ʒu�T���v�����O��
			//��������̃��C���q�b�g���Ȃ��悤�ȃP�[�X���w�ǂɂȂ��Ă��܂�
			//���̏ꍇ�͍��̏Փ˓_����̃��C���΂��Č����Ƀq�b�g����\�������ꍇ�͈ȉ��̂悤�ɍl������B
			//�����̂悤�Ȍ��ۂ͂��肦��������nextEventEstimation�̌v�Z���Ԉ���Ă��Ă����Ȃ�̂����m��Ȃ�

			//���͂�͂�ԈႢ������(2014.9.10) ���ʌ����̍��W�ϊ��Ō����ʂ͕ϊ����Ă��Ȃ��������߃q�b�g���Ȃ������B
			//nextEventEstimation���w�肵�Ȃ��ƌ����Ƀq�b�g���Ă�����g���[�X�����邽�߂ǂ����̌o�H�Ńq�b�g����\�����傫��������
			//nextEventEstimation���w�肵�Ă���ƌ�������̃g���[�X�͖������߃q�b�g����\���͖w�ǖ�������
			if (nextEventEstimation && shadow_ray_hit_count == 0)
			{
				const bool flg = true;

				// ���ڌ��̕]��
				Intersection lid;
				if (intersect_scene(Ray(hitpoint.position, dir), &lid))
				{
					const Entity *now_object2 = EntList->List[lid.object_id];
					if (now_object2->light)
					{
						if (flg)
						{
							Color emit = env->light_list.list[now_object2->light_id].light->material.emission;
							//���˗�����Z����Ă��Ȃ��̂Ō�Ŕ��˗����|����
							//incoming_radiance =  incoming_radiance + emit / (SHADOW_RAY_SAMPLING+1);

							direct_light = emit / (SHADOW_RAY_SAMPLING + 1) / russian_roulette_probability / (1.0 - scattering_probability);
						}
						else
						{
							const Vector3d light_pos = lid.hitpoint.position;
							const Light& light = env->light_list.list[now_object2->light_id];

							direct_light = eval_direct_radiance(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, light_pos, depth, wavelength);
							direct_light = direct_light / (SHADOW_RAY_SAMPLING + 1) / russian_roulette_probability / (1.0 - scattering_probability);
						}
					}
				}
			}
#endif


		}

	FIN:;
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia)
		{
			//�����i���̂̊O��)
			//���݂̌�_�ł̓��C�͕��̂̒��ɓ��낤�Ƃ�������Ȃ畨�̂̊O������
			if (inside_dir)
			{
				direct_light = transmittance_ratio*direct_light;
				emission = transmittance_ratio * emission;
				incoming_radiance = transmittance_ratio * incoming_radiance;
			}
		}
		return RGB2Spectrum(env->Environment_Light, wavelength) + RGB2Spectrum(emission, wavelength) + direct_light + MULTIPLY(weight, incoming_radiance);

	}

	inline Spectrum BRDF_Specular(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		const Material& hitpos_material = hitpoint.material;
		Color emission = ZERO();					//���Ȕ���
		Spectrum direct_light = ZERO();				//���ڌ�(����)
		Spectrum direct_light_refraction = ZERO();	//���ڌ�(����)
		Spectrum incoming_radiance = ZERO();
		Spectrum weight = 1.0;

		// no absorption inside glass
		bool inside_dir = const_cast<Entity*>(now_object)->isInsideDir2(ray, hitpoint);


		emission = hitpos_material.emission;		// Init color to emissive value

		Spectrum transmittance_ratio = ONE();
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

		const Vector3d dir = normalize(ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir));
		{

			// ���S���ʂȂ̂Ń��C�̔��˕����͌���I�B
			// ���V�A�����[���b�g�̊m���ŏ��Z����̂͏�Ɠ����B
			incoming_radiance = radiance(0, env, Ray(hitpoint.position, dir), rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);

			weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability / (1.0 - scattering_probability);


			if (nextEventEstimation)
			{
				direct_light = direct_light_perfect_refrection(env, hitpoint, Ray(hitpoint.position, dir), depth, wavelength, participatingMedia);

				// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
				incoming_radiance = incoming_radiance + direct_light;
			}
		}

		//light intensity weight
		double lightIntensity = 1.0;
		if (hitpoint.bump)
		{
			lightIntensity = std::max(0.1, dot(dir, hitpoint.bump_new_normal)) / dot(dir, hitpoint.normal);;
		}
		weight = weight*lightIntensity;

		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia)
		{
			//�����i���̂̊O��)
			//���݂̌�_�ł̓��C�͕��̂̒��ɓ��낤�Ƃ�������Ȃ畨�̂̊O������
			if (inside_dir)
			{
				emission = transmittance_ratio * emission;
				incoming_radiance = transmittance_ratio * incoming_radiance;
			}
		}
		return RGB2Spectrum(env->Environment_Light, wavelength) + RGB2Spectrum(emission, wavelength) + MULTIPLY(weight, incoming_radiance);

	}

	inline Spectrum BRDF_Refraction(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		const Material& hitpos_material = hitpoint.material;
		Color emission = ZERO();					//���Ȕ���
		Spectrum direct_light = ZERO();				//���ڌ�(����)
		Spectrum direct_light_refraction = ZERO();	//���ڌ�(����)
		Spectrum incoming_radiance = 0.0;
		Spectrum weight = 1.0;

		// no absorption inside glass
		bool inside_dir = const_cast<Entity*>(now_object)->isInsideDir2(ray, hitpoint);

		emission = hitpos_material.emission;		// Init color to emissive value

		Spectrum transmittance_ratio = ONE();
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

		{

			const Ray reflection_ray = Ray(hitpoint.position, normalize(ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir)));
			//const Ray reflection_ray = Ray(hitpoint.position, normalize(ray.dir - orienting_normal * 2.0 * dot(orienting_normal, ray.dir)));
			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // ���C���I�u�W�F�N�g����o��̂��A����̂�

			// ���˕����̒��ڌ��̕]��
			if (nextEventEstimation)
			{
				direct_light = direct_light_perfect_refrection(env, hitpoint, reflection_ray, depth, wavelength, participatingMedia);
			}

			// Snell�̖@��
			const double nc = 1.0;								// �^��̋��ܗ�
			double nt = hitpos_material.refractive_index;		// �I�u�W�F�N�g�̋��ܗ�(kIor)

			nt = hitpos_material.refractive(wavelength);

			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = dot(ray.dir, orienting_normal);
			const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

			//printf("now_object->refractive_index %f\n", now_object->refractive_index);
			if (cos2t < 0.0)
			{
				// �S����
				incoming_radiance = radiance(0, env, reflection_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);
				weight = RGB2Spectrum(hitpos_material.color, wavelength) / russian_roulette_probability / (1.0 - scattering_probability);

				// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
				incoming_radiance = incoming_radiance + direct_light;

				//light intensity weight
				double lightIntensity = 1.0;
				if (hitpoint.bump)
				{
					lightIntensity = std::max(0.1, dot(reflection_ray.dir, hitpoint.bump_new_normal)) / dot(reflection_ray.dir, hitpoint.normal);
				}
				weight = weight*lightIntensity;

				goto Rtn;
			}

			// ���܂̕���
			const Ray tdir = Ray(hitpoint.position,
				normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(fabs(cos2t)))));
			//const Ray tdir = Ray(hitpoint.position,
			//	normalize(ray.dir * nnt - orienting_normal * /*(into ? 1.0 : -1.0) **/ (ddn * nnt + sqrt(fabs(cos2t)))));

			// Schlick�ɂ��Fresnel�̔��ˌW���̋ߎ����g��
			const double a = nt - nc, b = nt + nc;
			const double R0 = (a * a) / (b * b);

			const double c = 1.0 - (into ? -ddn : dot(tdir.dir, -1.0 * orienting_normal));
			double Re = R0 + (1.0 - R0) * pow(c, 5.0); // ���˕����̌������˂���ray.dir�̕����ɉ^�Ԋ����B�����ɋ��ܕ����̌������˂�������ɉ^�Ԋ����B
			const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // ���C�̉^�ԕ��ˋP�x�͋��ܗ��̈قȂ镨�̊Ԃ��ړ�����Ƃ��A���ܗ��̔�̓��̕������ω�����B
			double Tr = (1.0 - Re) * nnt2; // ���ܕ����̌������܂���ray.dir�̕����ɉ^�Ԋ���

			// ���ȏヌ�C��ǐՂ�������܂Ɣ��˂̂ǂ��炩�����ǐՂ���B�i�����Ȃ��Ǝw���I�Ƀ��C��������j
			// ���V�A�����[���b�g�Ō��肷��B
			double probability = 0.25 + 0.5 * Re;

			if (hitpos_material.reflection_type == REFLECTION_TYPE_REFRACTION_FRESNEL)
			{
				if (into)
				{
					probability = 1.0;
					Re = 1.0;
					Tr = 0.0;
				}
				else
				{
					probability = 0.0;
					Re = 0.0;
					Tr = 1.0;
				}
			}

			// ���ܕ����̒��ڌ��̕]��
			if (nextEventEstimation)
			{
				direct_light_refraction = direct_light_perfect_refrection(env, hitpoint, tdir, depth, wavelength, participatingMedia);
			}

			if (depth > 2) {
				if (rnd->next01() < probability) { // ����
					incoming_radiance = (radiance(0, env, reflection_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)) * Re;
					weight = RGB2Spectrum(hitpos_material.color, wavelength) / (probability * russian_roulette_probability) / (1.0 - scattering_probability);

					// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
					incoming_radiance = incoming_radiance + direct_light;

					//light intensity weight
					double lightIntensity = 1.0;
					if (hitpoint.bump)
					{
						lightIntensity = std::max(0.1, dot(reflection_ray.dir, hitpoint.bump_new_normal)) / dot(reflection_ray.dir, hitpoint.normal);
					}
					weight = weight*lightIntensity;

				}
				else { // ����
					incoming_radiance = (radiance(0, env, tdir, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)) * Tr;
					weight = RGB2Spectrum(hitpos_material.color, wavelength) / ((1.0 - probability) * russian_roulette_probability) / (1.0 - scattering_probability);

					// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
					incoming_radiance = incoming_radiance + direct_light_refraction;

					//light intensity weight
					double lightIntensity = 1.0;
					if (hitpoint.bump)
					{
						lightIntensity = std::max(0.1, dot(tdir.dir, hitpoint.bump_new_normal)) / dot(tdir.dir, hitpoint.normal);;;
					}
					weight = weight*lightIntensity;

				}
			}
			else { // ���܂Ɣ��˂̗�����ǐ�

				if (hitpoint.bump)
				{
					weight = RGB2Spectrum(hitpos_material.color, wavelength) / (russian_roulette_probability) / (1.0 - scattering_probability);

					Spectrum weight1 = weight;
					Spectrum weight2 = weight;

					double lightIntensity1 = 1.0;
					double lightIntensity2 = 1.0;
					if (hitpoint.bump)
					{
						lightIntensity1 = std::max(0.1, dot(reflection_ray.dir, hitpoint.bump_new_normal));
						lightIntensity2 = std::max(0.1, dot(tdir.dir, hitpoint.bump_new_normal));
					}
					weight1 = weight1*lightIntensity1;
					weight2 = weight2*lightIntensity2;

					incoming_radiance =
						(radiance(0, env, reflection_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)) * Re * weight1 +
						(radiance(0, env, tdir, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)) * Tr * weight2;

					// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
					incoming_radiance = incoming_radiance + direct_light*weight1 + direct_light_refraction*weight2;
					weight = ONE();
				}
				else
				{
					incoming_radiance =
						(radiance(0, env, reflection_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)) * Re +
						(radiance(0, env, tdir, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia)) * Tr;
					weight = RGB2Spectrum(hitpos_material.color, wavelength) / (russian_roulette_probability) / (1.0 - scattering_probability);

					// emission �͔��˗�����Z����Ă��Ȃ��̂Ō��weight���|����K�v�����邽��incoming_radiance�ɉ��Z���Ă���
					incoming_radiance = incoming_radiance + direct_light + direct_light_refraction;
				}
			}
		}

	Rtn:;
		//�֗^�}�����l������ꍇ�͌��̌������l������
		if (participatingMedia)
		{
			//�����i���̂̊O��)
			//���݂̌�_�ł̓��C�͕��̂̒��ɓ��낤�Ƃ�������Ȃ畨�̂̊O������
			if (inside_dir)
			{
				emission = transmittance_ratio * emission;
				incoming_radiance = transmittance_ratio * incoming_radiance;
			}
		}
		return RGB2Spectrum(env->Environment_Light, wavelength) + RGB2Spectrum(emission, wavelength) + MULTIPLY(weight, incoming_radiance);

	}

	// �� = Relative index of refraction
	// Fdr = Average Fresnel reflectance
	// Ft = Fresnel transmittance
	// Fr = Fresnel reflectance
	// Fdt = Average Fresnel transmittance

	/* Fdr Numerically approximate the diffuse Fresnel reflectance */
	inline double fresnelDiffuseReflectance(double eta)
	{
		/* Fast mode: the following code approximates the
		* diffuse Frensel reflectance for the eta<1 and
		* eta>1 cases. An evalution of the accuracy led
		* to the following scheme, which cherry-picks
		* fits from two papers where they are best.
		*/
		if (eta < 1) {
			/* Fit by Egan and Hilgeman (1973). Works
			reasonably well for "normal" IOR values (<2).

			Max rel. error in 1.0 - 1.5 : 0.1%
			Max rel. error in 1.5 - 2   : 0.6%
			Max rel. error in 2.0 - 5   : 9.5%
			*/
			return -1.4399f * (eta * eta)
				+ 0.7099f * eta
				+ 0.6681f
				+ 0.0636f / eta;
		}
		else {
            //return -0.4399 + 0.7099 / eta - 0.3319 / (eta * eta) + 0.0636 / (eta * eta * eta);

			/* Fit by d'Eon and Irving (2011)
			*
			* Maintains a good accuracy even for
			* unrealistic IOR values.
			*
			* Max rel. error in 1.0 - 2.0   : 0.1%
			* Max rel. error in 2.0 - 10.0  : 0.2%
			*/
			double invEta = 1.0f / eta,
				invEta2 = invEta*invEta,
				invEta3 = invEta2*invEta,
				invEta4 = invEta3*invEta,
				invEta5 = invEta4*invEta;

			return 0.919317f - 3.4793f * invEta
				+ 6.75335f * invEta2
				- 7.80989f * invEta3
				+ 4.98554f * invEta4
				- 1.36881f * invEta5;
		}
	}


};

#endif
