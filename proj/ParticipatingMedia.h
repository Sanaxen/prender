#ifndef _PARTICIPATING_MEDIA_H

#define _PARTICIPATING_MEDIA_H

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

	Spectrum radiance(int* direct_light_hit, const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const double wavelength, const int nextEventEstimation, const int participatingMedia);

	//const int kDepth = 4; // ���V�A�����[���b�g�őł��؂�Ȃ��ő�[�x
	//const int kDpethLimit = 64;



	inline Spectrum eval_direct_radiance_media(const SceneEnv* env, Random *rnd, const Light& light, const IntersectionPos &v0, const Vector3d light_pos, const int depth, const double wavelength)
	{
		Vector3d light_normal;
		double pdf;
		Color emission;

		if (!light.infinity_light)
		{
			// ������̈�_���T���v�����O����
			switch (light.light->type)
			{
			case ENTITY_TYPE_SPHERE:
			{
				Sphere* s_light = (Sphere*)light.light;
				light_normal = normalize(light_pos - s_light->position);
				pdf = 1.0 / s_light->area;

			}
			break;
			case ENTITY_TYPE_UVPLANE:
			{
				UVPlane* s_light = (UVPlane*)light.light;
				light_normal = s_light->normal;
				pdf = 1.0 / s_light->area;
			}
			break;
			case ENTITY_TYPE_CIRCLE:
			{
				Circle* s_light = (Circle*)light.light;
				light_normal = s_light->normal;
				pdf = 1.0 / s_light->area;
			}
			break;
			};
			emission = light.light->material()->emission;
		}
		else
		{
			light_normal = light.dir;
			pdf = 1.0;
			emission = light.infinity_light_emitssin;
		}


		//������̂P�_���烌�C�̕��������߂�B
		Vector3d light_dir = normalize(light_pos - v0.position);
		double dist2 = (light_pos - v0.position).sqr();
		double dot1 = dot(light_normal, -1.0 * light_dir);

		if (light.infinity_light)
		{
			//���̕����͈��
			light_dir = light.dir*-1.0;
			pdf = 1.0;

			//���̑��̐ݒ�ύX
			dot1 = 1.0;
			dist2 = 1.0;//�����Ɉˑ����Ȃ�(���G���v�Z���鎞�ɏȗ����邽�߂�1.0�ɂ��Ă���)
		}
		//���s�����̏ꍇ
		if (light.parallel_light)
		{
			//���̕����͈��
			light_dir = light.dir*-1.0;
			pdf = 1.0;

			//���̑��̐ݒ�ύX
			dot1 = 1.0;
			dist2 = 1.0;//�����Ɉˑ����Ȃ�(���G���v�Z���鎞�ɏȗ����邽�߂�1.0�ɂ��Ă���)
		}

		//�X�|�b�g����
		if (light.spot_light)
		{
			//���̕����͈��
			light_dir = normalize(((Sphere*)light.light)->position - v0.position);
			pdf = 1.0;

			//���̑��̐ݒ�ύX
			dot1 = 1.0;
			dist2 = 1.0;//�����Ɉˑ����Ȃ�(���G���v�Z���鎞�ɏȗ����邽�߂�1.0�ɂ��Ă���)
		}

		if (dot1 >= 0)
		{
			//�W�I���g���^�[��
			const double G = dot1 / dist2;

			//�Փ˓_���烉�C�g�����ւ̃��C
			Ray light_ray = Ray(v0.position, light_dir);

			Intersection intersection;
			bool isect = intersect_scene(light_ray, &intersection, depth);
			if (isect || light.infinity_light)
			{
				bool cond = false;
				cond = fabs(sqrt(dist2) - intersection.hitpoint.distance) < 1e-3;

				if (light.infinity_light)
				{
					cond = true;
				}
				if (light.parallel_light)
				{
					cond = true;
				}
				if (light.spot_light)
				{
					cond = true;
				}

				if (cond)
				{
					//��Q���Ȃ璼�ڌ��͖���
					if (intersection.object_id >= 0 && EntList->List[intersection.object_id] != light.light)
					{
						if (light.infinity_light && intersection.hitpoint.material.IBL())
						{
							/* empty */
						}
						else
						{
							return ZERO();
						}
					}
					//printf("%p %p\n", (void*)EntList->List[intersection.object_id], (void*)light.light);

					//���ڌ��̑Ώۂ͊��S�g�U�Ƃ��Ă���̂� BRD F= color/��(��Ԓ��̓_�������󂯂�)
					if (light.light && light.light->material()->IBL())
					{
						//return RGB2Spectrum(intersection.hitpoint.material.emission, wavelength)* (1.0 / PS_PI)* G / pdf;
						return RGB2Spectrum(emission, wavelength)* (1.0 / PS_PI)* G / pdf;
					}
					if (light.infinity_light)
					{
						//����diffuse�ł͖����̂�BRDF=color��(1/��)�̌W���͂�����Ȃ�
						return RGB2Spectrum(emission, wavelength)* G / pdf;
					}
					else
						//���s�����̏ꍇ
						if (light.parallel_light)
						{
							//����diffuse�ł͖����̂�BRDF=color��(1/��)�̌W���͂�����Ȃ�
							return RGB2Spectrum(emission, wavelength)* G / pdf;
						}
						else
							//Spot�����̏ꍇ
							if (light.spot_light)
							{
								emission = light.SpotLightEmission(emission, ((Sphere*)light.light)->position, v0.position);

								//����diffuse�ł͖����̂�BRDF=color��(1/��)�̌W���͂�����Ȃ�
								return RGB2Spectrum(emission, wavelength)* G / pdf;
							}
							else
							{
								//���ڌ��̑Ώۂ͊��S�g�U�Ƃ��Ă���̂� BRDF= color/��
								return RGB2Spectrum(emission, wavelength)* (1.0 / PS_PI)* G / pdf;
							}
				}
			}
		}

		return ZERO();
	}

	// ������̓_���T���v�����O���Ē��ڌ����v�Z����
	inline Spectrum direct_radiance_sample_media(const SceneEnv* env, Random *rnd, const Light& light, const IntersectionPos &v0, const int depth, const double wavelength)
	{
		Vector3d light_pos;
		Vector3d light_normal;

		double pdf;

		if (!light.infinity_light)
		{
			// ������̈�_���T���v�����O����
			switch (light.light->type)
			{
			case ENTITY_TYPE_SPHERE:
			{
				Sphere* s_light = (Sphere*)light.light;
				//const double t = rnd->next01();
				//const double s = 4.0*PS_PI * rnd->next01();
				//const double tt = 1 - t;
				//const double tt2 = sqrt(1.0 - tt*tt);

				//const double r = s_light->radius;
				//light_pos = normalize(Vector3d(cos(s) * tt2, sin(s) * tt2, tt))*(r+PS_EPS) + s_light->position;
				//light_normal = normalize( light_pos - s_light->position);
				//pdf = 1.0/(4.0 * PS_PI * r*r);

				light_normal = normalize(light_pos - s_light->position);
				const double r1 = 2 * PS_PI * rnd->next01();
				const double r2 = 1.0 - 2.0 * rnd->next01();
				light_pos = s_light->position + ((s_light->radius + PS_EPS) * Vector3d(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2));
			}
			break;
			case ENTITY_TYPE_UVPLANE:
			{
				UVPlane* s_light = (UVPlane*)light.light;
				light_normal = s_light->normal;
				const double t = s_light->v_length * rnd->next01();
				const double s = s_light->u_length * rnd->next01();

				light_pos = s_light->org + s_light->u_axis*s + s_light->v_axis*t;
			}
			break;
			case ENTITY_TYPE_CIRCLE:
			{
				Circle* s_light = (Circle*)light.light;

				const double r1 = rnd->next01();
				const double r2 = rnd->next01();
				const double r = sqrt(r1);
				const double th = r2*2.0*PS_PI;
				light_pos = s_light->org + s_light->u_axis*(r*cos(th)) + s_light->v_axis*(r*sin(th));
			}
			break;
			};
		}
		else
		{
			light_normal = light.dir;
			pdf = 1.0;
		}

		double t = (v0.position - light_pos).length();
		if (light.spot_light)
		{
			t = (v0.position - ((Sphere*)light.light)->position).length();
		}
		if (light.infinity_light)
		{
			t = 0.0;
		}
		const Spectrum transmittance_ratio = env->participatingMediaParam.transmittanceRatio(t, wavelength);

		Spectrum radiance = eval_direct_radiance_media(env, rnd, light, v0, light_pos, depth, wavelength);

		return radiance * transmittance_ratio;
	}

	//���s����
	inline Spectrum parallel_light_direct_media(const IntersectionPos &v0, const Intersection& isect, const Vector3d& refdir, const Light& light, const double wavelength)
	{
		Spectrum direct_light = ZERO();

		if (light.parallel_light)
		{
			direct_light = RGB2Spectrum(isect.hitpoint.material.emission, wavelength);
			const Vector3d light_dir = light.dir*-1.0;
			const double dot0 = dot(refdir, light_dir);
			const double dot1 = 1.0;	//���s�����Ƃ��Ă���̂�

			const double dist2 = (isect.hitpoint.position - v0.position).sqr();
			if (dot0 >= 0.0 && dot1 >= 0.0)
			{
				//�W�I���g���^�[��
				//const double G = dot0 * dot1 / dist2;
				const double G = dot0 * dot1;	//�����ɖ��֌W

				//�Ƃ炳��Ă���
				direct_light = direct_light*G;
			}
			else
			{
				direct_light = ZERO();
			}
		}
		return direct_light;
	}

	inline Spectrum participatingMedia_radiance(const SceneEnv* env, Random* rnd, const Ray& ray, const IntersectionPos &hitpoint, const Entity* now_object, const Vector3d& orienting_normal, const double russian_roulette_probability, const double scattering_probability, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		const Material& hitpos_material = hitpoint.material;
		Color emission = ZERO();					//���Ȕ���
		Spectrum direct_light = ZERO();				//���ڌ�(����)
		Spectrum direct_light_refraction = ZERO();	//���ڌ�(����)
		Spectrum incoming_radiance = 0.0;
		Spectrum weight = 1.0;

		emission = hitpos_material.emission;		// Init color to emissive value

		//Spectrum transmittance_ratio = ONE();
		////�֗^�}�����l������ꍇ�͌��̌������l������
		//if (participatingMedia) transmittance_ratio = env->participatingMediaParam.transmittanceRatio(hitpoint.distance, wavelength);

		{
			Vector3d next_dir;
			Ray next_ray = ray;
			double d;
			double pdf;
			Spectrum transmittance_ratio;
			double dirPdf = PS_INV_FORPI;
			double phase = PS_INV_FORPI; // �����U��
			const int select = std::min(2, (int)(3.0 * rnd->next01()));

			int stat = -1;

			for (int k = 0; k < 1; k++)
			{
#if 10
				// ���͂���̉e�����l���i�܂�܂�肩�痈�������U������ray�����ɗ����j
				// �����ɉ������d�_�T���v�����O
#if 10
				const double u = rnd->next01();
				d = -log(1.0 + u * (exp(-env->participatingMediaParam.transmittance * hitpoint.distance) - 1.0)) / env->participatingMediaParam.transmittance;
				pdf = exp(-env->participatingMediaParam.transmittance * d) * (-env->participatingMediaParam.transmittance / (exp(-env->participatingMediaParam.transmittance * hitpoint.distance) - 1.0));
				//if (d > hitpoint.distance) fprintf(stderr, "error.\n");
#else		
				const double eps = 1.0e-10;
				double d = 0.0;
				do {
					d = -log(rnd->next01()) / env->participatingMediaParam.transmittance;
				} while (d >= hitpoint.distance*(1.0 - eps));
				const double pdf = exp(-env->participatingMediaParam.transmittance * d);
#endif

				// �i�n�_���狗��d�̈ʒu����������ς��)
				if (env->participatingMediaParam.phase_prm[select] == 0.0 || fabs(env->participatingMediaParam.phase_prm[select]) > 1.0)
				{
					//�����U��
					//const double r1 = 2 * PS_PI * rnd->next01();
					//const double r2 = 1.0 - 2.0 * rnd->next01();
					//next_dir = normalize(Vector3d(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2));
					sampleSphere(rnd, ray.dir, next_dir);
					phase = PhaseFunction_Isotropic();
					dirPdf = phase;
				}
				else
				{
					sampleHG(rnd, ray.dir, env->participatingMediaParam.phase_prm[select], next_dir);
					phase = PhaseFunction_HenyeyGreenstein(env, env->participatingMediaParam.phase_prm[select], dot(next_dir, ray.dir));
					dirPdf = phase;
				}

				//���C���_�ƏՓ˓_�̓r���ʒu����next_dir�ւ̃��C
				next_ray = Ray(ray.org + d * ray.dir, next_dir);
				transmittance_ratio = env->participatingMediaParam.transmittanceRatio(d, wavelength);

				if (1)
				{
					//Intersection intersectionWrk;
					//if (!intersect_scene(next_ray, &intersectionWrk))
					//{
					//	stat = 0;
					//	break;
					//}
					//Entity* object = env->EntList.List[intersectionWrk.object_id];
					//bool insidedir = object->isInsideDir(next_ray, intersectionWrk.hitpoint);

					//if (!insidedir)
					//{
					//	continue;
					//}
					//else
					{
						stat = 0;
						break;
					}
				}
			}
			if (stat != 0)
			{
				return ZERO();
			}

			IntersectionPos hitpoint_media = hitpoint;
			hitpoint_media.distance = d;
			hitpoint_media.position = next_ray.org;

			if (!now_object->light)
			{
				double light_probability = 1.0;
				if (env->light_list.list.size())
				{
					for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
					{
						const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
						direct_light = direct_light + direct_radiance_sample_media(env, rnd, light, hitpoint_media, depth, wavelength) / light_probability;
					}
					direct_light = direct_light / SHADOW_RAY_SAMPLING;
				}
			}


			if (pdf == 0.0) {
				return ZERO();
			}
			else {


				Spectrum incoming_radiance = ZERO();

				if (now_object->light)
				{
					if (depth == 0)
					{
						incoming_radiance = RGB2Spectrum(hitpoint.material.emission, wavelength);
						return transmittance_ratio*incoming_radiance;
					}
					else
					{
						return ZERO();
					}
				}
				Spectrum L = direct_light + radiance(0, env, next_ray, rnd, depth + 1, wavelength, nextEventEstimation, participatingMedia);

				L = transmittance_ratio * RGB2Spectrum(env->participatingMediaParam.scatteringColor, wavelength) * L
					* phase / pdf
					/ dirPdf
					/ scattering_probability
					/ russian_roulette_probability;

				//const Spectrum Tr = env->participatingMediaParam.transmittanceRatio(hitpoint.distance,wavelength);
				//L = L + Tr * radiance(0, env, next_ray, rnd, depth+1, wavelength, nextEventEstimation, participatingMedia)
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

	Rtn:;
		return ZERO();

	}

};

#endif
