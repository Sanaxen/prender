#ifndef _DIRECTLIGHT_RADIANCE_H

#define _DIRECTLIGHT_RADIANCE_H

#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "circle.h"
#include "intersection.h"
#include "random.h"

#include "entity.h"
#include "scene_env.h"
#include "Spectrum.h"

#include "ParticipatingMedia.h"

namespace prender {
	inline Spectrum eval_direct_radiance(const SceneEnv* env, const Vector3d &ray_dir, Random *rnd, const Light& light, const Entity* now_object,
		const IntersectionPos &v0, const Vector3d refraction_normal, const Vector3d next_dir,
		const Vector3d light_pos,
		const int depth, const double wavelength)
	{
		// �����ʒu�̖@���i���̂���̃��C�̓��o���l���j
		const Vector3d orienting_normal = refraction_normal;

		Vector3d light_normal;
		double pdf;
		Color emission;

		if (!light.infinity_light)
		{
			// ������̈�_
			switch (light.light->type)
			{
			case ENTITY_TYPE_SPHERE:
			{
				Sphere* s_light = (Sphere*)light.light;
				light_normal = normalize(light_pos - s_light->position);
				if (s_light->normal_vector_inverse < 0)
				{
					light_normal = light_normal*-1.0;
				}
				pdf = 1.0 / s_light->area;

			}
			break;
			case ENTITY_TYPE_UVPLANE:
			{
				UVPlane* s_light = (UVPlane*)light.light;
				light_normal = s_light->normal;
				if (s_light->normal_vector_inverse < 0)
				{
					light_normal = light_normal*-1.0;
				}
				pdf = 1.0 / s_light->area;
			}
			break;
			case ENTITY_TYPE_CIRCLE:
			{
				Circle* s_light = (Circle*)light.light;
				light_normal = s_light->normal;
				if (s_light->normal_vector_inverse < 0)
				{
					light_normal = light_normal*-1.0;
				}
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
		double dot0 = dot(orienting_normal, light_dir);
		double dot1 = dot(light_normal, -1.0 * light_dir);

		Spectrum radiance = ZERO();

		if (light.infinity_light)
		{
			//���̕����͈��
			light_dir = light.dir*-1.0;
			pdf = 1.0;

			//���̑��̐ݒ�ύX
			dot1 = 1.0;
			dist2 = 1.0;//�����Ɉˑ����Ȃ�(���G���v�Z���鎞�ɏȗ����邽�߂�1.0�ɂ��Ă���)
			dot0 = dot(orienting_normal, light_dir);
			//fprintf(stderr, "infinity_light dir %f %f %f\n", light.dir.x, light.dir.y, light.dir.z);
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
			dot0 = dot(orienting_normal, light_dir);
			//fprintf(stderr, "parallel_light dir %f %f %f\n", light.dir.x, light.dir.y, light.dir.z);
		}

		//Spot�����̏ꍇ
		if (light.spot_light)
		{
			//���̕����͈��
			light_dir = normalize(((Sphere*)light.light)->position - v0.position);
			pdf = 1.0;

			//���̑��̐ݒ�ύX
			dot1 = 1.0;
			dist2 = 1.0;//�����Ɉˑ�����(���emittion���v�Z����Ƃ��ɐ��荞�܂��̂ł����ł͏ȗ����邽�߂�1.0�ɂ��Ă���)
		}

		if (dot0 >= 0 && dot1 >= 0)
		{
			//�W�I���g���^�[��
			const double G = dot0 * dot1 / dist2;

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
						//�T���v�������������������Ńq�b�g�����̂�IBL�������ꍇ
						//IBL��薳����������D�悷��
						if (light.infinity_light && intersection.hitpoint.material.IBL())
						{
							/* empty */
						}
						else
						{
							return ZERO();
						}
					}
					//printf("%p %p\n", (void*)EntList->List[intersection.object_id], (void*)light);

					if (light.light && light.light->material()->IBL())
					{
						emission = intersection.hitpoint.material.emission;
					}
					if (light.spot_light)
					{
						emission = light.SpotLightEmission(emission, ((Sphere*)light.light)->position, v0.position);
					}

					if (light.infinity_light)
					{
						emission = light.infinity_light_emitssin;
					}

					if (env->participatingMedia /*&& !light.parallel_light*/ && (light.light && !light.light->material()->IBL()))
					{
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

						//�����̌���
						emission = emission * transmittance_ratio;
					}

					//const IntersectionPos &hitpoint = intersection.hitpoint;

					//���C�g����Փ˓_�ւ̃x�N�g��
					Vector3d orienting_light_dir = (-1.0)*light_dir;

					switch (v0.material.reflection_type)
					{
					case REFLECTION_TYPE_SSS_WARD_BRDF:
					case REFLECTION_TYPE_WARD_BRDF:
					{
						const double alpha_x = v0.material.brdfParameter.ward_brdf.alp_x;
						const double alpha_y = v0.material.brdfParameter.ward_brdf.alp_y;

						Vector3d in;
						Vector3d halfv;
						Vector3d dir;

						// orienting_normal�̕�������Ƃ������K�������(w, u, v)�����B
						Vector3d w, u, v;
						OrthonormalBasis(orienting_normal, w, u, v);

						in = (-1.0)*ray_dir;
						halfv = normalize(in + light_dir);

						dir = light_dir;

						double ww0 = pow(dot(halfv, orienting_normal), 2.0);
						if (fabs(ww0) < 1.0e-16) return ZERO();
						double ww1 = (pow(dot(halfv, u) / alpha_x, 2.0) + pow(dot(halfv, v) / alpha_y, 2.0)) / ww0;
						double ww2 = 4.0*PS_PI*alpha_x*alpha_y*sqrt(dot(dir, orienting_normal)* dot(in, orienting_normal));
						
						if ( fabs(ww2) < 1.0e-16 ) return ZERO();
						const double ww = exp(-ww1) / ww2;
						if ( ww < 0.0 ) return ZERO();

						Spectrum weight = RGB2Spectrum(v0.material.specular, wavelength) * ww;

						radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), weight)* G / pdf;
					}
					break;
					case REFLECTION_TYPE_PHONG_BRDF:
					{
						//brdf = ((n + 2) / 2��) cos(��)^n
						Vector3d ref = normalize(ray_dir - orienting_normal*2.0*dot(orienting_normal, ray_dir));
						if (dot(ref, light_dir) < 0.0) return ZERO();
						const double ww1 = pow(dot(ref, light_dir), v0.material.brdfParameter.phong_brdf.specular_exponent);

						const double ww = ((v0.material.brdfParameter.phong_brdf.specular_exponent + 2) / (PS_TWOPI))* ww1;
						if ( ww < 0.0 ) return ZERO();
						radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.specular*ww, wavelength))* G / pdf;
						break;
					}
#if 0	//�v��Ȃ������I�I
					case REFLECTION_TYPE_SPECULAR:
					{			
						radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.color, wavelength)* (1.0 / PS_PI))* G / pdf;
						break;
					}
					break;
					case REFLECTION_TYPE_SSS_REFRACTION:
					case REFLECTION_TYPE_REFRACTION:
					{
						radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.color, wavelength)* (1.0 / PS_PI))* G / pdf;
						break;
					}
					break;
#endif
					case REFLECTION_TYPE_DIFFUSE:
					case REFLECTION_TYPE_SSS_REFRACTION:

						if (light.infinity_light)
						{
							//����diffuse�ł͖����̂�BRDF=color��(1/��)�̌W���͂�����Ȃ�
							radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.color, wavelength))* G / pdf;
							//fprintf(stderr, "radiance %f %f %f\n", radiance.x, radiance.y, radiance.z);
						}else
						//���s�����̏ꍇ
						if (light.parallel_light)
						{
							//����diffuse�ł͖����̂�BRDF=color��(1/��)�̌W���͂�����Ȃ�
							radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.color, wavelength))* G / pdf;
						}
						else
						if (light.spot_light)
						{
							//����diffuse�ł͖����̂�BRDF=color��(1/��)�̌W���͂�����Ȃ�
							radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.color, wavelength))* G / pdf;
						}
						else
						{
							//���ڌ��̑Ώۂ͊��S�g�U�Ƃ��Ă���̂� BRDF= color/��
							radiance = MULTIPLY(RGB2Spectrum(emission, wavelength), RGB2Spectrum(v0.material.color, wavelength)* (PS_INV_PI))* G / pdf;
						}
					}
				}
			}
		}

		//light intensity weight
		double lightIntensity = 1.0;
		if (v0.bump)
		{
			lightIntensity = std::max(0.0, dot(light_dir, v0.bump_new_normal)) / dot(light_dir, v0.normal);
		}

		radiance = lightIntensity * radiance;
		return radiance;
	}



	// ������̓_���T���v�����O���Ē��ڌ����v�Z����
	inline Spectrum direct_radiance_sample(const SceneEnv* env, const Vector3d &ray_dir, Random *rnd, const Light& light, const Entity* now_object,
		const IntersectionPos &v0, const Vector3d refraction_normal, const Vector3d next_dir, const int depth, const double wavelength)
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
				//const double t = rnd->Next01();
				//const double s = 4.0*PS_PI * rnd->Next01();
				//const double tt = 1 - t;
				//const double tt2 = sqrt(1.0 - tt*tt);

				//const double r = s_light->radius;
				//light_pos = normalize(Vector3d(cos(s) * tt2, sin(s) * tt2, tt))*(r+PS_EPS) + s_light->position;
				//pdf = 1.0/(4.0 * PS_PI * r*r);

				const double r1 = PS_TWOPI * rnd->Next01();
				const double r2 = 1.0 - 2.0 * rnd->Next01();
				const double sqrtr22 = sqrt(1.0 - r2*r2);
				light_pos = s_light->position + ((s_light->radius + PS_EPS) * Vector3d(sqrtr22 * cos(r1), sqrtr22 * sin(r1), r2));

				if (s_light->hemisphere)
				{
					const double r1 = PS_TWOPI * rnd->Next01();
					const double r2 = rnd->Next01();
					light_pos = s_light->position + ((s_light->radius + PS_EPS) * Vector3d(sqrtr22 * cos(r1), sqrtr22 * sin(r1), r2));
				}
			}
			break;
			case ENTITY_TYPE_UVPLANE:
			{
				UVPlane* s_light = (UVPlane*)light.light;
				const double t = s_light->v_length * rnd->Next01();
				const double s = s_light->u_length * rnd->Next01();
				light_normal = s_light->normal;

				light_pos = s_light->org + s_light->u_axis*s + s_light->v_axis*t + light_normal*PS_EPS;
			}
			break;
			case ENTITY_TYPE_CIRCLE:
			{
				Circle* s_light = (Circle*)light.light;

				const double r1 = rnd->Next01();
				const double r2 = rnd->Next01();
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

		return eval_direct_radiance(env, ray_dir, rnd, light, now_object,
			v0, refraction_normal, next_dir, light_pos, depth, wavelength);
	}

	//����������
	inline Spectrum infinity_light_direct(const SceneEnv* env, const Intersection& isect, const Vector3d& refdir, int &exist_infinity_light, const double wavelength)
	{
		Spectrum direct_light = ZERO();

		//��_�������i�������j����������Hit�@�܂��́@IBL��Hit�������������ɂ�Hit�����Ȃ疳������D��
		if ((isect.object_id >= 0 && EntList->List[isect.object_id]->material()->IBL()) || isect.object_id < 0)
		{
			direct_light = const_cast<SceneEnv*>(env)->light_list.InfinityLight(refdir, exist_infinity_light, wavelength);
		}
		return direct_light;
	}

	//���s����
	inline Spectrum parallel_light_direct(const IntersectionPos &v0, const Intersection& isect, const Vector3d& refdir, const Light& light, const double wavelength)
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

	//Spot����
	inline Spectrum spot_light_direct(const IntersectionPos &v0, const Intersection& isect, const Vector3d& refdir, const Light& light, const double wavelength)
	{
		Spectrum direct_light = ZERO();

		if (light.spot_light)
		{
			Color& Intensity = light.SpotLightEmission(isect.hitpoint.material.emission, ((Sphere*)light.light)->position, v0.position);
			direct_light = RGB2Spectrum(Intensity, wavelength);

			const Vector3d light_dir = normalize(((Sphere*)light.light)->position - v0.position);

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

#if 0
	Spectrum direct_light_driver(const int direct_light_meshod, const Entity* now_object, const Vector3d& dir, const IntersectionPos& hitpoint, const Vector3d& orienting_normal, const SceneEnv* env, const Ray &ray, Random *rnd, const int depth, const double wavelength, const int nextEventEstimation = 0, const int participatingMedia = 0)
	{
		Spectrum direct_light = ZERO();
		if (direct_light_meshod && nextEventEstimation && env->light_list.list.size())
		{
			// ���ڌ��̕]��
			//�����ƌ������T���v�����O���Ă��̃T���v�����O�ʒu���猻�݈ʒu�����Ԍo�H�Ԃŕ��ˋP�x�����߂�
			double light_probability = 1.0;
			for (int i = 0; i < SHADOW_RAY_SAMPLING; i++)
			{
				const Light& light = const_cast<SceneEnv*>(env)->light_list.randomLight(rnd, light_probability);
				direct_light = direct_light + direct_radiance_sample(env, ray.dir, rnd, light, now_object, hitpoint, orienting_normal, dir, depth, wavelength) / light_probability;

			}
			direct_light = direct_light / SHADOW_RAY_SAMPLING;
			//if ( shadow_ray_hit_count ) fprintf(stderr, "shadow_ray_hit_count %d\n", shadow_ray_hit_count);

			//�֗^�}���ɂ�錸�����l��
			if (participatingMedia)
			{
				direct_light = env->participatingMediaParam.transmittanceRatio(hitpoint.distance,wavelength) * direct_light;
			}
		}
		//direct_light �͂��łɔ��˗�����Z�ς݂Ȃ̂ŁAweight���|����K�v�͂Ȃ�
		//�Ȃ̂Ō��weight*incoming_radiance�����邪incoming_radiance�ɂ͉��Z���Ȃ�
		//�܂�@direct_light + weight*incoming_radiance

		if (!direct_light_meshod && nextEventEstimation)
		{
			// ���ڌ��̕]��
			Intersection lid;
			if (intersect_scene(Ray(hitpoint.position, dir), &lid))
			{
				const Entity *now_object2 = EntList->List[lid.object_id];
				if (now_object2->light)
				{
					direct_light = RGB2Spectrum(env->light_list.list[now_object2->light_id].light->material()->emission, wavelength);

					//�֗^�}���ɂ�錸�����l��
					if (participatingMedia)
					{
						direct_light = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance,wavelength) * direct_light;
					}
				}
			}
		}
		return direct_light;
	}
#endif

	inline Spectrum direct_light_perfect_refrection(const SceneEnv* env, const IntersectionPos &hitpoint, const Ray& ray, const int depth, const double wavelength, const int participatingMedia = 0)
	{
		Spectrum direct_light_refraction = ZERO();
		Intersection lid;
		Entity *now_object2 = 0;

		//��_�����邩
		bool isect = intersect_scene(ray, &lid, depth);

		if (isect) now_object2 = EntList->List[lid.object_id];

		int exist_infinity_light = 0;

		//��_�������i�������j����������Hit�@�܂��́@IBL��Hit�������������ɂ�Hit�����Ȃ疳������D��
		direct_light_refraction = infinity_light_direct(env, lid, ray.dir, exist_infinity_light, wavelength);

		//�����������͂Ȃ������Ƀq�b�g
		if (now_object2 && now_object2->light && !exist_infinity_light)
		{

			direct_light_refraction = RGB2Spectrum(lid.hitpoint.material.emission, wavelength);

			//���s��������̊�^�����邩
			const Light& light = env->light_list.list[now_object2->light_id];
			if (light.parallel_light)
			{
				direct_light_refraction = parallel_light_direct(hitpoint, lid, ray.dir, light, wavelength);
			}
			if (light.spot_light)
			{
				direct_light_refraction = spot_light_direct(hitpoint, lid, ray.dir, light, wavelength);
			}
			//�֗^�}���ɂ�錸�����l��
			if (participatingMedia)
			{
				direct_light_refraction = env->participatingMediaParam.transmittanceRatio(lid.hitpoint.distance, wavelength) * direct_light_refraction;
			}
		}

		return direct_light_refraction;
	}

};

#endif
