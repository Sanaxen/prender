#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include <cmath>
#include "entity.h"

#include "ObjLoader.h"
#include <vector>

extern "C" int intersect_triangle(double orig[3], double dir[3],
	double vert0[3], double vert1[3], double vert2[3],
	double *t, double *u, double *v);
extern "C" int intersect_triangle3(const double orig[3], const double dir[3],
	const double vert0[3], const double vert1[3], const double vert2[3],
	double *t, double *u, double *v);

#pragma warning(disable:4244)
#pragma warning(disable:4996)

namespace prender {
	class Triangle : public Entity
	{

	public:
		::Obj* mesh;					//元のメッシュデータ
		int face_id;					//FaceのID
		int normal_id;				//Faceの法線ベクトル
		int vert_id[3];						//三角形の頂点座標ID
		int vert_normal_id[3];				//三角形の頂点法線ID
		int texture_idx[3];				//テクスチャ座標ID
		int material_texture_idx;		//テクスチャID

		inline Triangle()
		{
			material_texture_idx = -1;
		}

#if 0	//このアルゴリズムはNG!!（結果が思わしくない)
		int CheckIntersection(const Ray &ray, double& u, double& v, double& t) const
		{
			Vector3d edge1(vert[1] - vert[0]);
			Vector3d edge2(vert[2] - vert[0]);

			Vector3d P = cross(ray.dir, edge2);
			double det = dot(P, edge1);

			if (det > PS_EPS)
			{
				// solve u
				Vector3d T(ray.org - (vert[0]));
				u = dot(P, T);

				if (u > 0 && u < det)
				{
					// solve v
					Vector3d Q = cross(T, edge1);
					v = dot(Q, ray.dir);

					if (v>0 && u + v < det)
					{
						t = dot(Q, edge2) / det;

						if (t > PS_EPS)
						{
							u /= det;
							v /= det;
							return 1;
						}
					}
				}
			}
			return 0;
		}
#endif

		//Moller-Trumbore intersection alogorithm
		inline int CheckIntersection2(const Ray& ray, double& u, double& v, double& w, double& t) const
		{
			Vector3d  tVector3d, pVector3d, qVector3d;
			double det, inv_det;

			/* find Vector3dtors for two edges sharing vert0 */
			const Vector3d edge1 = mesh->vertex(vert_id[1]) - mesh->vertex(vert_id[0]);
			const Vector3d edge2 = mesh->vertex(vert_id[2]) - mesh->vertex(vert_id[0]);

			/* begin calculating determinant - also used to calculate U parameter */
			pVector3d = cross(ray.dir, edge2);

			/* if determinant is near zero, ray lies in plane of triangle */
			det = dot(edge1, pVector3d);

			/* calculate distance from vert0 to ray origin */
			tVector3d = ray.org - mesh->vertex(vert_id[0]);

			qVector3d = cross(tVector3d, edge1);

			if (det > PS_EPS)
			{
				inv_det = 1.0 / det;
				u = dot(tVector3d, pVector3d);
				if (u < 0.0 || u > det)	return 0;

				/* calculate V parameter and test bounds */
				v = dot(ray.dir, qVector3d);
				if (v < 0.0 || u + v > det)	return 0;

			}
			else if (det < -PS_EPS)
			{
				inv_det = 1.0 / det;
				/* calculate U parameter and test bounds */
				u = dot(tVector3d, pVector3d);
				if (u > 0.0 || u < det)	return 0;

				/* calculate V parameter and test bounds */
				v = dot(ray.dir, qVector3d);
				if (v > 0.0 || u + v < det)	return 0;
			}
			else return 0;  /* ray is parallell to the plane of the triangle */

			t = dot(edge2, qVector3d) * inv_det;
			u *= inv_det;
			v *= inv_det;
			w = 1.0 - u - v;

			return 1;
		}

		inline bool intersect(const Ray &ray, IntersectionPos *hitpoint) const
		{
			//printf("\n\npolygon-Ray-Triangle Intersection\n");
			bool status = false;

			IntersectionPos hit;
			double t, u, v, w;


			//Moller-Trumbore intersection alogorithm
#if 10
			//const double org[3] = {ray.org.x, ray.org.y, ray.org.z}; 
			//const double dir[3] = {ray.dir.x, ray.dir.y, ray.dir.z};

			//double vv[3][3];
			//
			//vv[0][0] = vert[0].x;
			//vv[0][1] = vert[0].y;
			//vv[0][2] = vert[0].z;

			//vv[1][0] = vert[1].x;
			//vv[1][1] = vert[1].y;
			//vv[1][2] = vert[1].z;

			//vv[2][0] = vert[2].x;
			//vv[2][1] = vert[2].y;
			//vv[2][2] = vert[2].z;
			//int stat = intersect_triangle3(org, dir, vv[0], vv[1], vv[2], &t, &u, &v);

			int stat = CheckIntersection2(ray, u, v, w, t);
#else
			int stat = CheckIntersection(ray, u, v, t);
#endif
			if (stat == 1 && t > PS_EPS)
			{
				//printf("%d %f\n", i, t);
				status = true;
				hit.distance = t;
				hit.position = ray.org + ray.dir*t;
				hit.normal = mesh->fnormal(normal_id);
				hit.material = *((Triangle*)this)->material();

				if (mesh->smooth)
				{
					const Vector3d vt[3] = {
						mesh->normal(vert_normal_id[0]),
						mesh->normal(vert_normal_id[1]),
						mesh->normal(vert_normal_id[2])
					};

					//const Vector3d normal2 = normalize(vt[0] + (vt[1]-vt[0])*u + (vt[2]-vt[0])*v);
					const Vector3d normal2 = normalize(vt[0] * w + vt[1] * u + vt[2] * v);	//上の式と同じ

					hit.normal = normal2;
				}
				if (normal_vector_inverse < 0)
				{
					hit.normal = hit.normal*-1.0;
				}

				if (mesh->uv.size())
				{
					//hit.u = mesh->uv[vtx[0]][0]*w + mesh->uv[vtx[1]][0]*u + mesh->uv[vtx[2]][0]*v; 
					//hit.v = mesh->uv[vtx[0]][1]*w + mesh->uv[vtx[1]][1]*u + mesh->uv[vtx[2]][1]*v;
					hit.u = mesh->uv[texture_idx[0]].g[0] * w + mesh->uv[texture_idx[1]].g[0] * u + mesh->uv[texture_idx[2]].g[0] * v;
					hit.v = mesh->uv[texture_idx[0]].g[1] * w + mesh->uv[texture_idx[1]].g[1] * u + mesh->uv[texture_idx[2]].g[1] * v;
				}

				if (mesh->set_color)
				{
					Color c[3] = {
						mesh->vertex_color(vert_id[0]),
						mesh->vertex_color(vert_id[1]),
						mesh->vertex_color(vert_id[2])
					};

					hit.material.color.x = c[0].x*w + c[1].x*u + c[2].x*v;
					hit.material.color.y = c[0].y*w + c[1].y*u + c[2].y*v;
					hit.material.color.z = c[0].z*w + c[1].z*u + c[2].z*v;
				}
				if (mesh->texture_material.size() && material_texture_idx >= 0)
				{
					if (mesh->texture_material[material_texture_idx].textureImage)
					{
						BitMap* image = (BitMap*)mesh->texture_material[material_texture_idx].textureImage;

						double dmy;
						double uu = hit.u;
						double vv = hit.v;

						if (uu < 0 || uu > 1) uu = modf(uu, &dmy);
						if (vv < 0 || vv > 1) vv = modf(vv, &dmy);
						if (uu < 0) uu = 1.0 - fabs(uu);
						if (vv < 0) vv = 1.0 - fabs(vv);

						int j = (image->W() - 1)*uu;
						int i = (image->H() - 1)*vv;

						//const Color& diffuse = Color(mesh->texture_material[material_texture_idx].diffuse);
						//hit.material.color.x = diffuse.x*(double)image->cell(i,j).r/255.0;
						//hit.material.color.y = diffuse.y*(double)image->cell(i,j).g/255.0;
						//hit.material.color.z = diffuse.z*(double)image->cell(i,j).b/255.0;

						const Color& diffuse = ((Triangle*)this)->material()->color;
						hit.material.color.x = diffuse.x*(double)image->cell(i, j).r / 255.0;
						hit.material.color.y = diffuse.y*(double)image->cell(i, j).g / 255.0;
						hit.material.color.z = diffuse.z*(double)image->cell(i, j).b / 255.0;

						if (mesh->texture_material[material_texture_idx].user_ReflectionType >= 0)
						{
							//hit.material.color = Color(1,0,0);	//<== Debug用
							hit.material.reflection_type = (ReflectionType)mesh->texture_material[material_texture_idx].user_ReflectionType;
							hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0] * (double)image->cell(i, j).r / 255.0;;
							hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1] * (double)image->cell(i, j).g / 255.0;;
							hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2] * (double)image->cell(i, j).b / 255.0;;

							hit.material.roughness = mesh->texture_material[material_texture_idx].roughness;

							if (hit.material.reflection_type == 3)
							{
								hit.material.ward_brdf.alp_x = mesh->texture_material[material_texture_idx].ward[0];
								hit.material.ward_brdf.alp_y = mesh->texture_material[material_texture_idx].ward[1];
								hit.material.ward_brdf.specular.x = mesh->texture_material[material_texture_idx].specular.g[0] * (double)image->cell(i, j).r / 255.0;;
								hit.material.ward_brdf.specular.y = mesh->texture_material[material_texture_idx].specular.g[1] * (double)image->cell(i, j).g / 255.0;;
								hit.material.ward_brdf.specular.z = mesh->texture_material[material_texture_idx].specular.g[2] * (double)image->cell(i, j).b / 255.0;;

								hit.material.ward_brdf.diffuse.x = mesh->texture_material[material_texture_idx].diffuse.g[0] * (double)image->cell(i, j).r / 255.0;;
								hit.material.ward_brdf.diffuse.y = mesh->texture_material[material_texture_idx].diffuse.g[1] * (double)image->cell(i, j).g / 255.0;;
								hit.material.ward_brdf.diffuse.z = mesh->texture_material[material_texture_idx].diffuse.g[2] * (double)image->cell(i, j).b / 255.0;;
							}
						}
					}
					else
					{
						const Color& diffuse = Color(mesh->texture_material[material_texture_idx].diffuse.g);
						hit.material.color.x = diffuse.x;
						hit.material.color.y = diffuse.y;
						hit.material.color.z = diffuse.z;
						if (mesh->texture_material[material_texture_idx].user_ReflectionType >= 0)
						{
							//hit.material.color = Color(1,0,0);	//<== Debug用
							hit.material.reflection_type = (ReflectionType)mesh->texture_material[material_texture_idx].user_ReflectionType;
							hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2];

							hit.material.roughness = mesh->texture_material[material_texture_idx].roughness;

							if (hit.material.reflection_type == 3)
							{
								hit.material.color = Color(1, 0, 0);	//<== Debug用
								hit.material.ward_brdf.alp_x = mesh->texture_material[material_texture_idx].ward[0];
								hit.material.ward_brdf.alp_y = mesh->texture_material[material_texture_idx].ward[1];
								hit.material.ward_brdf.specular.x = mesh->texture_material[material_texture_idx].specular.g[0];
								hit.material.ward_brdf.specular.y = mesh->texture_material[material_texture_idx].specular.g[1];
								hit.material.ward_brdf.specular.z = mesh->texture_material[material_texture_idx].specular.g[2];

								hit.material.ward_brdf.diffuse.x = mesh->texture_material[material_texture_idx].diffuse.g[0];
								hit.material.ward_brdf.diffuse.y = mesh->texture_material[material_texture_idx].diffuse.g[1];
								hit.material.ward_brdf.diffuse.z = mesh->texture_material[material_texture_idx].diffuse.g[2];
							}
						}
					}
				}

				if (((Triangle*)this)->material()->texture)
				{
					int j = (((Triangle*)this)->material()->texture->W() - 1)*hit.u;
					int i = (((Triangle*)this)->material()->texture->H() - 1)*hit.v;

					//const Color& diffuse = material.color;
					//hit.material.color.x = diffuse.x*(double)material.texture->cell(i,j).r/255.0;
					//hit.material.color.y = diffuse.y*(double)material.texture->cell(i,j).g/255.0;
					//hit.material.color.z = diffuse.z*(double)material.texture->cell(i,j).b/255.0;
					hit.material.color.x = (double)((Triangle*)this)->material()->texture->cell(i, j).r / 255.0;
					hit.material.color.y = (double)((Triangle*)this)->material()->texture->cell(i, j).g / 255.0;
					hit.material.color.z = (double)((Triangle*)this)->material()->texture->cell(i, j).b / 255.0;
				}
				if (((Triangle*)this)->material()->IBL())
				{
					int j = (((Triangle*)this)->material()->IBL_W() - 1)*hit.u;
					int i = (((Triangle*)this)->material()->IBL_H() - 1)*hit.v;

					Material* ma = const_cast<Material*>(((Triangle*)this)->material());
					hit.material.color = ma->IBL_Color(i, j);
					hit.material.emission = ((Triangle*)this)->material()->ibl_texture_coef*ma->IBL_Color(i, j);
				}
				if (((Triangle*)this)->material()->bump_texture)
				{
					int j = (((Triangle*)this)->material()->bump_texture->W() - 1)*hit.u;
					int i = (((Triangle*)this)->material()->bump_texture->H() - 1)*hit.v;

					Color c((double)((Triangle*)this)->material()->bump_texture->cell(i, j).r,
						(double)((Triangle*)this)->material()->bump_texture->cell(i, j).g,
						(double)((Triangle*)this)->material()->bump_texture->cell(i, j).b);
					c = c / 255.0;

					const double cc = c.length();
					hit.material.color = hit.material.color * cc;
				}
			}
			if (hit.material.reflection_type == -1)
			{
				return false;
			}
			if (status)
			{
				*hitpoint = hit;
			}
			return status;
		}

		void CalcArea()
		{
		}
		Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0)
		{}

		Entity* ConstructBVH()
		{
			return this;
		}
		Entity* ConstructQBVH()
		{
			return this;
		}

		void CreateBoundingBox()
		{
			Vector3d min(PS_INF, PS_INF, PS_INF);
			Vector3d max(-PS_INF, -PS_INF, -PS_INF);
			for (int i = 0; i < 3; i++)
			{
				Vector3d& v = mesh->vertex(vert_id[i]);
				if (min.x > v.x) min.x = v.x;
				if (min.y > v.y) min.y = v.y;
				if (min.z > v.z) min.z = v.z;

				if (max.x < v.x) max.x = v.x;
				if (max.y < v.y) max.y = v.y;
				if (max.z < v.z) max.z = v.z;
			}
			Vector3d eps(0.01, 0.01, 0.01);
			boundingBox = BoundingBox(min - eps, max + eps);
		}

		void MatrixTransformation(std::vector<Matrix4x4>& matrix)
		{
		}
	};

#ifdef USE_STXXL
	typedef 	stxxl::VECTOR_GENERATOR<Triangle>::result					Triangle2;
	typedef 	stxxl::VECTOR_GENERATOR<Triangle>::result::iterator			Entity2;
	typedef		stxxl::VECTOR_GENERATOR<Entity2>::result					Entity2List;
#else
	typedef		std::vector<Triangle>				Triangle2;
	typedef		std::vector<Triangle>::iterator		Entity2;
	typedef		std::vector<Entity2>				Entity2List;
#endif

};
#endif
