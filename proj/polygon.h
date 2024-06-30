#ifndef _POLYGON_H_
#define _POLYGON_H_

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

class Triangle: public Entity
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

		Vector3d P = cross(ray.dir,edge2);
		double det = dot(P,edge1);

		if ( det > PS_EPS) 
		{
			// solve u
			Vector3d T(ray.org - (vert[0]));
			u = dot(P,T);

			if (u>0 && u< det) 
			{
				// solve v
				Vector3d Q = cross(T,edge1);
				v = dot(Q,ray.dir);

				if (v>0 && u+v<det) 
				{
					t = dot(Q,edge2) / det;

					if (t>PS_EPS) 
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
		double det,inv_det;

		/* find Vector3dtors for two edges sharing vert0 */
		const Vector3d edge1 = mesh->vertex(vert_id[1]) - mesh->vertex(vert_id[0]);
		const Vector3d edge2 = mesh->vertex(vert_id[2]) - mesh->vertex(vert_id[0]);

		/* begin calculating determinant - also used to calculate U parameter */
		pVector3d = cross( ray.dir, edge2);

		/* if determinant is near zero, ray lies in plane of triangle */
		det = dot(edge1, pVector3d);

		/* calculate distance from vert0 to ray origin */
		tVector3d =  ray.org - mesh->vertex(vert_id[0]);
   
		qVector3d = cross( tVector3d, edge1);
      
		const double eps = PS_EPS*0.1;
		if (det > eps)
		{
			inv_det = 1.0 / det;
			u = dot(tVector3d, pVector3d);
			if (u < 0.0 || u > det)	return 0;
            
			/* calculate V parameter and test bounds */
			v = dot(ray.dir, qVector3d);
			if (v < 0.0 || u + v > det)	return 0;
      
		}
		else if (det < -eps)
		{
			inv_det = 1.0 / det;
			/* calculate U parameter and test bounds */
			u = dot(tVector3d, pVector3d);
			if (u > 0.0 || u < det)	return 0;
      
			/* calculate V parameter and test bounds */
			v = dot(ray.dir, qVector3d) ;
			if (v > 0.0 || u + v < det)	return 0;
		}
		else return 0;  /* ray is parallell to the plane of the triangle */

		t = dot(edge2, qVector3d) * inv_det;
		u *= inv_det;
		v *= inv_det;
		w = 1.0 - u -v;

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
		double alp = 1.0;	//αチャンネル

		if ( stat == 1 && t > PS_EPS)
		{
			//printf("%d %f\n", i, t);
			status = true;
			hit.distance = t;
			hit.position = ray.org + ray.dir*t;
			hit.normal   = mesh->fnormal(normal_id);
			hit.material = *((Triangle*)this)->material();

			//メッシュ全体をスムーズ化
			int smooth = mesh->smooth;

			if ( mesh->texture_material.size() && material_texture_idx >= 0)
			{
				//この部分だけスムーズ化のON,OFFが指定されている
				int s = mesh->texture_material[material_texture_idx].smooth;
				if ( s )
				{
					if ( s == 1 ) smooth = 1;
					else smooth = 0;
				}
			}
			if (smooth)
			{
				const Vector3d vt[3] ={
					mesh->normal(vert_normal_id[0]),
					mesh->normal(vert_normal_id[1]),
					mesh->normal(vert_normal_id[2])
					};

				//const Vector3d normal2 = normalize(vt[0] + (vt[1]-vt[0])*u + (vt[2]-vt[0])*v);
				const Vector3d normal2 = normalize(vt[0]*w + vt[1]*u + vt[2]*v);	//上の式と同じ
				
				hit.normal = normal2;
			}
			if ( normal_vector_inverse < 0)
			{
				hit.InversNormal();
			}

			if ( mesh->uv.size() )
			{
				//hit.u = mesh->uv[vtx[0]][0]*w + mesh->uv[vtx[1]][0]*u + mesh->uv[vtx[2]][0]*v; 
				//hit.v = mesh->uv[vtx[0]][1]*w + mesh->uv[vtx[1]][1]*u + mesh->uv[vtx[2]][1]*v;
				hit.u = mesh->uv[texture_idx[0]].g[0] * w + mesh->uv[texture_idx[1]].g[0] * u + mesh->uv[texture_idx[2]].g[0] * v;
				hit.v = mesh->uv[texture_idx[0]].g[1] * w + mesh->uv[texture_idx[1]].g[1] * u + mesh->uv[texture_idx[2]].g[1] * v;
			}

			if ( mesh->set_color )
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

			//マテリアルの設定がある場合(newmtlキーワード)
			if ( mesh->texture_material.size() && material_texture_idx >= 0)
			{
				if (mesh->texture_material[material_texture_idx].diffuse.g[3] >= 0.0)
				{
					alp = mesh->texture_material[material_texture_idx].diffuse.g[3];
				}

				//テクスチャが設定されている
				if ( mesh->texture_material[material_texture_idx].textureImage )
				{
					BitMap* image = (BitMap*)mesh->texture_material[material_texture_idx].textureImage;
				
					double dmy;
					double uu = hit.u * mesh->texture_material[material_texture_idx].texture_repeat;
					double vv = hit.v * mesh->texture_material[material_texture_idx].texture_repeat;

					if ( uu < 0 || uu > 1 ) uu = modf(uu, &dmy);
					if ( vv < 0 || vv > 1 ) vv = modf(vv, &dmy);
					if ( uu < 0 ) uu = 1.0 - fabs(uu);
					if ( vv < 0 ) vv = 1.0 - fabs(vv);

					int j = (image->W()-1)*uu;
					int i = (image->H()-1)*vv;

					//const Color& diffuse = Color(mesh->texture_material[material_texture_idx].diffuse);
					//hit.material.color.x = diffuse.x*(double)image->cell(i,j).r/255.0;
					//hit.material.color.y = diffuse.y*(double)image->cell(i,j).g/255.0;
					//hit.material.color.z = diffuse.z*(double)image->cell(i,j).b/255.0;

					const Color& diffuse = ((Triangle*)this)->material()->color;
					hit.material.color.x = diffuse.x*(double)image->cell(i,j).r/255.0;
					hit.material.color.y = diffuse.y*(double)image->cell(i,j).g/255.0;
					hit.material.color.z = diffuse.z*(double)image->cell(i,j).b/255.0;

					alp = (double)image->cell(i,j).alp/255.0;

					hit.material.emission = ((Triangle*)this)->material()->emission;
					if (mesh->texture_material[material_texture_idx].user_ReflectionType >= 0)
					{

						hit.material.emission.x = mesh->texture_material[material_texture_idx].emission.g[0] * (double)image->cell(i, j).r / 255.0;;
						hit.material.emission.y = mesh->texture_material[material_texture_idx].emission.g[1] * (double)image->cell(i, j).g / 255.0;;
						hit.material.emission.z = mesh->texture_material[material_texture_idx].emission.g[2] * (double)image->cell(i, j).b / 255.0;;

						//hit.material.color = Color(1,0,0);	//<== Debug用
						hit.material.reflection_type = (ReflectionType)mesh->texture_material[material_texture_idx].user_ReflectionType;
						hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0] * (double)image->cell(i, j).r / 255.0;;
						hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1] * (double)image->cell(i, j).g / 255.0;;
						hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2] * (double)image->cell(i, j).b / 255.0;;

						hit.material.refractive_index = mesh->texture_material[material_texture_idx].refraction_index;
						hit.material.roughness = mesh->texture_material[material_texture_idx].roughness;

						if (hit.material.reflection_type == REFLECTION_TYPE_SPECULAR)
						{
							hit.material.color.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.color.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.color.z = mesh->texture_material[material_texture_idx].specular.g[2];
						}
						if (hit.material.reflection_type == REFLECTION_TYPE_REFRACTION)
						{
							hit.material.color.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.color.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.color.z = mesh->texture_material[material_texture_idx].specular.g[2];
						}
						if (hit.material.reflection_type == REFLECTION_TYPE_WARD_BRDF)
						{
							hit.material.brdfParameter.ward_brdf.alp_x = mesh->texture_material[material_texture_idx].ward[0];
							hit.material.brdfParameter.ward_brdf.alp_y = mesh->texture_material[material_texture_idx].ward[1];
							hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0] * (double)image->cell(i, j).r / 255.0;;
							hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1] * (double)image->cell(i, j).g / 255.0;;
							hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2] * (double)image->cell(i, j).b / 255.0;;

							hit.material.color.x = mesh->texture_material[material_texture_idx].diffuse.g[0] * (double)image->cell(i, j).r / 255.0;;
							hit.material.color.y = mesh->texture_material[material_texture_idx].diffuse.g[1] * (double)image->cell(i, j).g / 255.0;;
							hit.material.color.z = mesh->texture_material[material_texture_idx].diffuse.g[2] * (double)image->cell(i, j).b / 255.0;;
						}
						if (hit.material.reflection_type == REFLECTION_TYPE_PHONG_BRDF)
						{
							hit.material.brdfParameter.phong_brdf.specular_exponent = mesh->texture_material[material_texture_idx].shininess;
							hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0] * (double)image->cell(i, j).r / 255.0;;
							hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1] * (double)image->cell(i, j).g / 255.0;;
							hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2] * (double)image->cell(i, j).b / 255.0;;

							hit.material.color.x = mesh->texture_material[material_texture_idx].diffuse.g[0] * (double)image->cell(i, j).r / 255.0;;
							hit.material.color.y = mesh->texture_material[material_texture_idx].diffuse.g[1] * (double)image->cell(i, j).g / 255.0;;
							hit.material.color.z = mesh->texture_material[material_texture_idx].diffuse.g[2] * (double)image->cell(i, j).b / 255.0;;
						}
					}
				}
				else
				{
					//テクスチャは設定されていない
					const Color& diffuse = Color(mesh->texture_material[material_texture_idx].diffuse.g);
					const Color& specular = Color(mesh->texture_material[material_texture_idx].specular.g);
					hit.material.color = diffuse;
					hit.material.specular = specular;

					if (hit.material.reflection_type == REFLECTION_TYPE_DIFFUSE)
					{
						hit.material.specular = Color(0.0);
					}
					if (mesh->texture_material[material_texture_idx].diffuse.g[3] >= 0.0)
					{
						alp = mesh->texture_material[material_texture_idx].diffuse.g[3];
					}

					hit.material.emission = ((Triangle*)this)->material()->emission;

					if (mesh->texture_material[material_texture_idx].user_ReflectionType >= 0)
					{
						hit.material.emission.x = mesh->texture_material[material_texture_idx].emission.g[0];
						hit.material.emission.y = mesh->texture_material[material_texture_idx].emission.g[1];
						hit.material.emission.z = mesh->texture_material[material_texture_idx].emission.g[2];

						//hit.material.color = Color(1,0,0);	//<== Debug用
						hit.material.reflection_type = (ReflectionType)mesh->texture_material[material_texture_idx].user_ReflectionType;
						hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0];
						hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1];
						hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2];

						hit.material.refractive_index = mesh->texture_material[material_texture_idx].refraction_index;
						hit.material.roughness = mesh->texture_material[material_texture_idx].roughness;

						if (hit.material.reflection_type == REFLECTION_TYPE_SPECULAR)
						{
							hit.material.color.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.color.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.color.z = mesh->texture_material[material_texture_idx].specular.g[2];
						}
						if (hit.material.reflection_type == REFLECTION_TYPE_REFRACTION)
						{
							hit.material.color.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.color.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.color.z = mesh->texture_material[material_texture_idx].specular.g[2];
						}
						if (hit.material.reflection_type == REFLECTION_TYPE_WARD_BRDF)
						{
							hit.material.brdfParameter.ward_brdf.alp_x = mesh->texture_material[material_texture_idx].ward[0];
							hit.material.brdfParameter.ward_brdf.alp_y = mesh->texture_material[material_texture_idx].ward[1];
							hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2];

							hit.material.color.x = mesh->texture_material[material_texture_idx].diffuse.g[0];
							hit.material.color.y = mesh->texture_material[material_texture_idx].diffuse.g[1];
							hit.material.color.z = mesh->texture_material[material_texture_idx].diffuse.g[2];
						}
						if (hit.material.reflection_type == REFLECTION_TYPE_PHONG_BRDF)
						{
							hit.material.brdfParameter.phong_brdf.specular_exponent = mesh->texture_material[material_texture_idx].shininess;
							hit.material.specular.x = mesh->texture_material[material_texture_idx].specular.g[0];
							hit.material.specular.y = mesh->texture_material[material_texture_idx].specular.g[1];
							hit.material.specular.z = mesh->texture_material[material_texture_idx].specular.g[2];

							hit.material.color.x = mesh->texture_material[material_texture_idx].diffuse.g[0];
							hit.material.color.y = mesh->texture_material[material_texture_idx].diffuse.g[1];
							hit.material.color.z = mesh->texture_material[material_texture_idx].diffuse.g[2];
						}
					}
				}

				//バンプマッピングの指定がある
				if ( mesh->texture_material[material_texture_idx].bump_textureImage )
				{
					BitMap* image = (BitMap*)mesh->texture_material[material_texture_idx].bump_textureImage;
					double dmy;
					double uu = hit.u * mesh->texture_material[material_texture_idx].texture_repeat;
					double vv = hit.v * mesh->texture_material[material_texture_idx].texture_repeat;

					if ( uu < 0 || uu > 1 ) uu = modf(uu, &dmy);
					if ( vv < 0 || vv > 1 ) vv = modf(vv, &dmy);
					if ( uu < 0 ) uu = 1.0 - fabs(uu);
					if ( vv < 0 ) vv = 1.0 - fabs(vv);

					int j = (image->W()-1)*uu;
					int i = (image->H()-1)*vv;
#if 0
					Color c((double)image->cell(i,j).r,
							(double)image->cell(i,j).g,
							(double)image->cell(i,j).b);
					c = c / 255.0;

					const double cc = c.length();
					hit.material.color = hit.material.color * cc;
#else

					//Mathematics for 3D Game Programming and Computer Graphics, Third Edition
					// begin-CalculateTangent
					Vector3d v1 = mesh->vertex(vert_id[0]);
					Vector3d v2 = mesh->vertex(vert_id[1]);
					Vector3d v3 = mesh->vertex(vert_id[2]);
					Vector3d w1 = mesh->texcoord(texture_idx[0]);
					Vector3d w2 = mesh->texcoord(texture_idx[1]);
					Vector3d w3 = mesh->texcoord(texture_idx[2]);

					double x1 = v2.x - v1.x;
					double x2 = v3.x - v1.x;
					double y1 = v2.y - v1.y;
					double y2 = v3.y - v1.y;
					double z1 = v2.z - v1.z;
					double z2 = v3.z - v1.z;

					double s1 = w2.x - w1.x;
					double s2 = w3.x - w1.x;
					double t1 = w2.y - w1.y;
					double t2 = w3.y - w1.y;

					double ww = (s1 * t2 - s2 * t1);
					double r = 0.0;
					if ( fabs(ww) > 1.0e-12 ) r = 1.0 / (s1 * t2 - s2 * t1);
					Vector3d sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r,	(t2 * z1 - t1 * z2) * r);
					Vector3d tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r,	(s1 * z2 - s2 * z1) * r);
					// end-CalculateTangent

					//compute them using central difference
					int iu0 = j - 1;
					int iu1 = j + 1;

					int iv0 = i - 1;
					int iv1 = i + 1;

					if (iu0 < 0) iu0 = image->W() - 1;
					if (iv0 < 0) iv0 = image->H() - 1;
					if (iu1 < 0) iu1 = image->W() - 1;
					if (iv1 < 0) iv1 = image->H() - 1;
					if (iu0 >= image->W()) iu0 = 0;
					if (iu1 >= image->W()) iu1 = 0;
					if (iv0 >= image->H()) iv0 = 0;
					if (iv1 >= image->H()) iv1 = 0;

					Color px0(
						image->cell(i, iu0).r,
						image->cell(i, iu0).g,
						image->cell(i, iu0).b);

					Color px1(
						image->cell(i, iu1).r,
						image->cell(i, iu1).g,
						image->cell(i, iu1).b);

					Color py0(
						image->cell(iv0, j).r,
						image->cell(iv0, j).g,
						image->cell(iv0, j).b);

					Color py1(
						image->cell(iv1, j).r,
						image->cell(iv1, j).g,
						image->cell(iv1, j).b);

					//グレースケール（白黒画像）に変換
					Color YCrCb_Y(0.299, 0.587, 0.114);

					//central difference
					double x_gradient = dot((px1 - px0), YCrCb_Y)*0.5;
					double y_gradient = dot((py1 - py0), YCrCb_Y)*0.5;

					//法線をすでに反転している
					if (normal_vector_inverse < 0)
					{
						//一時的に元に戻す
						hitpoint->normal = hitpoint->normal*-1.0;
					}
					//BumpMapping new Normal vector 
					hitpoint->bump_new_normal = hitpoint->normal + cross(hitpoint->normal, tdir)*x_gradient - cross(hitpoint->normal, sdir)*y_gradient;
					hitpoint->bump = 1;

					if (normal_vector_inverse < 0)
					{
						//反転する
						hitpoint->InversNormal();
					}

#endif
				}

			}

			if ( ((Triangle*)this)->material()->texture )
			{
				int j = (((Triangle*)this)->material()->texture->W()-1)*hit.u;
				int i = (((Triangle*)this)->material()->texture->H()-1)*hit.v;

				//const Color& diffuse = material.color;
				//hit.material.color.x = diffuse.x*(double)material.texture->cell(i,j).r/255.0;
				//hit.material.color.y = diffuse.y*(double)material.texture->cell(i,j).g/255.0;
				//hit.material.color.z = diffuse.z*(double)material.texture->cell(i,j).b/255.0;
				hit.material.color.x = (double)((Triangle*)this)->material()->texture->cell(i,j).r/255.0;
				hit.material.color.y = (double)((Triangle*)this)->material()->texture->cell(i,j).g/255.0;
				hit.material.color.z = (double)((Triangle*)this)->material()->texture->cell(i,j).b/255.0;

				alp = (double)((Triangle*)this)->material()->texture->cell(i,j).alp/255.0;
			}
			if ( ((Triangle*)this)->material()->IBL() )
			{
				int j = (((Triangle*)this)->material()->IBL_W()-1)*hit.u;
				int i = (((Triangle*)this)->material()->IBL_H()-1)*hit.v;

				Material* ma = const_cast<Material*>(((Triangle*)this)->material());
				hit.material.color = ma->IBL_Color(i,j);
				hit.material.emission = ((Triangle*)this)->material()->ibl_texture_coef*ma->IBL_Color(i,j);
			}
			if ( ((Triangle*)this)->material()->bump_texture )
			{
				int j = (((Triangle*)this)->material()->bump_texture->W()-1)*hit.u;
				int i = (((Triangle*)this)->material()->bump_texture->H()-1)*hit.v;

#if 0
				Color c((double)((Triangle*)this)->material()->bump_texture->cell(i,j).r,
					    (double)((Triangle*)this)->material()->bump_texture->cell(i,j).g,
					    (double)((Triangle*)this)->material()->bump_texture->cell(i,j).b);
				c = c / 255.0;

				const double cc = c.length();
				hit.material.color = hit.material.color * cc;
#else

				//Mathematics for 3D Game Programming and Computer Graphics, Third Edition
				// begin-CalculateTangent
				Vector3d v1 = mesh->vertex(vert_id[0]);
				Vector3d v2 = mesh->vertex(vert_id[1]);
				Vector3d v3 = mesh->vertex(vert_id[2]);
				Vector3d w1 = mesh->texcoord(texture_idx[0]);
				Vector3d w2 = mesh->texcoord(texture_idx[1]);
				Vector3d w3 = mesh->texcoord(texture_idx[2]);

				double x1 = v2.x - v1.x;
				double x2 = v3.x - v1.x;
				double y1 = v2.y - v1.y;
				double y2 = v3.y - v1.y;
				double z1 = v2.z - v1.z;
				double z2 = v3.z - v1.z;

				double s1 = w2.x - w1.x;
				double s2 = w3.x - w1.x;
				double t1 = w2.y - w1.y;
				double t2 = w3.y - w1.y;

				double ww = (s1 * t2 - s2 * t1);
				double r = 0.0;
				if ( fabs(ww) > 1.0e-12 ) r = 1.0 / (s1 * t2 - s2 * t1);
				Vector3d sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r,	(t2 * z1 - t1 * z2) * r);
				Vector3d tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r,	(s1 * z2 - s2 * z1) * r);
				// end-CalculateTangent



				//BumpMapping new Normal vector 
				int ii = i-1;
				int jj = j;

				if (i - 1 < 0) ii = 0;
				Color yy1((double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).r,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).g,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).b);
				
				ii = i + 1;
				if (ii >= (((Triangle*)this)->material()->bump_texture->H())) ii = ((Triangle*)this)->material()->bump_texture->H() - 1;
				Color yy2((double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).r,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).g,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).b);
				
				ii = i;
				jj = j - 1;
				if (j - 1 < 0) jj = 0;
				Color xx1((double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).r,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).g,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).b);

				jj = j + 1;
				if (jj >= (((Triangle*)this)->material()->bump_texture->W())) jj = ((Triangle*)this)->material()->bump_texture->W() - 1;
				Color xx2((double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).r,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).g,
					(double)((Triangle*)this)->material()->bump_texture->cell(ii, jj).b);

				double x_gradient = dot((xx2 - xx1),Vector3d(1,1,1))*0.5/3.0;
				double y_gradient = dot((yy2 - yy1),Vector3d(1,1,1))*0.5/3.0;

				//法線をすでに反転している
				if (normal_vector_inverse < 0)
				{
					hitpoint->normal = normalize(hitpoint->normal*-1.0 + cross(tdir, hitpoint->normal)*x_gradient - cross(sdir, hitpoint->normal)*y_gradient)*-1.0;
				}
				else
				{
					hitpoint->normal = normalize(hitpoint->normal + cross(tdir, hitpoint->normal)*x_gradient - cross(sdir, hitpoint->normal)*y_gradient);
				}
#endif
			}
		}
		if (hit.material.reflection_type == -1)
		{
			return false;
		}

		//αチャンネル
		if ( alp == 0.0)
		{
			return false;
		}
		if (alp < 1.0)
		{
			if (alp < rnd->next01()) return false;
		}
		if (status)
		{
			*hitpoint = hit;
		}
		return status;
	}

	void CalcArea()
	{
		const Vector3d edge1 = mesh->vertex(vert_id[1]) - mesh->vertex(vert_id[0]);
		const Vector3d edge2 = mesh->vertex(vert_id[2]) - mesh->vertex(vert_id[0]);

		area = cross(edge1, edge2).length()*0.5;
	}

	Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0)
	{
#if 0
		const double u = rnd->next01();
		const double v = rnd->next01();
		const double w = 1.0 - u - v;
#else
		const double u1 = rnd->next01();
		const double u2 = rnd->next01();

		const double z = sqrt(u1);
		const double u = 1.0 - z;
		const double v =(1.0 -  u2)*z;
		const double w = u2*z;
#endif
		if (nrm)
		{
			if (normal_vector_inverse < 0)
			{
				*nrm = -1.0*normalize(mesh->normal(vert_normal_id[0])*u + mesh->normal(vert_normal_id[1])*v + mesh->normal(vert_normal_id[2])*w);
			}
			else
			{
				*nrm = normalize(mesh->normal(vert_normal_id[0])*u + mesh->normal(vert_normal_id[1])*v + mesh->normal(vert_normal_id[2])*w);
			}
		}
		return mesh->vertex(vert_id[0])*u + mesh->vertex(vert_id[1])*v + mesh->vertex(vert_id[2])*w;
	}

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
		for ( int i = 0; i < 3; i++ )
		{
			Vector3d& v = mesh->vertex( vert_id[i] );
			if ( min.x > v.x ) min.x = v.x;
			if ( min.y > v.y ) min.y = v.y;
			if ( min.z > v.z ) min.z = v.z;

			if ( max.x < v.x ) max.x = v.x;
			if ( max.y < v.y ) max.y = v.y;
			if ( max.z < v.z ) max.z = v.z;
		}
		Vector3d eps(0.01, 0.01, 0.01);
		boundingBox = BoundingBox(min-eps, max+eps);
	}

	void MatrixTransformation(std::vector<Matrix4x4>& matrix)
	{
	}
};

class Polygon: public Entity 
{
public:

	::Obj* mesh;		//元のメッシュデータ

	Polygon(Obj* obj, int m)
	{
		type = ENTITY_TYPE_POLYGON;
		mesh = obj;
		mesh->smooth = 0;
		material_id = m;
	}

	~Polygon()
	{
		if ( triangles.size() )
		{
			const int sz = triangles.size();
			for ( int i = 0; i < sz; i++ ) delete triangles[i];
		}
		triangles.clear();
		
		for ( int i = 0; i < mesh->texture_material.size(); i++ )
		{
			if ( mesh->texture_material[i].textureImage )
			{
				delete (BitMap*)mesh->texture_material[i].textureImage;
				mesh->texture_material[i].textureImage = 0;
			}
			if ( mesh->texture_material[i].bump_textureImage )
			{
				delete (BitMap*)mesh->texture_material[i].bump_textureImage;
				mesh->texture_material[i].bump_textureImage = 0;
			}
		}
		delete mesh;
	}

	void MatrixTransformation(std::vector<Matrix4x4>& matrix)
	{
		Vector3d org(0,0,0);
		for ( int j =0; j < matrix.size() ; j++ )
		{
			for ( int i = 0; i < mesh->nv; i++ )
			{
				Vector3d v = mesh->vertex(i);
				v = matrix[j] * v;
				mesh->vert[i].g[0] = v.x;
				mesh->vert[i].g[1] = v.y;
				mesh->vert[i].g[2] = v.z;
			}
			for ( int i = 0; i < mesh->nf; i++ )
			{
				Matrix4x4 g = matrix[j].vector_coordinate_transformation();
				Vector3d fn = mesh->fnormal(i);
				fn = g * fn;

				fn = normalize(fn);
				mesh->fnorm[i].g[0] = fn.x;
				mesh->fnorm[i].g[1] = fn.y;
				mesh->fnorm[i].g[2] = fn.z;
			}
			for ( int i = 0; i < mesh->nvn; i++ )
			{
				Matrix4x4 g = matrix[j].vector_coordinate_transformation();
				Vector3d n = mesh->normal(i);
				n = g * n;
				n = normalize(n);
				mesh->norm[i].g[0] = n.x;
				mesh->norm[i].g[1] = n.y;
				mesh->norm[i].g[2] = n.z;
			}
		}
		//mesh->CalcNormalVector();
	}

	void translate(double dx, double dy, double dz)
	{
#pragma omp parallel for schedule(dynamic, 1)
		for ( int i = 0; i < mesh->nv; i++ )
		{
			mesh->vert[i].g[0] += dx;
			mesh->vert[i].g[1] += dy;
			mesh->vert[i].g[2] += dz;
		}
	}
	void scale(double sx, double sy, double sz)
	{
#pragma omp parallel for schedule(dynamic, 1)
		for ( int i = 0; i < mesh->nv; i++ )
		{
			mesh->vert[i].g[0] *= sx;
			mesh->vert[i].g[1] *= sy;
			mesh->vert[i].g[2] *= sz;
		}
	}

	bool intersect(const Ray &ray, IntersectionPos *hitpoint) const 
	{
		double org[3] = {ray.org.x, ray.org.y, ray.org.z}; 
		double dir[3] = {ray.dir.x, ray.dir.y, ray.dir.z};

		//printf("\n\npolygon-Ray-Triangle Intersection\n");
		bool status = false;

		IntersectionPos hit;
		const int sz = mesh->nf;
		int index = -1;
		for ( int i = 0; i < sz; i++ )
		{
			if ( !triangles[i]->boundingBox.CheckIntersection(ray))
			{
				continue;
			}
			double vert[3][3];

			vert[0][0] = mesh->vert[mesh->face[i].g[0]].g[0];
			vert[0][1] = mesh->vert[mesh->face[i].g[0]].g[1];
			vert[0][2] = mesh->vert[mesh->face[i].g[0]].g[2];

			vert[1][0] = mesh->vert[mesh->face[i].g[1]].g[0];
			vert[1][1] = mesh->vert[mesh->face[i].g[1]].g[1];
			vert[1][2] = mesh->vert[mesh->face[i].g[1]].g[2];

			vert[2][0] = mesh->vert[mesh->face[i].g[2]].g[0];
			vert[2][1] = mesh->vert[mesh->face[i].g[2]].g[1];
			vert[2][2] = mesh->vert[mesh->face[i].g[2]].g[2];

			double t, u, v;
			int stat = intersect_triangle3(org, dir, vert[0], vert[1], vert[2], &t, &u, &v);
			if ( stat == 1 && hit.distance > t && t > PS_EPS)
			{
				//printf("%d %f\n", i, t);
				status = true;
				hit.distance = t;
				hit.position = ray.org + t * ray.dir;
				hit.normal = Vector3d(mesh->fnorm[i].g[0], mesh->fnorm[i].g[1], mesh->fnorm[i].g[2]);
				hit.material = *triangles[i]->material();

				if ( normal_vector_inverse < 0)
				{
					hit.normal =  hit.normal*-1.0;
				}
				index = i;
			}
		}
		if ( status )
		{
			*hitpoint = hit;
			//if ( this->bvh )
			//{
			//	printf("総あたり    <%d> %f %f %f  %f %f %f\n", index,
			//					hit.position.x,
			//					hit.position.y,
			//					hit.position.z,
			//					hit.normal.x,
			//					hit.normal.y,
			//					hit.normal.z);
			//				fflush(stdout);
			//}
		}
		return status;
	}

	void CreateBoundingBox()
	{
		Vector3d min(PS_INF, PS_INF, PS_INF);
		Vector3d max(-PS_INF, -PS_INF, -PS_INF);
		for ( int i = 0; i < mesh->nv; i++ )
		{
			if ( min.x > mesh->vert[i].g[0] ) min.x = mesh->vert[i].g[0];
			if ( min.y > mesh->vert[i].g[1] ) min.y = mesh->vert[i].g[1];
			if ( min.z > mesh->vert[i].g[2] ) min.z = mesh->vert[i].g[2];
			if ( max.x < mesh->vert[i].g[0] ) max.x = mesh->vert[i].g[0];
			if ( max.y < mesh->vert[i].g[1] ) max.y = mesh->vert[i].g[1];
			if ( max.z < mesh->vert[i].g[2] ) max.z = mesh->vert[i].g[2];
		}
		Vector3d eps(0.01, 0.01, 0.01);
		boundingBox = BoundingBox(min-eps, max+eps);

		printf("boundingBox=[%f,%f,%f]->[%f,%f,%f]\n", 
			boundingBox.min().x,
			boundingBox.min().y,
			boundingBox.min().z,
			boundingBox.max().x,
			boundingBox.max().y,
			boundingBox.max().z);
	}


	//三角形データリストの作成
#ifdef USE_STXXL
	typedef stxxl::VECTOR_GENERATOR<Entity*>::result			EntityList_t;
	typedef stxxl::VECTOR_GENERATOR<Triangle>::result			TriangleList_t;
#else
	typedef std::vector<Entity*>								EntityList_t;
	typedef std::vector<Triangle>								TriangleList_t;
#endif

	TriangleList_t TriangleList;
	EntityList_t triangles;
	void MakeTrianglesEntity()
	{
		printf("MakeTrianglesEntity start\n");
		fflush(stdout);
		const int sz = mesh->nf;

		TriangleList.clear();
		triangles.clear();

		area = 0.0;
		TriangleList.resize(sz);
		printf("mesh %p triangle %d\n", mesh, sz); fflush(stdout);

		for ( int i = 0; i < sz; i++ )
		{
			Triangle *tri = &(TriangleList[i]);

			tri->mesh = mesh;
			tri->normal_vector_inverse = normal_vector_inverse;
			tri->face_id = i;

			tri->material_id = material_id;
			tri->materialList= materialList;

			tri->vert_id[0] = mesh->face[i].g[0];
			tri->vert_id[1] = mesh->face[i].g[1];
			tri->vert_id[2] = mesh->face[i].g[2];

			tri->normal_id   = i;


			tri->vert_normal_id[0] = mesh->norm_id[i].g[0];
			tri->vert_normal_id[1] = mesh->norm_id[i].g[1];
			tri->vert_normal_id[2] = mesh->norm_id[i].g[2];

			if ( mesh->texture_id.size() )
			{
				tri->texture_idx[0] = mesh->texture_id[i].g[0];
				tri->texture_idx[1] = mesh->texture_id[i].g[1];
				tri->texture_idx[2] = mesh->texture_id[i].g[2];
			}
			if ( mesh->texture_material.size())
			{
				tri->material_texture_idx = mesh->texture_material_id[i];
			}

			//mesh->CalcNormalVector();
			tri->id = i;
			tri->CreateBoundingBox();
			tri->CalcArea();
			tri->media_ignore_front_back = media_ignore_front_back;
			tri->rnd = this->rnd;
			triangles.push_back((Entity*)tri);
			area += tri->area;
		}
		printf("MakeTrianglesEntity end\n");
		fflush(stdout);
	}

	void LoadTexture()
	{
		char drive[_MAX_DRIVE];	// ドライブ名
		char dir[_MAX_DIR];		// ディレクトリ名
		char fname[_MAX_FNAME];	// ファイル名
		char ext[_MAX_EXT];		// 拡張子

		const int sz = mesh->texture_material.size();
		for ( int i = 0; i < sz; i++ )
		{
			mesh->texture_material[i].textureImage = 0;
			char filename[1024];
			strcpy(filename, mesh->texture_material[i].texture_name.c_str());

			if (filename[0] == '\0')
			{
				continue;
			}

			FILE* fp = fopen(filename, "r");
			if ( fp )
			{
				printf("texture=[%s]\n", filename);
				fclose(fp);
			}else
			{
				_splitpath( mesh->materialName.c_str(), drive, dir, fname, ext );
				sprintf(  filename, "%s%s%s", drive, dir, mesh->texture_material[i].texture_name.c_str());
				printf("texture=[%s]\n", filename);
				fp = fopen(filename, "r");
				if ( fp )
				{
					fclose(fp);
				}else
				{
					printf("ファイルが開けません[%s]\n", filename);
					continue ;
				}
				printf("textureImage[%s]Load.\n", filename);
			}
			BitMap* bmp = new BitMap();
			bmp->Read(filename);
			if (bmp->GetImage())
			{
				mesh->texture_material[i].textureImage = (void*)bmp;
			}else
			{
				delete bmp;
			}
		}
	}
	void LoadBumpTexture()
	{
		char drive[_MAX_DRIVE];	// ドライブ名
		char dir[_MAX_DIR];		// ディレクトリ名
		char fname[_MAX_FNAME];	// ファイル名
		char ext[_MAX_EXT];		// 拡張子

		const int sz = mesh->texture_material.size();
		for ( int i = 0; i < sz; i++ )
		{
			mesh->texture_material[i].bump_textureImage = 0;
			char filename[1024];
			strcpy(filename, mesh->texture_material[i].bump_texture_name.c_str());

			if (filename[0] == '\0')
			{
				continue;
			}

			FILE* fp = fopen(filename, "r");
			if ( fp )
			{
				printf("bump_texture=[%s]\n", filename);
				fclose(fp);
			}else
			{
				_splitpath( mesh->materialName.c_str(), drive, dir, fname, ext );
				sprintf(  filename, "%s%s%s", drive, dir, mesh->texture_material[i].bump_texture_name.c_str());
				printf("bump_texture=[%s]\n", filename);
				fp = fopen(filename, "r");
				if ( fp )
				{
					fclose(fp);
				}else
				{
					printf("ファイルが開けません[%s]\n", filename);
					continue ;
				}
			}			
			BitMap* bmp = new BitMap();
			bmp->Read(filename);
			if (bmp->GetImage())
			{
				mesh->texture_material[i].bump_textureImage = (void*)bmp;
			}else
			{
				delete bmp;
			}
		}
	}

	void CalcArea()
	{
		//MakeTrianglesEntityで計算済
	}
	Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0)
	{
		int face_index = (int)(mesh->nf * rnd->Next01());
		if (face_index >= mesh->nf) face_index = mesh->nf - 1;
		Triangle* tri = (Triangle*)triangles[face_index];

		if (face_index_p) *face_index_p = face_index;
		if (nrm)
		{
			*nrm = mesh->fnormal(face_index);
		}
		return tri->randomPoint(nrm, face_index_p);
	}

	Entity* ConstructBVH()
	{
		if (bvh) delete bvh;
		
		MakeTrianglesEntity();
		
		bvh = new BVH();
		
		//bvh->Construct(BVH::CONSTRUCTION_OBJECT_MEDIAN, triangles);
		bvh->Construct(BVH::CONSTRUCTION_OBJECT_SAH, triangles);
		
		return this;
	}
	Entity* ConstructQBVH()
	{
		if (qbvh) delete qbvh;
		
		MakeTrianglesEntity();
		
		qbvh = new QBVH();
		
		qbvh->Construct(triangles);
		
		return this;
	}
};

};

#endif
