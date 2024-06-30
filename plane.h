#ifndef _PLANE_H
#define _PLANE_H

#include "entity.h"

namespace prender {

class Plane: public Entity {
	Vector3d local[2];
	Vector3d p[2];
public:
	Vector3d org;
	Vector3d normal;
	bool shadow;

	Plane(const Vector3d& o, const Vector3d& n,  int m,  const Vector3d* uv=0):org(o),normal(n)
	{
		type = EntityType::ENTITY_TYPE_PLANE;
		material_id = m;
		shadow = false;
		
		if ( !uv )
		{
			Vector3d x(1,0,0);
			Vector3d y(0,1,0);

			local[1] = cross(normal, x);
			if ( local[1].length() < PS_EPS )
			{
				local[0] = normalize( cross(y, normal) );
				local[1] = normalize( cross(normal, local[0]) );
			}else
			{
				local[1] = normalize( local[1] );
				local[0] = normalize( cross(local[1], normal) );
			}
		}else
		{
			local[0] = uv[0];
			local[1] = uv[1];
		}
		p[0] = org + local[0];
		p[1] = org + local[1];
	}

	void CalcArea()
	{
	}
	Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0)
	{
		return Vector3d(0,0,0);
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
		const double s1 =  1000000000.0;
		const double s2 = -1000000000.0;
		Vector3d q[4];
		q[0] = org + s1*local[0] + s1*local[1];
		q[1] = org + s1*local[0] + s2*local[1];
		q[2] = org + s2*local[0] + s1*local[1];
		q[3] = org + s2*local[0] + s2*local[1];

		Vector3d min = q[0];
		Vector3d max = q[0];
		for ( int i = 1; i < 4; i++ )
		{
			if ( min.x > q[i].x ) min.x = q[i].x;
			if ( min.y > q[i].y ) min.y = q[i].y;
			if ( min.z > q[i].z ) min.z = q[i].z;
			if ( max.x < q[i].x ) max.x = q[i].x;
			if ( max.y < q[i].y ) max.y = q[i].y;
			if ( max.z < q[i].z ) max.z = q[i].z;
		}
		boundingBox = BoundingBox(min, max);
	}

	bool intersect(const Ray &ray, IntersectionPos *hitpoint) const 
	{
		Vector3d normal2 = normal;
		if ( normal_vector_inverse < 0 )
		{
			normal2 =  normal*-1.0;
		}

		double d = dot(ray.dir, normal2);

		if ( fabs(d) < PS_EPS14)
		{
			return false;
		}
		hitpoint->distance = dot((org - ray.org), normal2)/d;
		hitpoint->normal   = normal2;
		hitpoint->position = ray.dir*hitpoint->distance + ray.org;
		hitpoint->material = *(((Plane*)this)->material());
	
		if ( hitpoint->distance < PS_EPS )
		{
			return false;
		}
		return true;
	}

	void MatrixTransformation(std::vector<Matrix4x4>& matrix)
	{
		for ( int j =0; j < matrix.size() ; j++ )
		{
			org = matrix[j] * org;
			p[0] = matrix[j] * p[0];
			p[1] = matrix[j] * p[1];

			Matrix4x4 g = matrix[j].vector_coordinate_transformation();
			local[0] = g*local[0];
			local[1] = g*local[1];
			normal = g*normal;
		}
	}

};

class UVPlane: public Entity {

public:
	Vector3d p[2];
	Vector3d org;
	Vector3d u_axis;
	Vector3d v_axis;
	double u_length;
	double v_length;
	Vector3d normal;
	Plane* base_plane;
	bool shadow;

	bool circle;

	UVPlane()
	{}
	UVPlane(const Vector3d& org, const Vector3d& X, const Vector3d& Y,  int m)
		:org(org),u_axis(X),v_axis(Y)
	{
		type = EntityType::ENTITY_TYPE_UVPLANE;
		material_id = m;
		shadow = false;
		circle = false;

		u_length = u_axis.length();
		v_length = v_axis.length();
		u_axis = normalize(u_axis);
		v_axis = normalize(v_axis);
		normal = normalize( cross(u_axis, v_axis) );

		//printf("%f %f %f\n", u_axis.x, u_axis.y, u_axis.z);
		//printf("%f %f %f\n", v_axis.x, v_axis.y, v_axis.z);
		//printf("%f %f %f\n", normal.x, normal.y, normal.z);
		Vector3d uv[2] = {u_axis, v_axis};
		base_plane = new Plane(org, normal,m, uv);
		p[0] = org + u_axis;
		p[1] = org + v_axis;
	}

	~UVPlane()
	{
		delete base_plane;
	}

	void CalcArea()
	{
		area = (u_axis*u_length).length()*(v_axis*v_length).length()*(1.0 - pow(dot(u_axis,v_axis),2.0) );
	}

	Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0)
	{
		const double t = v_length * rnd->next01();
		const double s = u_length * rnd->next01();

		if (nrm) *nrm = normal;
		return org + u_axis*s + v_axis*t;
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
		const double s1 =  0.0;
		const double s2 =  1.0;
		Vector3d q[4];
		q[0] = org + u_length*s1*u_axis + v_length*s1*v_axis;
		q[1] = org + u_length*s1*u_axis + v_length*s2*v_axis;
		q[2] = org + u_length*s2*u_axis + v_length*s1*v_axis;
		q[3] = org + u_length*s2*u_axis + v_length*s2*v_axis;

		Vector3d min = q[0];
		Vector3d max = q[0];
		for ( int i = 1; i < 4; i++ )
		{
			if ( min.x > q[i].x ) min.x = q[i].x;
			if ( min.y > q[i].y ) min.y = q[i].y;
			if ( min.z > q[i].z ) min.z = q[i].z;
			if ( max.x < q[i].x ) max.x = q[i].x;
			if ( max.y < q[i].y ) max.y = q[i].y;
			if ( max.z < q[i].z ) max.z = q[i].z;
		}
		boundingBox = BoundingBox(min, max);
	}

	bool intersect(const Ray &ray, IntersectionPos *hitpoint) const 
	{
		base_plane->normal_vector_inverse = normal_vector_inverse;
		hitpoint->material = *(((UVPlane*)this)->material());

		if ( !base_plane->intersect(ray, hitpoint))
		{
			return false;
		}
		
		double alp =1.0;

		hitpoint->u = dot(Vector3d(hitpoint->position - org), u_axis) / u_length;
		hitpoint->v = dot(Vector3d(hitpoint->position - org), v_axis) / v_length;

		if (!circle)
		{
			if ( hitpoint->u < 0 || hitpoint->u > u_length ) return false;
			if ( hitpoint->v < 0 || hitpoint->v > v_length ) return false;
		}

		double dmy;
		double uu = hitpoint->u * (double)((UVPlane*)this)->material()->repeat;
		double vv = hitpoint->v * (double)((UVPlane*)this)->material()->repeat;

		if (uu < 0 || uu > 1) uu = modf(uu, &dmy);
		if (vv < 0 || vv > 1) vv = modf(vv, &dmy);
		if (uu < 0) uu = 1.0 - fabs(uu);
		if (vv < 0) vv = 1.0 - fabs(vv);

		if ( ((UVPlane*)this)->material()->texture )
		{
			int j = (((UVPlane*)this)->material()->texture->W()-1)*uu;
			int i = (((UVPlane*)this)->material()->texture->H()-1)*vv;

			hitpoint->material.color.x = (double)((UVPlane*)this)->material()->texture->cell(i,j).r/255.0;
			hitpoint->material.color.y = (double)((UVPlane*)this)->material()->texture->cell(i,j).g/255.0;
			hitpoint->material.color.z = (double)((UVPlane*)this)->material()->texture->cell(i,j).b/255.0;

			hitpoint->material.specular.x = (double)((UVPlane*)this)->material()->texture->cell(i,j).r/255.0;
			hitpoint->material.specular.y = (double)((UVPlane*)this)->material()->texture->cell(i,j).g/255.0;
			hitpoint->material.specular.z = (double)((UVPlane*)this)->material()->texture->cell(i,j).b/255.0;

			alp = (double)((UVPlane*)this)->material()->texture->cell(i,j).alp/255.0;
		}
		if ( ((UVPlane*)this)->material()->IBL() )
		{
			int j = (((UVPlane*)this)->material()->IBL_W()-1)*uu;
			int i = (((UVPlane*)this)->material()->IBL_H()-1)*vv;

			Material* ma = const_cast<Material*>(((UVPlane*)this)->material());
			hitpoint->material.color = ma->IBL_Color(i,j);
			hitpoint->material.emission = ((UVPlane*)this)->material()->ibl_texture_coef*ma->IBL_Color(i,j);
		}
		if ( ((UVPlane*)this)->material()->bump_texture )
		{
			int j = (((UVPlane*)this)->material()->bump_texture->W()-1)*uu;
			int i = (((UVPlane*)this)->material()->bump_texture->H()-1)*vv;

#if 0
			Color c((double)((UVPlane*)this)->material()->bump_texture->cell(i,j).r,
					(double)((UVPlane*)this)->material()->bump_texture->cell(i,j).g,
					(double)((UVPlane*)this)->material()->bump_texture->cell(i,j).b);
			c = c / 255.0;

			const double cc = c.length();
			hitpoint->material.color = hitpoint->material.color * cc;
#else

			//Tangent
			Vector3d sdir = u_axis;
			Vector3d tdir = v_axis;

			//compute them using central difference
			int u0 = j - 1;
			int u1 = j + 1;

			int v0 = i - 1;
			int v1 = i + 1;

#if 0
			if (u0 < 0) u0 = 0;
			if (v0 < 0) v0 = 0;
			if (u1 < 0) u1 = 0;
			if (v1 < 0) v1 = 0;
			if (u0 >= ((UVPlane*)this)->material()->bump_texture->W()) u0 = ((UVPlane*)this)->material()->bump_texture->W() - 1;
			if (u1 >= ((UVPlane*)this)->material()->bump_texture->W()) u1 = ((UVPlane*)this)->material()->bump_texture->W() - 1;
			if (v0 >= ((UVPlane*)this)->material()->bump_texture->H()) v0 = ((UVPlane*)this)->material()->bump_texture->H() - 1;
			if (v1 >= ((UVPlane*)this)->material()->bump_texture->H()) v1 = ((UVPlane*)this)->material()->bump_texture->H() - 1;
#else
			if (u0 < 0) u0 = ((UVPlane*)this)->material()->bump_texture->W() - 1;
			if (v0 < 0) v0 = ((UVPlane*)this)->material()->bump_texture->W() - 1;
			if (u1 < 0) u1 = ((UVPlane*)this)->material()->bump_texture->H() - 1;
			if (v1 < 0) v1 = ((UVPlane*)this)->material()->bump_texture->H() - 1;;
			if (u0 >= ((UVPlane*)this)->material()->bump_texture->W()) u0 = 0;
			if (u1 >= ((UVPlane*)this)->material()->bump_texture->W()) u1 = 0;
			if (v0 >= ((UVPlane*)this)->material()->bump_texture->H()) v0 = 0;
			if (v1 >= ((UVPlane*)this)->material()->bump_texture->H()) v1 = 0;
#endif

			Color px0(
				(double)((UVPlane*)this)->material()->bump_texture->cell(i, u0).r,
				(double)((UVPlane*)this)->material()->bump_texture->cell(i, u0).g,
				(double)((UVPlane*)this)->material()->bump_texture->cell(i, u0).b);

			Color px1(
				(double)((UVPlane*)this)->material()->bump_texture->cell(i, u1).r,
				(double)((UVPlane*)this)->material()->bump_texture->cell(i, u1).g,
				(double)((UVPlane*)this)->material()->bump_texture->cell(i, u1).b);

			Color py0(
				(double)((UVPlane*)this)->material()->bump_texture->cell(v0, j).r,
				(double)((UVPlane*)this)->material()->bump_texture->cell(v0, j).g,
				(double)((UVPlane*)this)->material()->bump_texture->cell(v0, j).b);


			Color py1(
				(double)((UVPlane*)this)->material()->bump_texture->cell(v1, j).r,
				(double)((UVPlane*)this)->material()->bump_texture->cell(v1, j).g,
				(double)((UVPlane*)this)->material()->bump_texture->cell(v1, j).b);

			//グレースケール（白黒画像）に変換
			Color YCrCb_Y(0.299, 0.587, 0.114);

			//central difference
			double x_gradient = dot((px1 - px0), YCrCb_Y)*0.5;
			double y_gradient = dot((py1 - py0), YCrCb_Y)*0.5;

			//法線をすでに反転している
			if (normal_vector_inverse < 0)
			{
				//一時的に元に戻す
				hitpoint->InversNormal();
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
		if ( alp == 0.0 )
		{
			return false;
		}
		if (alp < 1.0)
		{
			if (alp >= rnd->next01()) return false;
		}
		return true;
	}

	void MatrixTransformation(std::vector<Matrix4x4>& matrix)
	{
		base_plane->MatrixTransformation(matrix);
		for ( int j =0; j < matrix.size() ; j++ )
		{
			org = matrix[j] * org;
			p[0] = matrix[j] * p[0];
			p[1] = matrix[j] * p[1];

			Matrix4x4 g = matrix[j].vector_coordinate_transformation();
			u_axis = g*u_axis;
			v_axis = g*v_axis;
			normal = g*normal;
		}
	}

};
};


#endif
