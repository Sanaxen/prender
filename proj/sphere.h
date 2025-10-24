#ifndef	_SPHERE_H_
#define	_SPHERE_H_

#include <cmath>
#include "entity.h"

namespace prender {

class Sphere: public Entity {
public:
	double radius;
	Vector3d position;
	int hemisphere;	//半球


	Sphere(const double radius, const Vector3d &position, int m) :
	  radius(radius), position(position)
	 {
		 type = ENTITY_TYPE_SPHERE;
		 material_id = m;
	}

	void CalcArea()
	{
		area = (PS_FORPI * radius*radius);
		if (hemisphere) area = 0.5*area;
	}

	Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0)
	{
		const double r1 = PS_TWOPI * rnd->next01();
		const double r2 = 1.0 - 2.0 * rnd->next01();
		const double sqrtr22 = sqrt(1.0 - r2*r2);
		if (nrm)
		{

			if (normal_vector_inverse < 0)
			{

				*nrm = -1.0*Vector3d(sqrtr22 * cos(r1), sqrtr22 * sin(r1), r2);
			}
			else
			{
				*nrm = Vector3d(sqrtr22 * cos(r1), sqrtr22 * sin(r1), r2);
			}
		}
		return position + (radius * Vector3d(sqrtr22 * cos(r1), sqrtr22 * sin(r1), r2));

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
		Vector3d ext(radius,radius,radius);
		Vector3d eps(0.01, 0.01, 0.01);
		boundingBox = BoundingBox(position-ext-eps, position+ext+eps);
	}


	// 入力のrayに対する交差点までの距離を返す。交差しなかったら0を返す。
	// rayとの交差判定を行う。交差したらtrue,さもなくばfalseを返す。
	bool intersect(const Ray &ray, IntersectionPos *hitpoint) const 
	{
		const Vector3d p_o = position - ray.org;
		const double b = dot(p_o, ray.dir);
		const double D4 = b * b - dot(p_o, p_o) + radius * radius;

		if (D4 < 0.0)
			return false;
		
		const double sqrt_D4 = sqrt(D4);
		const double t1 = b - sqrt_D4, t2 = b + sqrt_D4;
	
		if (t1 < PS_EPS && t2 < PS_EPS)
			return false;

		if (t1 > PS_EPS) {
			hitpoint->distance = t1;
		} else {
			hitpoint->distance = t2;
		}

		hitpoint->position = ray.org + hitpoint->distance * ray.dir;
		hitpoint->normal = normalize(hitpoint->position - position);
	 
		if (normal_vector_inverse < 0)
		{
			hitpoint->InversNormal();
		}
		hitpoint->material = *((Sphere*)this)->material();

		Vector3d p = hitpoint->position - position;

		if ( hemisphere )
		{
			if ( p.y < 0.0 ) return false;
		}
		const double pp[3] = {p.x, p.y, p.z};
		int X = 0;
		int Y = 1;
		int Z = 2;
		if (  hitpoint->material.hemisphere_map )
		{
			X = 2;
			Y = 0;
			Z = 1;
		}
		const Vector3d vv(pp[X]/radius, pp[Y]/radius, pp[Z]/radius);

		double alp = 1.0;	//αチャンネル
		if ( ((Sphere*)this)->material()->texture )
		{
			double th, phi;

			if ( fabs(pp[Z] - radius ) > PS_EPS16 )
			{
				//th = acos( pp[Z]/radius );	// 0～π
				//phi = acos(pp[X]/sqrt( pp[X]*pp[X] + pp[Y]*pp[Y]));	// 0～2π
				th = SphericalTheta(vv);
				phi = SphericalPhi(vv);

				hitpoint->u = phi / PS_TWOPI;
				hitpoint->v = th / PS_PI;

				if (hitpoint->material.equirectangular_map)
				{
					double len = std::sqrt(vv.x * vv.x + vv.y * vv.y + vv.z * vv.z);
					if (len == 0.0)
					{
						hitpoint->u = 0.0;
						hitpoint->v = 0.5;
					}
					else
					{
						double x = vv.x / len;
						double y = vv.y / len;
						double z = vv.z / len;

						double lambda = std::atan2(z, x);          // [-pi, pi]
						double phi = std::asin(std::fmax(-1.0, std::fmin(1.0, y))); // [-pi/2, pi/2]

						double u = (lambda + PS_PI) / (2.0 * PS_PI); // [0,1)
						double v = (0.5 * PS_PI - phi) / PS_PI;      // [0,1]

						// wrap/clamp
						u = u - std::floor(u);
						if (v < 0.0) v = 0.0; else if (v > 1.0) v = 1.0;

						hitpoint->u = u;
						hitpoint->v = v;
					}
				}
				if (  hitpoint->material.angular_map )
				{
					double r = (1.0 / PS_PI) * acos(vv.z);
					r /= sqrt(vv.x * vv.x + vv.y * vv.y);
					hitpoint->u = vv.x * r;  // u are in range [-1.0, 1.0]
					hitpoint->v = vv.y * r;  // v are in range [-1.0, 1.0]
					hitpoint->u = 0.5 * (hitpoint->u) + 0.5;
					hitpoint->v = 0.5 * (hitpoint->v) + 0.5;
				}
				if (  hitpoint->material.panoramic_map )
				{
					double s = vv.y;
					if ( vv.y < -1.0 ) s = -1.0;
					if ( vv.y >  1.0 ) s =  1.0;
					hitpoint->u = 1.0 + atan2(vv.x, -vv.z)/PS_PI;  // u are in range [0. 2]
					hitpoint->v = acos(s)/PS_PI;				// v are in range [0, 1.0]
					double a = hitpoint->u;
					double b = hitpoint->v;
					hitpoint->u = 0.5*a;
					hitpoint->v = 1.0 - b;
				}

				if ( hemisphere ) hitpoint->u *= 2.0;

				double dmy;
				double uu = hitpoint->u * (double)((Sphere*)this)->material()->repeat;
				double vv = hitpoint->v * (double)((Sphere*)this)->material()->repeat;

				if (uu < 0 || uu > 1) uu = modf(uu, &dmy);
				if (vv < 0 || vv > 1) vv = modf(vv, &dmy);
				if (uu < 0) uu = 1.0 - fabs(uu);
				if (vv < 0) vv = 1.0 - fabs(vv);

				int j = (((Sphere*)this)->material()->texture->W()-1)*uu;
				int i = (((Sphere*)this)->material()->texture->H()-1)*vv;

				if ( hitpoint->material.hemisphere_map )
				{
					j = (((Sphere*)this)->material()->texture->W()-1)/2 + (((Sphere*)this)->material()->texture->W()-1)*0.5*pp[X]/radius;
					i = (((Sphere*)this)->material()->texture->H()-1)/2 + (((Sphere*)this)->material()->texture->H()-1)*0.5*pp[Y]/radius;
				}

				if (i < 0) i = 0;
				if (i >= ((Sphere*)this)->material()->texture->H() - 1) i = ((Sphere*)this)->material()->texture->H() - 1;
				if (j < 0) j = 0;
				if (j >= ((Sphere*)this)->material()->texture->W() - 1) j = ((Sphere*)this)->material()->texture->W() - 1;

				hitpoint->material.color.x = (double)((Sphere*)this)->material()->texture->cell(i,j).r/255.0;
				hitpoint->material.color.y = (double)((Sphere*)this)->material()->texture->cell(i,j).g/255.0;
				hitpoint->material.color.z = (double)((Sphere*)this)->material()->texture->cell(i,j).b/255.0;

				hitpoint->material.specular.x = (double)((Sphere*)this)->material()->texture->cell(i,j).r/255.0;
				hitpoint->material.specular.y = (double)((Sphere*)this)->material()->texture->cell(i,j).g/255.0;
				hitpoint->material.specular.z = (double)((Sphere*)this)->material()->texture->cell(i,j).b/255.0;

				alp = (double)((Sphere*)this)->material()->texture->cell(i,j).alp/255.0;
			}
		}

		if ( ((Sphere*)this)->material()->IBL() )
		{
			double th, phi;

			if ( fabs(pp[Z] - radius ) > PS_EPS16 )
			{
				//th = acos( pp[Z]/radius );	// 0～π
				//phi = acos(pp[X]/sqrt( pp[X]*pp[X] + pp[Y]*pp[Y]));	// 0～2π
				th = SphericalTheta(vv);
				phi = SphericalPhi(vv);

				hitpoint->u = th / PS_PI;
				hitpoint->v = phi / PS_TWOPI;

				if (hitpoint->material.equirectangular_map)
				{
					double len = std::sqrt(vv.x * vv.x + vv.y * vv.y + vv.z * vv.z);
					if (len == 0.0)
					{
						hitpoint->u = 0.0;
						hitpoint->v = 0.5;
					}
					else
					{
						double x = vv.x / len;
						double y = vv.y / len;
						double z = vv.z / len;

						double lambda = std::atan2(z, x);          // [-pi, pi]
						double phi = std::asin(std::fmax(-1.0, std::fmin(1.0, y))); // [-pi/2, pi/2]

						double u = (lambda + PS_PI) / (2.0 * PS_PI); // [0,1)
						double v = (0.5 * PS_PI - phi) / PS_PI;      // [0,1]

						// wrap/clamp
						u = u - std::floor(u);
						if (v < 0.0) v = 0.0; else if (v > 1.0) v = 1.0;

						hitpoint->u = u;
						hitpoint->v = v;
					}
				}
				// http://www.pauldebevec.com/Probes/
				// Thus, if we consider the images to be normalized to have coordinates u=[-1,1], v=[-1,1], 
				// we have theta=atan2(v,u), phi=pi*sqrt(u*u+v*v). 
				// The unit vector pointing in the corresponding direction is obtained by rotating (0,0,-1) 
				// by phi degrees around the y (up) axis and then theta degrees around the -z (forward) axis.
				// If for a direction vector in the world (Dx, Dy, Dz), 
				// the corresponding (u,v) coordinate in the light probe image is (Dx*r,Dy*r) 
				// where r=(1/pi)*acos(Dz)/sqrt(Dx^2 + Dy^2).
				if (  hitpoint->material.angular_map )
				{
					double s = vv.z;
					if ( vv.z < -1.0 ) s = -1.0;
					if ( vv.z >  1.0 ) s =  1.0;
					double r = (1.0 / PS_PI) * acos(s);
					r /= sqrt(vv.x * vv.x + vv.y * vv.y);
					hitpoint->u = vv.x * r;  // u are in range [-1.0, 1.0]
					hitpoint->v = vv.y * r;  // v are in range [-1.0, 1.0]
					hitpoint->u = 0.5 * (hitpoint->u) + 0.5;
					hitpoint->v = 0.5 * (hitpoint->v) + 0.5;
				}

				//http://gl.ict.usc.edu/Data/HighResProbes/
				//Panoramic Format
				// if we consider the images to have a rectangular image domain of u=[0,2], v=[0,1], 
				// we have theta= pi*(u-1), phi=pi*v. 
				// The unit vector pointing in the corresponding direction is obtained by 
				// (Dx,Dy,Dz) = (sin(phi)*sin(theta), cos(phi), -sin(phi)*cos(theta)). 
				// For the reverse mapping from the direction vector in the world (Dx, Dy, Dz), 
				// the corresponding (u,v) coordinate in the light probe image is 
				// ( 1 + atan2(Dx,-Dz) / pi, arccos(Dy) / pi).
				if (  hitpoint->material.panoramic_map )
				{
					double s = vv.y;
					if ( vv.y < -1.0 ) s = -1.0;
					if ( vv.y >  1.0 ) s =  1.0;
					hitpoint->u = 1.0 + atan2(vv.x, -vv.z)/PS_PI;  // u are in range [0. 2]
					hitpoint->v = acos(s)/PS_PI;				// v are in range [0, 1.0]
					double a = hitpoint->u;
					double b = hitpoint->v;
					hitpoint->u = 0.5*a;
					hitpoint->v = 1.0 - b;
				}

				if ( hemisphere ) hitpoint->u *= 2.0;

				int j = (((Sphere*)this)->material()->IBL_W()-1)*hitpoint->u;
				int i = (((Sphere*)this)->material()->IBL_H()-1)*hitpoint->v;

				if ( hitpoint->material.hemisphere_map )
				{
					j = (((Sphere*)this)->material()->IBL_W()-1)/2 + (((Sphere*)this)->material()->IBL_W()-1)*0.5*pp[X]/radius;
					i = (((Sphere*)this)->material()->IBL_H()-1)/2 + (((Sphere*)this)->material()->IBL_H()-1)*0.5*pp[Y]/radius;
				}
				if (i < 0) i = 0;
				if (i >= ((Sphere*)this)->material()->IBL_H() - 1) i = ((Sphere*)this)->material()->IBL_H() - 1;
				if (j < 0) j = 0;
				if (j >= ((Sphere*)this)->material()->IBL_W() - 1) j = ((Sphere*)this)->material()->IBL_W() - 1;


				Material* ma = const_cast<Material*>(((Sphere*)this)->material());
				hitpoint->material.color = ma->IBL_Color(i,j);
				hitpoint->material.emission = ((Sphere*)this)->material()->ibl_texture_coef*ma->IBL_Color(i,j);
			}
		}
		if ( ((Sphere*)this)->material()->bump_texture )
		{
			double th, phi;

			if ( fabs(p.z - radius ) > PS_EPS16 )
			{
				//th = acos( pp[Z]/radius );	// 0～π
				//phi = acos(pp[X]/sqrt( pp[X]*pp[X] + pp[Y]*pp[Y]));	// 0～2π
				th = SphericalTheta(vv);
				phi = SphericalPhi(vv);

				hitpoint->u = th/PS_PI;
				hitpoint->v = phi/PS_TWOPI;

				if (hitpoint->material.equirectangular_map)
				{
					double len = std::sqrt(vv.x * vv.x + vv.y * vv.y + vv.z * vv.z);
					if (len == 0.0)
					{
						hitpoint->u = 0.0;
						hitpoint->v = 0.5;
					}
					else
					{
						double x = vv.x / len;
						double y = vv.y / len;
						double z = vv.z / len;

						double lambda = std::atan2(z, x);          // [-pi, pi]
						double phi = std::asin(std::fmax(-1.0, std::fmin(1.0, y))); // [-pi/2, pi/2]

						double u = (lambda + PS_PI) / (2.0 * PS_PI); // [0,1)
						double v = (0.5 * PS_PI - phi) / PS_PI;      // [0,1]

						// wrap/clamp
						u = u - std::floor(u);
						if (v < 0.0) v = 0.0; else if (v > 1.0) v = 1.0;

						hitpoint->u = u;
						hitpoint->v = v;
					}
				}
				if (  hitpoint->material.angular_map )
				{
					double r = (1.0 / PS_PI) * acos(vv.z);
					r /= sqrt(vv.x * vv.x + vv.y * vv.y);
					hitpoint->u = vv.x * r;  // u are in range [-1.0, 1.0]
					hitpoint->v = vv.y * r;  // v are in range [-1.0, 1.0]
					hitpoint->u = 0.5 * (hitpoint->u) + 0.5;
					hitpoint->v = 0.5 * (hitpoint->v) + 0.5;
				}
				if (  hitpoint->material.panoramic_map )
				{
					double s = vv.y;
					if ( vv.y < -1.0 ) s = -1.0;
					if ( vv.y >  1.0 ) s =  1.0;
					hitpoint->u = 1.0 + atan2(vv.x, -vv.z)/PS_PI;  // u are in range [0. 2]
					hitpoint->v = acos(s)/PS_PI;				// v are in range [0, 1.0]
					double a = hitpoint->u;
					double b = hitpoint->v;
					hitpoint->u = 0.5*a;
					hitpoint->v = 1.0 - b;
				}
				if ( hemisphere ) hitpoint->u *= 2.0;

				double dmy;
				double uu = hitpoint->u * (double)((Sphere*)this)->material()->repeat;
				double vv = hitpoint->v * (double)((Sphere*)this)->material()->repeat;

				if (uu < 0 || uu > 1) uu = modf(uu, &dmy);
				if (vv < 0 || vv > 1) vv = modf(vv, &dmy);
				if (uu < 0) uu = 1.0 - fabs(uu);
				if (vv < 0) vv = 1.0 - fabs(vv);


				int j = (((Sphere*)this)->material()->bump_texture->W()-1)*uu;
				int i = (((Sphere*)this)->material()->bump_texture->H()-1)*vv;

				if ( hitpoint->material.hemisphere_map )
				{
					j = (((Sphere*)this)->material()->bump_texture->W()-1)/2 + (((Sphere*)this)->material()->bump_texture->W()-1)*0.5*pp[X]/radius;
					i = (((Sphere*)this)->material()->bump_texture->H()-1)/2 + (((Sphere*)this)->material()->bump_texture->H()-1)*0.5*pp[Y]/radius;
				}

#if 0
				Color c((double)((Sphere*)this)->material()->bump_texture->cell(i,j).r,
						(double)((Sphere*)this)->material()->bump_texture->cell(i,j).g,
						(double)((Sphere*)this)->material()->bump_texture->cell(i,j).b);
				c = c / 255.0;
				const double cc = c.length();
				hitpoint->material.color = hitpoint->material.color * cc;
#else

				//CalculateTangent
				Vector3d sdir = Vector3d( cos(th)*cos(phi), cos(th)*sin(phi), -sin(th));
				Vector3d tdir = Vector3d( -sin(th)*sin(phi), sin(th)*cos(phi), 0.0);

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
				if (u0 >= ((Sphere*)this)->material()->bump_texture->W()) u0 = ((Sphere*)this)->material()->bump_texture->W() - 1;
				if (u1 >= ((Sphere*)this)->material()->bump_texture->W()) u1 = ((Sphere*)this)->material()->bump_texture->W() - 1;
				if (v0 >= ((Sphere*)this)->material()->bump_texture->H()) v0 = ((Sphere*)this)->material()->bump_texture->H() - 1;
				if (v1 >= ((Sphere*)this)->material()->bump_texture->H()) v1 = ((Sphere*)this)->material()->bump_texture->H() - 1;
#else
				if (u0 < 0) u0 = ((Sphere*)this)->material()->bump_texture->W() - 1;
				if (v0 < 0) v0 = ((Sphere*)this)->material()->bump_texture->W() - 1;
				if (u1 < 0) u1 = ((Sphere*)this)->material()->bump_texture->H() - 1;
				if (v1 < 0) v1 = ((Sphere*)this)->material()->bump_texture->H() - 1;;
				if (u0 >= ((Sphere*)this)->material()->bump_texture->W()) u0 = 0;
				if (u1 >= ((Sphere*)this)->material()->bump_texture->W()) u1 = 0;
				if (v0 >= ((Sphere*)this)->material()->bump_texture->H()) v0 = 0;
				if (v1 >= ((Sphere*)this)->material()->bump_texture->H()) v1 = 0;
#endif


				Color px0(
					(double)((Sphere*)this)->material()->bump_texture->cell(i, u0).r,
					(double)((Sphere*)this)->material()->bump_texture->cell(i, u0).g,
					(double)((Sphere*)this)->material()->bump_texture->cell(i, u0).b);

				Color px1(
					(double)((Sphere*)this)->material()->bump_texture->cell(i, u1).r,
					(double)((Sphere*)this)->material()->bump_texture->cell(i, u1).g,
					(double)((Sphere*)this)->material()->bump_texture->cell(i, u1).b);

				Color py0(
					(double)((Sphere*)this)->material()->bump_texture->cell(v0, j).r,
					(double)((Sphere*)this)->material()->bump_texture->cell(v0, j).g,
					(double)((Sphere*)this)->material()->bump_texture->cell(v0, j).b);


				Color py1(
					(double)((Sphere*)this)->material()->bump_texture->cell(v1, j).r,
					(double)((Sphere*)this)->material()->bump_texture->cell(v1, j).g,
					(double)((Sphere*)this)->material()->bump_texture->cell(v1, j).b);

				//グレースケール（白黒画像）に変換
				Color YCrCb_Y(0.299, 0.587, 0.114);

				//central difference
				double x_gradient = 0.1*dot((px1 - px0), YCrCb_Y)*0.5;
				double y_gradient = 0.1*dot((py1 - py0), YCrCb_Y)*0.5;

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
		for ( int j = 0; j < matrix.size() ; j++ )
		{
			position = matrix[j] * position;
		}
	}

	inline Vector3d SphericalCoordinatesToXYZ( const double th, const double phi)
	{
		return Vector3d(
			radius*sin(th)*cos(phi),
			radius*sin(th)*sin(phi),
			radius*cos(th)) + position;
	}

	inline int XYZToSphericalCoordinates(const Vector3d& pos, double& th, double& phi)
	{
		Vector3d p = pos - position;

		if ( fabs(p.z - radius ) < PS_EPS16 )
		{
			return -1;
		}

		th = acos( p.z/radius );	// 0～π

		phi = p.x/sqrt( p.x*p.x + p.y*p.y);	// 0～2π

		return 1;
	}

	//物体が球ならpntはその球の中にあるはず
	inline int isIn(const Vector3d& pnt) const
	{
		Vector3d dd = pnt - position;
		double r = sqrt(dot(dd, dd));

		if( r > radius)
		{
			return -1;
		}
		if( r < radius)
		{
			return 1;
		}
		return 0;
	}

};

};

#endif
