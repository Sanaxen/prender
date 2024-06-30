#ifndef _SAMPLER_H_
#define _SAMPLER_H_

#include "random.h"
#include "constant.h"
#include "vector3d.h"

namespace prender {

inline void sampleSphere(Random* rnd, const Vector3d& raydir, Vector3d& dir)
{
	const double r1 = PS_TWOPI * rnd->next01();
	const double r2 = 1.0 - 2.0 * rnd->next01();
	const double r22 = (r2 < 1.0) ? 1.0 - r2*r2 : 0.0;
	const double sqr22 = sqrt(r22);
	dir = normalize(Vector3d(sqr22* cos(r1), sqr22 * sin(r1), r2));
	
	//Vector3d w, u, v;
	//OrthonormalBasis(raydir, w, u, v);
	//dir =  normalize(u*sqrt(1.0 - r2*r2) * cos(r1) + v*sqrt(1.0 - r2*r2) * sin(r1) + w*r2);
}

//Henyey-Greenstein
inline void sampleHG(Random* rnd, const Vector3d& raydir, const double g, Vector3d& dir)
{
	if (g == 0.0) return sampleSphere(rnd, raydir, dir);
	
	//const double r1 = rnd->next01();
	//const double r2 = rnd->next01();
	//const double s = 1.0 - 2.0*r1;
	//const double cost = (s + 2.0*g*g*g * (-1.0 + r1) * r1 + g*g*s + 2.0*g*(1.0 - r1 + r1*r1)) / ((1.0 + g*s)*(1.0 + g*s));
	//const double sint = sqrt(1.0 - cost*cost);
	//dir = normalize(Vector3d(cos(2.0 * PS_PI * r2) * sint, sin(2.0 * PS_PI * r2) * sint, cost));
	

	//const double r1 = 2 * PS_PI * rnd->next01();
	//const double r0 = rnd->next01();
	//const double r2 = (1.0 + g*g) / 2.0 - pow((1.0 - g*g) / (1 + g - 2.0*g*r0), 2.0) / (2.0*g);
	//dir =  normalize(Vector3d(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2));

	//const double r1 = rnd->next01();
	//const double r2 = rnd->next01();
	//const double cost = (1.0 / 2.0*g)*(1.0 + g*g - pow((1.0 - g*g) / (1.0 - g + 2.0*g*r1), 2.0));
	//const double sint = sqrt(1.0 - cost*cost);
	//const double ph = 2 * PS_PI*r2;

	//// (sin(th)*cos(ph), sin(th)*sin(ph), cos(th))
	//dir = normalize(Vector3d(sint*cos(ph), sint*sin(ph), cost));

	Vector3d w, u, v;
	OrthonormalBasis(raydir, w, u, v);

	const double r1 = PS_TWOPI * rnd->next01();
	const double r2 = rnd->next01();
	const double s = (1.0 - g*g) / (1.0 - g + 2.0*g*r2);
	const double cost = (1.0 / 2.0*g)*(1.0 + g*g - s*s);
	const double sint = sqrt( std::max(0.0, 1.0 - cost*cost));
	dir = normalize( SphericalDirection(sint, cost, r1, u, v, w) );
	//dir = normalize(u*sint * cos(r1) + v*sint * sin(r1) + w*cost);

	//const double r1 = 2 * PS_PI * rnd->next01();
	//const double r2 = (1.0 / 2.0*g)*(1.0 + g*g - pow((1.0 - g*g) / (1.0 - g + 2.0*g*r1), 2.0));
	//dir = normalize(u*sqrt(1.0 - r2*r2) * cos(r1) + v*sqrt(1.0 - r2*r2) * sin(r1) + w*r2);

	//const double r1 = rnd->next01();
	//const double r2 = rnd->next01();
	//const double cost = (1.0 / (2.0*g))*(1.0 + g*g - pow((1.0 - g*g) / (1.0 - g + 2.0*g*r1), 2.0));
	//const double sint = sqrt(1.0 - cost*cost);
	//const double ph = 2.0 * PS_PI*r2;
	//const double zz = 1.0-raydir.z*raydir.z;

	//const double cosp = cos(ph);
	//const double sinp = sin(ph);

	//if ( zz > 0.000001 )
	//{
	//	const double ww = sqrt(zz);
	//	dir.x = sint*(raydir.x*raydir.z*cosp - raydir.y*sinp)/ww + raydir.x*cost;
	//	dir.y = sint*(raydir.y*raydir.z*cosp + raydir.x*sinp)/ww + raydir.y*cost;
	//	dir.z = -ww*sint*cosp + raydir.z*cost;
	//	dir = normalize( dir );
	//}else
	//{
	//	// (sin(th)*cos(ph), sin(th)*sin(ph), cos(th))
	//	if ( raydir.z > 0.0 )
	//	{
	//		dir = normalize(Vector3d(sint*cosp, sint*sinp, cost));
	//	}else
	//	{
	//		dir = normalize(Vector3d(sint*cosp, -sint*sinp, -cost));
	//	}
	//}
}


};

#endif
