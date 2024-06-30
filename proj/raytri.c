/* Ray-Triangle Intersection Test Routines          */
/* Different optimizations of my and Ben Trumbore's */
/* code from journals of graphics tools (JGT)       */
/* http://www.acm.org/jgt/                          */
/* by Tomas Moller, May 2000                        */

#include <math.h>
#include "constant.h"

//#define EPSILON 0.000001
#define EPSILON PS_EPS

#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 

/* the original jgt code */
int intersect_triangle(double orig[3], double dir[3],
		       double vert0[3], double vert1[3], double vert2[3],
		       double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tVector3d[3], pVector3d[3], qVector3d[3];
   double det,inv_det;

   /* find Vector3dtors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pVector3d, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pVector3d);

   if (det > -EPSILON && det < EPSILON)
     return 0;
   inv_det = 1.0 / det;

   /* calculate distance from vert0 to ray origin */
   SUB(tVector3d, orig, vert0);

   /* calculate U parameter and test bounds */
   *u = DOT(tVector3d, pVector3d) * inv_det;
   if (*u < 0.0 || *u > 1.0)
     return 0;

   /* prepare to test V parameter */
   CROSS(qVector3d, tVector3d, edge1);

   /* calculate V parameter and test bounds */
   *v = DOT(dir, qVector3d) * inv_det;
   if (*v < 0.0 || *u + *v > 1.0)
     return 0;

   /* calculate t, ray intersects triangle */
   *t = DOT(edge2, qVector3d) * inv_det;

   return 1;
}


/* code rewritten to do tests on the sign of the determinant */
/* the division is at the end in the code                    */
int intersect_triangle1(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tVector3d[3], pVector3d[3], qVector3d[3];
   double det,inv_det;

   /* find Vector3dtors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pVector3d, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pVector3d);

   if (det > EPSILON)
   {
      /* calculate distance from vert0 to ray origin */
      SUB(tVector3d, orig, vert0);
      
      /* calculate U parameter and test bounds */
      *u = DOT(tVector3d, pVector3d);
      if (*u < 0.0 || *u > det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qVector3d, tVector3d, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qVector3d);
      if (*v < 0.0 || *u + *v > det)
	 return 0;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate distance from vert0 to ray origin */
      SUB(tVector3d, orig, vert0);
      
      /* calculate U parameter and test bounds */
      *u = DOT(tVector3d, pVector3d);
/*      printf("*u=%f\n",(float)*u); */
/*      printf("det=%f\n",det); */
      if (*u > 0.0 || *u < det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qVector3d, tVector3d, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qVector3d) ;
      if (*v > 0.0 || *u + *v < det)
	 return 0;
   }
   else return 0;  /* ray is parallell to the plane of the triangle */


   inv_det = 1.0 / det;

   /* calculate t, ray intersects triangle */
   *t = DOT(edge2, qVector3d) * inv_det;
   (*u) *= inv_det;
   (*v) *= inv_det;

   return 1;
}

/* code rewritten to do tests on the sign of the determinant */
/* the division is before the test of the sign of the det    */
int intersect_triangle2(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tVector3d[3], pVector3d[3], qVector3d[3];
   double det,inv_det;

   /* find Vector3dtors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pVector3d, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pVector3d);

   /* calculate distance from vert0 to ray origin */
   SUB(tVector3d, orig, vert0);
   inv_det = 1.0 / det;
   
   if (det > EPSILON)
   {
      /* calculate U parameter and test bounds */
      *u = DOT(tVector3d, pVector3d);
      if (*u < 0.0 || *u > det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qVector3d, tVector3d, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qVector3d);
      if (*v < 0.0 || *u + *v > det)
	 return 0;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate U parameter and test bounds */
      *u = DOT(tVector3d, pVector3d);
      if (*u > 0.0 || *u < det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qVector3d, tVector3d, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qVector3d) ;
      if (*v > 0.0 || *u + *v < det)
	 return 0;
   }
   else return 0;  /* ray is parallell to the plane of the triangle */

   /* calculate t, ray intersects triangle */
   *t = DOT(edge2, qVector3d) * inv_det;
   (*u) *= inv_det;
   (*v) *= inv_det;

   return 1;
}

/* code rewritten to do tests on the sign of the determinant */
/* the division is before the test of the sign of the det    */
/* and one CROSS has been moved out from the if-else if-else */
int intersect_triangle3(const double orig[3], const double dir[3],
			const double vert0[3], const double vert1[3], const double vert2[3],
			double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tVector3d[3], pVector3d[3], qVector3d[3];
   double det,inv_det;

   /* find Vector3dtors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pVector3d, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pVector3d);

   /* calculate distance from vert0 to ray origin */
   SUB(tVector3d, orig, vert0);
   inv_det = 1.0 / det;
   
   CROSS(qVector3d, tVector3d, edge1);
      
   if (det > EPSILON)
   {
      *u = DOT(tVector3d, pVector3d);
      if (*u < 0.0 || *u > det)
	 return 0;
            
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qVector3d);
      if (*v < 0.0 || *u + *v > det)
	 return 0;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate U parameter and test bounds */
      *u = DOT(tVector3d, pVector3d);
      if (*u > 0.0 || *u < det)
	 return 0;
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qVector3d) ;
      if (*v > 0.0 || *u + *v < det)
	 return 0;
   }
   else return 0;  /* ray is parallell to the plane of the triangle */

   *t = DOT(edge2, qVector3d) * inv_det;
   (*u) *= inv_det;
   (*v) *= inv_det;

   return 1;
}
