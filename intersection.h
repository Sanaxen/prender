#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "Vector3d.h"
#include "constant.h"
#include "material.h"

namespace prender {

//交点位置情報
struct IntersectionPos 
{
	double distance;	//交点までの距離
	Vector3d normal;	//交点における法線ベクトル
	Vector3d position;	//交点座標
	double u;		//V座標
	double v;		//V座標
	Material material;

	char bump;
	Vector3d bump_new_normal;	//交点における法線ベクトル

	IntersectionPos() : distance(PS_INF), normal(), position(), u(-PS_INF), v(-PS_INF), bump_new_normal(), bump('\0'){}
	inline void InversNormal()
	{
		bump_new_normal = bump_new_normal*-1.0;
		normal = normal*-1.0;
	}

	inline Vector3d getNormal() const
	{
		if (bump) return bump_new_normal;
		
		return normal;
	}
};

//交点情報
struct Intersection {
	IntersectionPos hitpoint;	//交点位置情報
	int object_id;				//交点の対象要素ID

	Intersection() : object_id(-1) {}
};

};

#endif
