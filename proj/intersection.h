#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "Vector3d.h"
#include "constant.h"
#include "material.h"

namespace prender {

//��_�ʒu���
struct IntersectionPos 
{
	double distance;	//��_�܂ł̋���
	Vector3d normal;	//��_�ɂ�����@���x�N�g��
	Vector3d position;	//��_���W
	double u;		//V���W
	double v;		//V���W
	Material material;

	char bump;
	Vector3d bump_new_normal;	//��_�ɂ�����@���x�N�g��

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

//��_���
struct Intersection {
	IntersectionPos hitpoint;	//��_�ʒu���
	int object_id;				//��_�̑Ώۗv�fID

	Intersection() : object_id(-1) {}
};

};

#endif
