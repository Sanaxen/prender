#ifndef _ENTITY_H
#define _ENTITY_H


#include "Vector3d.h"
#include "matrix.h"
#include "ray.h"
#include "material.h"
#include "constant.h"
#include "intersection.h"
#include "BoundingBox.h"
#include "BVH.h"
#include "QBVH.h"

namespace prender {


enum EntityType {
	ENTITY_TYPE_INF_LIGHT,		// 
	ENTITY_TYPE_SPHERE,		// 球面
	ENTITY_TYPE_PLANE,		// 無限平面
	ENTITY_TYPE_UVPLANE,	// 有限平面
	ENTITY_TYPE_CIRCLE,		// 円盤
	ENTITY_TYPE_POLYGON		// ポリゴン
};

class BVH;
class QBVH;
class Entity
{
public:
	EntityType type;			//形状タイプ
	int id;						//形状のID
	Matrix4x4 mat;				//形状の変換マトリクス
	int normal_vector_inverse;	//法線を逆にする(-1)しない(1)

	int	back;					//裏を反射面にする(=1)しない(0)
	BoundingBox boundingBox;
	BVH* bvh;
	QBVH* qbvh;

	int material_id;
	MaterialList* materialList;

	bool light;			// 光源か？
	int light_id;		// 光源なら光源のID(LightList.listのインデックス)

	double area;		// 面積

	int scattering;		//

	unsigned char media_ignore_front_back;	//表裏の無視(PaticipatingMedia)

	Random* rnd;

	Entity()
	{
		light = false;
		light_id = -1;
		back = 1;
		normal_vector_inverse = 1;
		bvh = 0;
		qbvh = 0;
		area = PS_INF;
		material_id = -1;
		materialList = 0;
		rnd = NULL;
		media_ignore_front_back = '\0';
	}

	virtual ~Entity()
	{
		delete bvh;
		delete qbvh;
		bvh = 0;
		qbvh = 0;
		if ( material()->texture ) free( material()->texture);
		 material()->texture = 0;
	}

	inline Material* material()
	{
		return materialList->getMaterial(material_id);
	}

	virtual void MatrixTransformation(std::vector<Matrix4x4>& matrix) = 0;

	virtual	Entity* ConstructBVH() = 0;
	virtual	Entity* ConstructQBVH() = 0;

	virtual bool intersect(const Ray &ray, IntersectionPos *hitpoint) const = 0;
	virtual void CreateBoundingBox() = 0;

	virtual void CalcArea() = 0;

	virtual Vector3d randomPoint(Vector3d* nrm = 0, int* face_index_p = 0) = 0;

	inline bool isInsideDir(const Ray& ray, const IntersectionPos &hitpoint) const
	{
		//if ( dot(hitpoint.normal , ray.dir) < 0.0 ) return true;
		//return false;

		// Rayは物体の中に入って行くか？

		// 交差位置の法線（物体からのレイの入出を考慮）
		const Vector3d orienting_normal = dot(hitpoint.normal , ray.dir) < 0.0 ? hitpoint.normal: (-1.0 * hitpoint.normal); 
		if ( dot(orienting_normal, hitpoint.normal) > 0 )
		{
			//物体の中に入って行く
			return true;
		}
		//物体の外側に出て行く
		return false;
	}
	inline bool isInsideDir2(const Ray& ray, const IntersectionPos &hitpoint) const
	{
		bool s = isInsideDir(ray, hitpoint);
		if ( !media_ignore_front_back ) return s;

		switch (hitpoint.material.reflection_type)
		{
		case REFLECTION_TYPE_REFRACTION:
		case REFLECTION_TYPE_SSS_REFRACTION:
		case REFLECTION_TYPE_SSS_DIFFUSE:
		case REFLECTION_TYPE_SSS_WARD_BRDF:
		case REFLECTION_TYPE_REFRACTION_FRESNEL:
			return s;
		}
		return true;
	}
};
};

#endif