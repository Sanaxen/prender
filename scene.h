#ifndef	_SCENE_H_
#define	_SCENE_H_

#include "constant.h"
#include "sphere.h"
#include "plane.h"
#include "circle.h"
#include "intersection.h"

#include <omp.h>

#include "polygon.h"
#include "sceneObject.h"
#include "scene_env.h"


namespace prender 
{
	// シーンを構成する物体等の全リスト
	extern ScenObjct* EntList;

	// シーンを構成する物体等の全リストの理解(生成)
	void CreateEntList(SceneEnv& env);

	// シーンとの交差判定関数
	bool intersect_scene_(const Ray &ray, Intersection *intersection, const int depth, int target_id=-1);

	bool intersect_scene(const Ray &ray, Intersection *intersection, const int depth, int target_id=-1);

	bool intersect_scene__(const Ray &ray, Intersection *intersection, const int depth, int target_id=-1);

};

#endif
