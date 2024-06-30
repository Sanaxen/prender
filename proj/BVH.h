#ifndef _BVH_
#define _BVH_


#include <vector>
#include "entity.h"
#include "intersection.h"

namespace prender {

class Entity;
#ifdef USE_STXXL
typedef stxxl::VECTOR_GENERATOR<Entity*>::result			EntityList_t;
#else
typedef std::vector<Entity*>								EntityList_t;
#endif

class BVH_structure 
{
public:
	//BoundingBox box;
	double box[2][3];
	int children[2];
	int axis;
	std::vector<Entity *> objects;

	BVH_structure(){}
};

class BVH 
{
public:
	enum CONSTRUCTION_TYPE 
	{
		CONSTRUCTION_OBJECT_MEDIAN,
		CONSTRUCTION_OBJECT_SAH,
	};

public:
	explicit BVH() : m_root() {}
	~BVH();

	void Construct(const CONSTRUCTION_TYPE type, const EntityList_t &targets);
	bool CheckIntersection(const Ray &ray, Intersection &info) const;

	void CollectBoundingBoxes(int depth, std::vector<BoundingBox> &result); // for Visualization

	const BVH_structure *GetRootNode() const;
	size_t GetBVHNodeCount() const;
	const BVH_structure *GetFirstChild(const BVH_structure *parent) const;
	const BVH_structure *GetSecondChild(const BVH_structure *parent) const;
	bool IsLeaf(const BVH_structure *node) const;

private:
	void Construct_internal(const CONSTRUCTION_TYPE type, const EntityList_t &targets, int index);
	void MakeLeaf_internal(const EntityList_t &targets, int index);

	void CollectBoundingBoxes_internal(int currentDepth, int targetDepth, int index, std::vector<BoundingBox> &result);
	std::vector<BVH_structure> m_root;
};

}

#endif