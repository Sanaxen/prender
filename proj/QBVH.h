#ifndef _QBCH_H
#define _QBCH_H

#include <vector>
#include <memory>
#include "BVH.h"

#include "BoundingBox.h"

//32bit => 0x80000000
#define QBVH_TOP_BIT	(1L << ((sizeof(long) * 8) - 1))

namespace prender 
{
class Entity;
class BVH;
class BVH_structure;

#ifdef USE_STXXL
typedef stxxl::VECTOR_GENERATOR<Entity*>::result			EntityList_t;
#else
typedef std::vector<Entity*>								EntityList_t;
#endif

class QBVH 
{
public:
	explicit QBVH() : m_root(NULL), m_allocatedQBVHNodeSize(0), m_usedNodeCount(0) { m_leafObjectArray.clear(); }
	~QBVH();

	void Construct(const EntityList_t &targets);
	bool CheckIntersection(const Ray &ray, Intersection &info) const;

	void CollectBoundingBoxes(int depth, std::vector<BoundingBox> &result); // for Visualization

private:
	void Construct_internal(size_t nextindex, const BVH &bvh, const BVH_structure *nextTarget);
	void MakeLeaf_internal(size_t index, const BVH_structure *leaf, size_t childindex);

	bool CheckIntersection_Leaf(const Ray &ray, size_t leafStartIndex, Intersection &hitResultDetail, Entity** hit_obj) const;

	static size_t SetChildindexAsLeaf(size_t childindex);
	static size_t GetInvalidChildIndex();
	static bool IsChildindexLeaf(size_t childindex);
	static bool IsValidIndex(size_t index);
	static size_t GetIndexOfObjectInChildLeaf(size_t childleafindex);

	void CollectBoundingBoxes_internal(int currentDepth, int targetDepth, int index, std::vector<BoundingBox> &result);

	void ReallocateQBVH_root(size_t addSize);

private:
	struct QBVH_structure  
	{
		__m128 bboxes[2][3];//4 float min-max xyz
		size_t children[4]; //4 children
		int axis_top;       //top axis
		int axis_left;      //left axis
		int axis_right;     //right axis
		int reserved;       //padding 
	};
	std::shared_ptr<QBVH_structure> m_root;
	size_t m_allocatedQBVHNodeSize, m_usedNodeCount;
	EntityList_t m_leafObjectArray;
};
}

#endif
