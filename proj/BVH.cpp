#include <iostream>
#include "BVH.h"
#include <stack>
#include <algorithm>
#include <emmintrin.h>

using namespace std;

//#define assert(x)	if ( !(x) ) abort()
#define assert(x)	/*Å@*/


namespace prender {

BVH::~BVH()
{
  //delete [] m_root;
  m_root.clear();
}

bool BVH::CheckIntersection(const Ray &ray, Intersection &info) const 
{
	info.hitpoint.distance = PS_INF;
	info.object_id = -1;

	std::vector<size_t> next_list;
	next_list.reserve(m_root.size());
	next_list.push_back(0);

	const double rayDir[3] = {ray.dir.x, ray.dir.y, ray.dir.z};
	const double rayOrig[3] = {ray.org.x, ray.org.y, ray.org.z};

	while (!next_list.empty()) 
	{
		const BVH_structure *next = &m_root[next_list.back()];
		next_list.pop_back();

		if (next->children[0] == -1) 
		{
			// leaf
			bool isHit = false;
			int index = -1;
			for (size_t i=0; next->objects[i]; i++) 
			{
				Intersection hit;
				if (next->objects[i]->intersect(ray, &hit.hitpoint)) 
				{
					//printf("%d ", i);
					if (info.hitpoint.distance > hit.hitpoint.distance && hit.hitpoint.distance > PS_EPS) 
					{
						isHit = true;
						info.hitpoint = hit.hitpoint;
						//index = next->objects[i]->id;
						//printf("--><BVH>    <%d> %f %f %f  %f %f %f\n", index,
						//	info.hitpoint.position.x,
						//	info.hitpoint.position.y,
						//	info.hitpoint.position.z,
						//	info.hitpoint.normal.x,
						//	info.hitpoint.normal.y,
						//	info.hitpoint.normal.z);
						//fflush(stdout);

					}
				}
			}
			//printf("\n");
			if (isHit) 
			{
						//printf("<BVH>    <%d> %f %f %f  %f %f %f\n", index,
						//	info.hitpoint.position.x,
						//	info.hitpoint.position.y,
						//	info.hitpoint.position.z,
						//	info.hitpoint.normal.x,
						//	info.hitpoint.normal.y,
						//	info.hitpoint.normal.z);
						//fflush(stdout);

				return true;
			}
		} else {
			// internal node
			// check intersection with children
			double dist1 = PS_INF, dist2 = PS_INF;
			const BVH_structure *next1 = &m_root[next->children[0]];

			bool hit1 = BoundingBox::CheckIntersection(rayDir, rayOrig, next1->box[0], next1->box[1], dist1);//next1->box.Intersect(ray, dist1);
			bool hit2 = false;

			if (next->children[1] >= 0) 
			{
				const BVH_structure *next2 = &m_root[next->children[1]];
				hit2 = BoundingBox::CheckIntersection(rayDir, rayOrig, next2->box[0], next2->box[1], dist2); //m_root[next->children[1]].box.Intersect(ray, dist2);
			}
			if (hit1 && hit2) 
			{
				if (dist1 < dist2) {
					// check child1 at first
					next_list.push_back(next->children[1]);
					next_list.push_back(next->children[0]);
				} else 
				{
					// check child2 at first
					next_list.push_back(next->children[0]);
					next_list.push_back(next->children[1]);
				}
			} else if (hit1) 
			{
				next_list.push_back(next->children[0]);
			} else if (hit2) 
			{
				next_list.push_back(next->children[1]);
			}
		}
	}

	return !(info.object_id == -1);
}

void BVH::Construct(const BVH::CONSTRUCTION_TYPE type, const EntityList_t &targets)
{
	m_root.clear();

	// how many nodes are required????
	// Should I use dynamic array (vector)???

	//m_root = new BVH_structure[8*targets.size()+1]; // max size is 2*N+1
	m_root.reserve(8*targets.size());
	//m_bvh_node_size = targets.size()*targets.size()+1;
	//m_bvh_node_size = 8*targets.size()+1;
	m_root.push_back(BVH_structure());

	Construct_internal(type, targets, 0);
}

namespace 
{
	void CalcBoundingBoxOfObjects(const EntityList_t &objects, BoundingBox &boxResult)
{
  Vector3d box_min(objects[0]->boundingBox.min());
  Vector3d box_max(objects[0]->boundingBox.max());

  const int sz = objects.size();
  for ( int i = 1; i < sz ; i++ )
  {
	  const Entity* a = objects[i];
      if (a->boundingBox.min().x < box_min.x) box_min.x = a->boundingBox.min().x;
      if (a->boundingBox.min().y < box_min.y) box_min.y = a->boundingBox.min().y;
      if (a->boundingBox.min().z < box_min.z) box_min.z = a->boundingBox.min().z;
      if (a->boundingBox.max().x > box_max.x) box_max.x = a->boundingBox.max().x;
      if (a->boundingBox.max().y > box_max.y) box_max.y = a->boundingBox.max().y;
      if (a->boundingBox.max().z > box_max.z) box_max.z = a->boundingBox.max().z;
  }
  boxResult.SetBox(box_min, box_max);
}

template <typename FLOATING>
void CalcBoundingBoxOfObjects(const EntityList_t &objects, FLOATING min[3], FLOATING max[3])
{
  min[0] = static_cast<FLOATING>(objects[0]->boundingBox.min().x);
  min[1] = static_cast<FLOATING>(objects[0]->boundingBox.min().y);
  min[2] = static_cast<FLOATING>(objects[0]->boundingBox.min().z);
  max[0] = static_cast<FLOATING>(objects[0]->boundingBox.max().x);
  max[1] = static_cast<FLOATING>(objects[0]->boundingBox.max().y);
  max[2] = static_cast<FLOATING>(objects[0]->boundingBox.max().z);

  const int sz = objects.size();
  for ( int i = 1; i < sz; i++ )
  {
	  const Entity* a = objects[i];
      if (a->boundingBox.min().x < min[0]) min[0] = static_cast<FLOATING>(a->boundingBox.min().x);
      if (a->boundingBox.min().y < min[1]) min[1] = static_cast<FLOATING>(a->boundingBox.min().y);
      if (a->boundingBox.min().z < min[2]) min[2] = static_cast<FLOATING>(a->boundingBox.min().z);
      if (a->boundingBox.max().x > max[0]) max[0] = static_cast<FLOATING>(a->boundingBox.max().x);
      if (a->boundingBox.max().y > max[1]) max[1] = static_cast<FLOATING>(a->boundingBox.max().y);
      if (a->boundingBox.max().z > max[2]) max[2] = static_cast<FLOATING>(a->boundingBox.max().z);
  }

}
}

void BVH::MakeLeaf_internal(const EntityList_t &targets, int index)
{
	BVH_structure *st = &m_root[index];
	st->children[0] = st->children[1] = -1;

	const int sz = targets.size();
	for (size_t i=0; i<sz; i++)	st->objects.push_back(targets[i]);

	st->objects.push_back(NULL);
	CalcBoundingBoxOfObjects<double>(targets, st->box[0], st->box[1]);
}


void BVH::Construct_internal(const CONSTRUCTION_TYPE type, const EntityList_t &targets, int index)
{
	//assert (index < m_bvh_node_size);

	const double T_aabb = 1.0; // cost of check intersection of AABB
	const double T_tri = 1.0; // cost of check intersection of Triangle

	// calculate this node's bounding box
	CalcBoundingBoxOfObjects(targets, m_root[index].box[0], m_root[index].box[1]);
	const double currentBoxSurface = BoundingBox::CalcSurfaceArea(m_root[index].box[0], m_root[index].box[1]);
	const double currentBoxSurfaceInverse = 1.0/currentBoxSurface;

	EntityList_t axisSortedLeft[3], axisSortedRight[3];

	int bestAxis = -1;
	int bestIndex = -1;
	double bestCost = -1;

	for (int axis=0; axis<3; axis++) 
	{

		axisSortedLeft[axis] = targets;

		// sort objects by the axis
		std::sort(axisSortedLeft[axis].begin(), axisSortedLeft[axis].end(),
			[&axis](Entity * const &a, Entity * const &b) -> bool {
			switch (axis)
			{
			case 0: return a->boundingBox.position().x < b->boundingBox.position().x;
			case 1: return a->boundingBox.position().y < b->boundingBox.position().y;
			case 2: return a->boundingBox.position().z < b->boundingBox.position().z;
			}
			assert(false);
			return a->boundingBox.position().x < b->boundingBox.position().x;
		});

		switch (type) {
		case CONSTRUCTION_OBJECT_MEDIAN:
			if (targets.size() > 2)
			{
				// select the media of the objects
				const int select = axisSortedLeft[axis].size()/2-1;
				const Entity *selectedObj = axisSortedLeft[axis][select];
				const Entity *leftest = axisSortedLeft[axis][0];
				const Entity *rightest = axisSortedLeft[axis][axisSortedLeft[axis].size()-1];
				const double dist = (selectedObj->boundingBox.position()-leftest->boundingBox.position()).sqr() + 
				(selectedObj->boundingBox.position()-rightest->boundingBox.position()).sqr();

				if (dist > bestCost) 
				{
					bestCost = dist;
					bestIndex = select;
					bestAxis = axis;
				}
			}
		break;

		case CONSTRUCTION_OBJECT_SAH: 
		{
			double *leftAreas = new double[axisSortedLeft[axis].size()];
			double *rightAreas = new double[axisSortedLeft[axis].size()];
			rightAreas[axisSortedLeft[axis].size()-1] = PS_INF;

			BoundingBox boxTmp;
			axisSortedRight[axis].push_back(axisSortedLeft[axis].back());
			axisSortedLeft[axis].pop_back();

			// calculate right area
			CalcBoundingBoxOfObjects(axisSortedRight[axis], boxTmp);
			const int sz = axisSortedLeft[axis].size();
			for (int i=sz-1; i>=0; i--) 
			{
				rightAreas[i] = boxTmp.CalcSurfaceArea();
				const Entity *obj = axisSortedLeft[axis][i];
				boxTmp.MergeAnotherBox(obj->boundingBox);

				axisSortedRight[axis].push_back(axisSortedLeft[axis].back());
				axisSortedLeft[axis].pop_back();
			}

			axisSortedLeft[axis].push_back(axisSortedRight[axis].back());
			axisSortedRight[axis].pop_back();
			CalcBoundingBoxOfObjects(axisSortedLeft[axis], boxTmp);
      
			const double leafCost = T_tri * targets.size();
			const int sz2 = axisSortedRight[axis].size();
			for (int i=sz2-1, cutIndex=0; i>=0; i--, cutIndex++) 
			{
				// calculate both surface area
				leftAreas[cutIndex] = boxTmp.CalcSurfaceArea();
				const Entity *nextObj = axisSortedRight[axis].back();
				boxTmp.MergeAnotherBox(nextObj->boundingBox);

				// move right to left
				axisSortedLeft[axis].push_back(axisSortedRight[axis].back());
				axisSortedRight[axis].pop_back();

				// calc SAH
				double cost = 2*T_aabb + (leftAreas[cutIndex]*axisSortedLeft[axis].size() + rightAreas[cutIndex]*axisSortedRight[axis].size())*currentBoxSurfaceInverse * T_tri;

				if ((cost < bestCost || bestCost == -1) && cost < leafCost)
				{
					bestCost = cost;
					bestIndex = cutIndex;
					bestAxis = axis;
				} else {
					int a = 10;
				}
			}
			delete [] leftAreas;
			delete [] rightAreas;
		}
		break;
		}
	}

	if (bestAxis == -1)
	{
		// make leaf
		MakeLeaf_internal(targets, index);
		return;
	}


	EntityList_t lefts = targets;
	std::sort(lefts.begin(), lefts.end(),
		[&bestAxis](Entity * const &a, Entity * const &b) -> bool {
		switch (bestAxis)
		{
		case 0: return a->boundingBox.position().x < b->boundingBox.position().x;
		case 1: return a->boundingBox.position().y < b->boundingBox.position().y;
		case 2: return a->boundingBox.position().z < b->boundingBox.position().z;
		}
		assert(false);
		return a->boundingBox.position().x < b->boundingBox.position().x;
	});

#ifdef USE_STXXL
	EntityList_t rights;
	for (EntityList_t::iterator k = lefts.begin()+bestIndex + 1; k != lefts.end(); ++k) rights.push_back(*k);
#else
	EntityList_t rights(lefts.begin() + bestIndex + 1, lefts.end());
#endif
	lefts.resize(bestIndex+1);

	BVH_structure *current = &m_root[index];
	current->axis = bestAxis;
	current->children[0] = m_root.size();
	current->children[1] = m_root.size()+1;
	m_root.push_back(BVH_structure()); m_root.push_back(BVH_structure());

	// constructs children
	Construct_internal(type, lefts, current->children[0]);
	Construct_internal(type, rights, current->children[1]);
}

void BVH::CollectBoundingBoxes(int depth, std::vector<BoundingBox> &result)
{
  result.clear();
  result.reserve(m_root.size());
  CollectBoundingBoxes_internal(0, depth, 0, result);
}

void BVH::CollectBoundingBoxes_internal(int currentDepth, int targetDepth, int index, std::vector<BoundingBox> &result)
{
  if (targetDepth < currentDepth) return;

  const BVH_structure *current = &m_root[index];

  if (targetDepth == currentDepth)
  {
    const BoundingBox box(current->box[0], current->box[1]);
    result.push_back(box);
    return;
  }

  if (current->children[0] == -1) return;
  CollectBoundingBoxes_internal(currentDepth+1, targetDepth, current->children[0], result);

  if (current->children[1] == -1) return;
  CollectBoundingBoxes_internal(currentDepth+1, targetDepth, current->children[1], result);
}

const BVH_structure *BVH::GetRootNode() const 
{
  return &m_root[0];
}

size_t BVH::GetBVHNodeCount() const 
{
  return m_root.size();
}

const BVH_structure *BVH::GetFirstChild(const BVH_structure *parent) const 
{
  assert (parent);
  if (IsLeaf(parent)) return NULL;
  assert (m_root.size() > parent->children[0]);
  return &m_root[parent->children[0]];
}

const BVH_structure *BVH::GetSecondChild(const BVH_structure *parent) const 
{
  assert (parent);
  if (IsLeaf(parent)) return NULL;
  assert (m_root.size() > parent->children[1]);
  return &m_root[parent->children[1]];
}

bool BVH::IsLeaf(const BVH_structure *node) const 
{
  assert (node);
  return node->children[0] == (-1);
}

}
