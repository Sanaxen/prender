#ifndef _SCENEOBJECT_H_

#define _SCENEOBJECT_H_

#include <vector>

#include "entity.h"

namespace prender {

class ScenObjct
{
public:

	ScenObjct(){}

	~ScenObjct()
	{
	}

	std::vector<Entity*> List;

	void Add(Entity* ent)
	{
		ent->CreateBoundingBox();
		ent->id = List.size();
		List.push_back(ent);
	}

};
}

#endif
