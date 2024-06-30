#ifndef _homogeneous_h_
#define _homogeneous_h_

#include "entity.h"

namespace prender {

class Homogeneous
{
public:

	Entity* scattering_obj;			//散乱を考慮する要素
	int homogeneous_id;				//HomogeneousListのインデックス
};

class HomogeneousList
{
public:
	std::vector<Homogeneous> list;

	//要素が散乱を扱う要素か？
	inline bool isHomogeneous( const Entity* ent ) const
	{
		return ent->scattering;
	}
};

};

#endif