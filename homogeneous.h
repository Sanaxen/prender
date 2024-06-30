#ifndef _homogeneous_h_
#define _homogeneous_h_

#include "entity.h"

namespace prender {

class Homogeneous
{
public:

	Entity* scattering_obj;			//�U�����l������v�f
	int homogeneous_id;				//HomogeneousList�̃C���f�b�N�X
};

class HomogeneousList
{
public:
	std::vector<Homogeneous> list;

	//�v�f���U���������v�f���H
	inline bool isHomogeneous( const Entity* ent ) const
	{
		return ent->scattering;
	}
};

};

#endif