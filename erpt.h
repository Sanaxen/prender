#ifndef _ERPT

#define _ERPT

#include "color.h"
#include "random.h"
#include "stack"

#include "def.h"

namespace prender {

	// Metropolis�@�Ȃ̂Ō��݂̏�Ԃ��玟�̏�ԂւƑJ�ڂ����邪���̑J�ڂ̋�Ԃ��]���̃��[���h��Ԃł͂Ȃ�
	// �p�X�𐶐�����̂Ɏg���闐���̋�ԂɂȂ��Ă���B
	// ������Ԃɂ�����ψفiMutate()�j�̌��ʂ��g���čĂѓ����悤�Ƀp�X�g���Ńp�X�𐶐�����Ƃ��̃p�X�͎��R�ɕψٌ�̃p�X�ɂȂ��Ă���B

#if 10
	struct PrimarySampleERPT {
		int  mutating;
		double value;
		PrimarySampleERPT()
		{
			mutating = 0;
			value = -1.0;
		}
	};
	struct ERPTSampler : public RandomMC{
		double MutateDistance;
		inline double mutate(double value) {
			value += MutateDistance * (2.0 * random()->next01() - 1.0);
			if (value > 1.0) value -= 1.0;
			if (value < 0.0) value += 1.0;
			return value;
		}


	public:
		int used_rand_coords;

		std::vector<PrimarySampleERPT> primary_samples;

		inline ERPTSampler()
		{
			used_rand_coords = 0;
			primary_samples.resize(32);
			for (int i = 0; i < 32; i++)
				primary_samples[i].value = random()->next01();
		}
		inline void reset() {
			used_rand_coords = 0;
			const int sz = primary_samples.size();
			for (int i = 0; i < sz; i++)
				primary_samples[i].mutating = 0;
		}

		// �ʏ�̗����̂����ɂ�����g��
		inline double NextSample()
		{
			if (primary_samples.size() <= used_rand_coords) {
				const int now_max = primary_samples.size();
				primary_samples.resize(primary_samples.size() * 1.5); // �g������
				for (int i = now_max; i < primary_samples.size(); i++)
					primary_samples[i].value = random()->next01();
			}
			used_rand_coords++;

			double value = primary_samples[used_rand_coords - 1].value;
			if (primary_samples[used_rand_coords - 1].mutating)
			{
				value = mutate(value);
			}
			return value;
		}

		inline void mutate()
		{
			const int sz = primary_samples.size();
			for (int i = 0; i < sz; i++)
				primary_samples[i].mutating = 1;
		}

		virtual double next01(){ return NextSample(); }
		virtual void seed(unsigned int initial_seed)
		{
			return;
		}
		virtual double Next01(){ return random()->next01(); }
	};
#else
	struct ERPTSampler :  public MLT{

		double MutateDistance;
		inline double mutate(double value) {
			value += MutateDistance * (2.0 * random()->next01() - 1.0);
			if (value > 1.0) value -= 1.0;
			if (value < 0.0) value += 1.0;
			return value;
		}

public:
	int used_rand_coords;

	std::vector<double> primary_samples;

	inline ERPTSampler() 
	{
		MutateDistance = 0.05;
		used_rand_coords = 0;
		primary_samples.resize(128);
		for (int i = 0; i < 128; i++)
			primary_samples[i] = random()->next01();
	}
	inline void reset() {
		used_rand_coords = 0;
	}

	// �ʏ�̗����̂����ɂ�����g��
	inline double NextSample() 
	{
		if (primary_samples.size() <= used_rand_coords) {
			const int now_max = primary_samples.size();
			primary_samples.resize(primary_samples.size() * 1.5); // �g������
			for (int i = now_max; i < primary_samples.size(); i++)
				primary_samples[i] = random()->next01();
		}
		used_rand_coords++;
		return primary_samples[used_rand_coords - 1];
	}

	inline void mutate()
	{
		const int sz = primary_samples.size();
		for (int i = 0; i < sz; i++)
			primary_samples[i] = Mutate(primary_samples[i]);
	}

	virtual double next01(){ return NextSample(); }
	virtual void seed(unsigned int initial_seed)
	{
		return;
	}
	virtual double Next01(){ return random()->next01(); }
};
#endif

};

#endif