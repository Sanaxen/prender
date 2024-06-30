#ifndef _MLT

#define _MLT

#include "color.h"
#include "random.h"
#include <stack>
#include <algorithm>
#include <functional>

#include "def.h"

namespace prender {

	extern std::vector<MTRandom*> PrimarySampleRondomVector;

	inline MTRandom* random()
	{
#ifdef _OPENMP
		return PrimarySampleRondomVector[omp_get_thread_num()];
#else
		return PrimarySampleRondomVector[0];
#endif

	}
	// A Simple and Robust Mutation Strategy for the Metropolis���Q�ƁB
	// Kelemen style MLT�p�f�[�^�\��
	// Kelemen style�ł̓p�X�����Ɏg�������̋�Ԃŕψق������肷��B
	// ���̈��̃T���v���̃f�[�^�\���B

	// Metropolis Local Declarations
	struct PrimarySample {
		int modify_time;
		double value_;
		PrimarySample() {
			modify_time = 0;
			value_ = random()->next01();
		}
	};

	// LuxRender����q�؂��Ă����ψي֐�
	inline double Mutate0(const double  x)
	{
		const double r = random()->next01();
		const double s1 = 1.0 / 512.0, s2 = 1.0 / 16.0;
		const double dx = s1 / (s1 / s2 + fabsf(2.0 * r - 1.0)) - s1 / (s1 / s2 + 1.0);
		if (r < 0.5) {
			const double x1 = x + dx;
			return (x1 < 1.0) ? x1 : x1 - 1.0;
		}
		else {
			const double x1 = x - dx;
			return (x1 < 0.0) ? x1 + 1.0 : x1;
		}
	}
	//�_���uSimple and Robust Mutation Strategy for the Metropolis Light Transport Algorithm�v�Ɠ���
	inline double Mutate1(double  x)
	{
		const float s1 = 1. / 1024, s2 = 1. / 64;
		const float dv = s2*expf(-logf(s2 / s1)*random()->next01());
		if (random()->next01() < 0.5) {
			x += dv; if (x > 1) x -= 1;
		}
		else {
			x -= dv; if (x < 0) x += 1;
		}
		return x;
	}

	class RandomMC : public Random
	{
	public:
		RandomMC(){}
		virtual ~RandomMC(){}
		virtual double NextSample(void) = 0;

	};



	// Kelemen MLT�ɂ����āA�p�X�����Ɏg���e�헐����primary space��������Ă���B
	// PrimarySample()��ʏ��rand01()�̑���Ɏg���ăp�X��������B����͕��ʂ̃p�X�g���[�X���g�����B�i�o�����p�X�g�������g����j
	// Metropolis�@�Ȃ̂Ō��݂̏�Ԃ��玟�̏�ԂւƑJ�ڂ����邪���̑J�ڂ̋�Ԃ��]���̃��[���h��Ԃł͂Ȃ�
	// �p�X�𐶐�����̂Ɏg���闐���̋�ԂɂȂ��Ă���B
	// ������Ԃɂ�����ψفiMutate()�j�̌��ʂ��g���čĂѓ����悤�Ƀp�X�g���Ńp�X�𐶐�����Ƃ��̃p�X�͎��R�ɕψٌ�̃p�X�ɂȂ��Ă���B
	struct KelemenMLT : public RandomMC{
	public:
		// �_���Ŏg���Ă�����̂Ƃقړ���
		int global_time;
		int large_step;
		int large_step_time;
		int used_rand_coords;

		std::vector<PrimarySample> primary_samples;
		std::stack<PrimarySample> primary_samples_stack;

		inline KelemenMLT()
		{
			global_time = large_step = large_step_time = used_rand_coords = 0;
			primary_samples.resize(128);
		}

		inline void InitUsedRandCoords() {
			used_rand_coords = 0;
		}

		// �ʏ�̗����̂����ɂ�����g��
		// �_���ɂ̂��Ă���R�[�h�Ƃقړ���	(PrimarySample)
		inline double NextSample() 
		{
			if (primary_samples.size() <= used_rand_coords) {
				primary_samples.resize(primary_samples.size() * 1.5); // �g������
			}

			PrimarySample& primary_sample = primary_samples[used_rand_coords];
			if (primary_sample.modify_time < global_time) {
				
				// large step
				if (large_step > 0) 
				{
					primary_samples_stack.push(primary_sample);
					primary_sample.modify_time = global_time;
					primary_sample.value_ = random()->next01();
				}
				else 
				{
					// small step
					if (primary_sample.modify_time < large_step_time) {
						primary_sample.modify_time = large_step_time;
						primary_sample.value_ = random()->next01();
					}

					// lazy evaluation of mutations
					while (primary_sample.modify_time < global_time - 1) {
						primary_sample.value_ = Mutate0(primary_sample.value_);
						primary_sample.modify_time++;
					}
					primary_samples_stack.push(primary_sample);
					primary_sample.value_ = Mutate0(primary_sample.value_);
					primary_sample.modify_time = global_time;
				}
			}

			used_rand_coords++;
			return  primary_samples[used_rand_coords - 1].value_;
		}

		inline void StackClear()
		{
			while (!primary_samples_stack.empty()) // �X�^�b�N��ɂ���
			  primary_samples_stack.pop();
		}
		inline void RestoreState()
		{
			int idx = used_rand_coords - 1;
			while (!primary_samples_stack.empty()) {
				primary_samples[idx --] = primary_samples_stack.top();
				primary_samples_stack.pop();
			}
		}

		virtual double next01(){ return NextSample(); }
		virtual void seed(unsigned int initial_seed)
		{
			return;
		}
		virtual double Next01(){ return random()->next01(); }
	};


	inline double luminance(const Color &color) {
		return dot(Vector3d(0.2126, 0.7152, 0.0722), color);
		//return color.Max();
	}

	// ��̃p�X�g���Ő��������p�X��ۑ����Ă���
	struct PathSample {
		int x, y;
		double weight;
		bool direct_hit;
		bool cancel;

		Color F;
		inline PathSample(const int x_ = 0, const int y_ = 0, const Color &F_ = ZERO(), const double weight_ = 0.0, const bool direct_hit_ = false) :
			x(x_), y(y_), F(F_), weight(weight_), direct_hit(direct_hit_)
		{
			if (F.x == 0.0 && F.y == 0.0 && F.z == 0.0) cancel = true;
			else cancel = false;
		}
	};
	inline bool operator<(const PathSample& left, const PathSample& right)
	{
		return luminance(left.F) < luminance(right.F) ;
	}

	inline bool operator>(const PathSample& left, const PathSample& right)
	{
		return luminance(left.F) > luminance(right.F) ;
	}

	inline int SPPMLTSamplingPath(KelemenMLT& mlt, std::vector<PathSample>& seed_paths, const double sumI, const int seedPathNum)
	{
		// �ŏ��̃p�X�����߂�B�P�x�l�Ɋ�Â��d�_�T���v�����O�ɂ���đI��ł���B
		const double r = mlt.Next01() * sumI;	//�P�x�l�T���v�����O
		int selecetd_path = (int)(mlt.Next01()*seedPathNum);
		if (selecetd_path >= seedPathNum) selecetd_path = seedPathNum - 1;
		double accumulated_importance = 0.0;
		for (int i = 0; i < seedPathNum; i++)
		{
			accumulated_importance += luminance(seed_paths[i].F);
			if (r <= accumulated_importance)
			{
				//�P�x�l�Ɋ�Â��d�_�T���v�����O�ōŏ��̃p�X(i)�����߂�
				selecetd_path = i;
				break;
			}
		}
		return selecetd_path;
	}
};

#endif