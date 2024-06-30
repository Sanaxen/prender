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
	// A Simple and Robust Mutation Strategy for the Metropolisを参照。
	// Kelemen style MLT用データ構造
	// Kelemen styleではパス生成に使う乱数の空間で変異させたりする。
	// その一つ一つのサンプルのデータ構造。

	// Metropolis Local Declarations
	struct PrimarySample {
		int modify_time;
		double value_;
		PrimarySample() {
			modify_time = 0;
			value_ = random()->next01();
		}
	};

	// LuxRenderから拝借してきた変異関数
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
	//論文「Simple and Robust Mutation Strategy for the Metropolis Light Transport Algorithm」と同じ
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



	// Kelemen MLTにおいて、パス生成に使う各種乱数はprimary spaceからもってくる。
	// PrimarySample()を通常のrand01()の代わりに使ってパス生成する。今回は普通のパストレースを使った。（双方向パストレ等も使える）
	// Metropolis法なので現在の状態から次の状態へと遷移をするがこの遷移の空間が従来のワールド空間ではなく
	// パスを生成するのに使われる乱数の空間になっている。
	// 乱数空間における変異（Mutate()）の結果を使って再び同じようにパストレでパスを生成するとそのパスは自然に変異後のパスになっている。
	struct KelemenMLT : public RandomMC{
	public:
		// 論文で使われているものとほぼ同じ
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

		// 通常の乱数のかわりにこれを使う
		// 論文にのっているコードとほぼ同じ	(PrimarySample)
		inline double NextSample() 
		{
			if (primary_samples.size() <= used_rand_coords) {
				primary_samples.resize(primary_samples.size() * 1.5); // 拡張する
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
			while (!primary_samples_stack.empty()) // スタック空にする
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

	// 上のパストレで生成したパスを保存しておく
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
		// 最初のパスを求める。輝度値に基づく重点サンプリングによって選んでいる。
		const double r = mlt.Next01() * sumI;	//輝度値サンプリング
		int selecetd_path = (int)(mlt.Next01()*seedPathNum);
		if (selecetd_path >= seedPathNum) selecetd_path = seedPathNum - 1;
		double accumulated_importance = 0.0;
		for (int i = 0; i < seedPathNum; i++)
		{
			accumulated_importance += luminance(seed_paths[i].F);
			if (r <= accumulated_importance)
			{
				//輝度値に基づく重点サンプリングで最初のパス(i)を求める
				selecetd_path = i;
				break;
			}
		}
		return selecetd_path;
	}
};

#endif