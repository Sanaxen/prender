#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <climits>
#include <time.h>
#include <vector>

//�X���b�h�Z�[�t���l������ꍇ
//#define MPTHREAD_CIRITICAL

extern std::vector<int> prime_numbers;

namespace prender {

class Random
{
public:
	Random(){}
	virtual ~Random(){}
	virtual double next01(void)=0;
	virtual double Next01(void)=0;
	virtual void seed(unsigned int s)=0;
};

// Xor-Shift�ɂ�闐���W�F�l���[�^
class XorShift : public Random
{
	unsigned int seed_[4];
public:
	inline unsigned int next(void) { 
		const unsigned int t = seed_[0] ^ (seed_[0] << 11);
		seed_[0] = seed_[1]; 
		seed_[1] = seed_[2];
		seed_[2] = seed_[3];
		return seed_[3] = (seed_[3] ^ (seed_[3] >> 19)) ^ (t ^ (t >> 8)); 
	}

	virtual double next01(void) {
		return (double)next() / UINT_MAX;
	}
	virtual double Next01(void) {
		return next01();
	}

	XorShift() {}

	XorShift(const unsigned int initial_seed) {
		seed(initial_seed);
	}
	virtual void seed(unsigned int initial_seed)
	{
		unsigned int s = initial_seed;
		for (int i = 1; i <= 4; i++){
			seed_[i - 1] = s = 1812433253U * (s ^ (s >> 30)) + i;
		}
	}
};

#if 0
#define USE_XorShift	1
typedef XorShift MTRandom;
#else
#define USE_XorShift	0

#define MT64

#ifndef MT64
#define MTRandom_N 624
#else
#define MTRandom_N 312
#endif

class MTRandom : public Random
{
	int N;
	unsigned long *mt; /* the array for the state Vector3dtor  */
	int mti; /* mti==N+1 means mt[N] is not initialized */

#ifndef MT64
	/* initializes mt[N] with a seed */
	void init_genrand(unsigned long s);

	/* initialize by an array with array-length */
	/* init_key is the array for initializing keys */
	/* key_length is its length */
	/* slight change for C++, 2004/2/26 */
	void init_by_array(unsigned long init_key[], int key_length);

	/* generates a random number on [0,0xffffffff]-interval */
	unsigned long genrand_int32(void);

	/* generates a random number on [0,0x7fffffff]-interval */
	long genrand_int31(void);

	/* generates a random number on [0,1]-real-interval */
	double genrand_real1(void);

	/* generates a random number on [0,1)-real-interval */
	double genrand_real2(void);

	/* generates a random number on (0,1)-real-interval */
	double genrand_real3(void);

	/* generates a random number on [0,1) with 53-bit resolution*/
	double genrand_res53(void);
	/* These real versions are due to Isaku Wada, 2002/01/09 added */
#else
	void init_genrand64(unsigned long long seed);

	/* initialize by an array with array-length */
	/* init_key is the array for initializing keys */
	/* key_length is its length */
	void init_by_array64(unsigned long long init_key[], unsigned long long key_length);

	/* generates a random number on [0, 2^64-1]-interval */
	unsigned long long genrand64_int64(void);

	/* generates a random number on [0, 2^63-1]-interval */
	long long genrand64_int63(void);

	/* generates a random number on [0,1]-real-interval */
	double genrand64_real1(void);

	/* generates a random number on [0,1)-real-interval */
	double genrand64_real2(void);

	/* generates a random number on (0,1)-real-interval */
	double genrand64_real3(void);
#endif

public:
	MTRandom()
	{
		N = MTRandom_N;
		mt = new unsigned long[N+3];
		mt[N] = '#';
		mt[N+1] = '#';
		mt[N+2] = '#';
		mti=N+1; /* mti==N+1 means mt[N] is not initialized */
#ifndef MT64
		init_genrand((unsigned)time(NULL));
#else
		init_genrand64((unsigned)time(NULL));
#endif
	}
	MTRandom(const unsigned int initial_seed)
	{
		N = MTRandom_N;
		mt = new unsigned long[N+3];
		mt[N] = '#';
		mt[N + 1] = '#';
		mt[N + 2] = '#';
		seed(initial_seed);
	}

	virtual void seed(unsigned int initial_seed)
	{
#ifndef MT64
		init_genrand((unsigned)initial_seed);
#else
		init_genrand64((unsigned)initial_seed);
#endif
	}

	virtual double next01(void)
	{
#ifndef MT64
		return genrand_real1();
#else
		return genrand64_real1();
#endif
	}
	virtual double Next01(void) {
		return next01();
	}

	~MTRandom()
	{
		if (mt[N] != '#' || mt[N + 1] != '#' || mt[N + 2] != '#')
		{
			fprintf(stderr, "fata error MEMORY BOUNDARY OVER!!.");
			printf( "fata error MEMORY BOUNDARY OVER!!.");
		}
		delete [] mt;
	}
};

#endif
//typedef XorShift Random;
//typedef MTRandom RandomSTD;



void init_prime_numbers();
class QMC : public Random{
private:
	MTRandom rnd;
public:
	int j; // Halton���j�ڂ̒l�𓾂�
	int ith; // i�Ԗڂ̃T���v��
	QMC(const int ith_) : ith(ith_){
		j = 0;
	}
	QMC()
	{
		j = 0;
	}

	virtual void seed(unsigned int initial_seed)
	{
		ith = (int)initial_seed;
	}

	// Halton��𓾂�
	// Radical inverse function 
	//   ��_b(i) = 0.d1d2d3...dk...
	// i�Ԗڂ�n����Halton�� 
	//   x_i = (��_2(i), ��_3(i),..., ��_pj(i),..., ��_pn(i)), pn = n�Ԗڂ̑f��
	//

	// i�Ԗڂ�n����Halton���j�ڂ̗v�f�𓾂�
	// �܂��_pj(i)���v�Z����
	virtual  double next01()
	{
		// �ȉ��A��_b(i)���v�Z����
		int base = 0;
#ifdef MPTHREAD_CIRITICAL
#pragma omp critical
#endif
		{
			if (j < prime_numbers.size()) base = prime_numbers[j];
			j++;
		}
		
		// �f���e�[�u���𒴂��鎟���̃T���v���𓾂邱�Ƃ͂ł��Ȃ��I
		// ���V�A�����[���b�g�����܂�ɂ����܂��������ꍇ�ȂǁA�����ɓ��邱�Ƃ����蓾��i���ɒ�m�������j
		// ��ɃX�^�b�N�I�[�o�[�t���[���邩������Ȃ�
		if (base == 0)
		{
			return rnd.next01();
		}

		double value = 0.0;
		double inv_base = 1.0 / (double)base;
		double factor = inv_base;
		int i = ith;

		while (i > 0) {
			// �K����d_k�����ւ���	
			//const int d_k = reverse(i % base, base);
			//value += d_k * factor;
			value += (double)(i % base) * factor;
			i /= base;
			factor *= inv_base;
		}
		return value;
	}
	virtual double Next01(void) {
		return next01();
	}

	inline int reverse(const int i, const int p) {
		if (i == 0)
			return i;
		else
			return p - i;
	}
};

//typedef QMC Random;
//typedef MTRandom Random0;


};

#endif
