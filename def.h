#ifndef __DEF_H_
#define __DEF_H_

//#define SPECTRUM_USE	"FullSpectral_"	//Fullスペクトルレンダリング

#define BVH_USE		0
#define QBVH_USE	1

#define USE_SHADOW_RAY	1
#define SHADOW_RAY_SAMPLING	3

#define USE_PARTICIPATING_MEDIA	1		//USE_SHADOW_RAY と　SHADOW_RAY_SAMPLING　を有効にする！！

#define IBL_LIGHT	0					//IBLを光源（ライト）として扱う:1

#ifndef SPECTRUM_USE

#define RGB2Spectrum(x,y)	x
#define Spectrum2RGB(x)		1.0
#define RGB_MAX(color)	(std::max(color.x, std::max(color.y, color.z)))
#define MULTIPLY(x,y)	(multiply((x),(y)))
#define ZERO()			(Color(0.0,0.0,0.0))
#define ONE()			(Color(1.0,1.0,1.0))
#define DIFF(x,y)		(((x)-(y)).lenght())
#define SIZE(x)			((x).length())
#else
#define RGB_MAX(color)	(color)
#define MULTIPLY(x,y)	(x*y)
#define ZERO()			(0.0)
#define ONE()			(1.0)
#define DIFF(x,y)		(fabs((x) - (y)))
#define SIZE(x)			(fabs(x))
#endif

//四捨五入
inline double round_(double r) 
{
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}
inline double round_(float r) 
{
    return (r > 0.0f) ? floor(r + 0.5f) : ceil(r - 0.5f);
}

inline void Abort(char* msg, int line)
{
	if ( msg )
	{
		fprintf(stderr, "Abort:%s LINE:%d\n", msg, line);
		printf("Abort:%s LINE:%d\n", msg, line);
	}
	abort();
}

namespace prender 
{

#ifdef SPECTRUM_USE
typedef double Spectrum;
#else
typedef Color Spectrum;
#endif
};
#endif	// __DEF_H_
