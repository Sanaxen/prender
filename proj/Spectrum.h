#ifndef _SPECTRUM_H__

#define _SPECTRUM_H__
//光のスペクトルに関するクラス
#include <algorithm>

#include "random.h"
#include "vector3d.h"
#include "material.h"
#include "constant.h"
#include "def.h"
#include <stdint.h>

/// スペクトルの範囲に注意する事!!

#define D_LINE_WAVELENGTH_SODIUM	0.5893//D-line wavelength of sodium[589.3nm]

#define WL_MAX	471

namespace prender 
{

//#ifdef SPECTRUM_USE
//typedef double Spectrum;
//#else
//typedef Color Spectrum;
//#endif

class spectrum
{

public:

	spectrum()
	{
	}

	//波長からXYZ値への変換
	double* Wavelength2xyzTable(void);

	// XYZ値をRGB値に変換する
	// CIE RGB
	inline Color xyz2rgb(const Color &xyz) const
	{
		return Color(dot(Color( 2.3655, -0.8971, -0.4683), xyz),
					 dot(Color(-0.5151,  1.4264,  0.0887), xyz),
					 dot(Color( 0.0052, -0.0144,  1.0089), xyz));
	}
	// RGB値をXYZ値に変換する
	// CIE RGB
	inline Color rgb2xyz(const Color &rgb) const
	{
		return Color(dot(Color( 0.4898, 0.3101, 0.2001), rgb),
					 dot(Color( 0.1769, 0.8124, 0.0107), rgb),
					 dot(Color( 0.0000, 0.0100, 0.9903), rgb));
	}

	inline double Wavelength(const Color &xyz) 
	{
		double *wavelength = Wavelength2xyzTable();//360nm - 830nm (5nm幅)
		double dmin = PS_INF;
		int index = -1;
		for ( int i = 0; i < WL_MAX; i++ )
		{
			Color xyz0(wavelength[3*i+1], wavelength[3*i+2], wavelength[3*i+3]);
			Color tmp = xyz - xyz0;
			double d = dot(tmp, tmp);
			if ( d < dmin )
			{
				dmin = d;
				index = i;
			}
		}
		double s = (360.0 + index*1.0) / 1000.0;;

		return s;
	}

	inline Color getColor(int wavelength_index) 
	{
		double *wavelength = Wavelength2xyzTable();//360nm - 830nm (1nm幅)

		return xyz2rgb(Color(
			wavelength[wavelength_index * 3 + 1], 
			wavelength[wavelength_index * 3 + 2],
			wavelength[wavelength_index * 3 + 3]));
	}

	inline Color getColor(double wavelength) 
	{
		double *wavelength_tbl = Wavelength2xyzTable();//360nm - 830nm (1nm幅)

		int wavelength_index = 470*(int)(std::min((wavelength - 0.36)/(0.83-0.360)+0.5,1.0));

		return xyz2rgb(Color(
			wavelength_tbl[wavelength_index * 3 + 1], 
			wavelength_tbl[wavelength_index * 3 + 2],
			wavelength_tbl[wavelength_index * 3 + 3]));
	}
};

#ifdef SPECTRUM_USE
Color Spectrum2RGB(const double lambda);
#endif

class spectrum_df
{
public:
	// 下で各波長についてパストレーシングするわけだが、波長の選択を完全にランダムにするのは効率が悪い。
	// 波長ごとに最終的な画像への寄与度は異なるからである。そこで、各波長の寄与度に基づいた重点サンプリングをする。
	// 具体的には各波長について、XYZ応答のうちY、すなわち輝度分に比例する確率密度関数を作り、それに基づいてサンプリングするようにする。
	// pdfはその密度関数でcdfはそこからサンプリングするための累積分布関数
	double cdf[WL_MAX];
	double pdf[WL_MAX];
	double luminance_table[WL_MAX];
	spectrum sp;

	spectrum_df()
	{

		double total = 0.0;
		const double* wavelength2xyz_table = sp.Wavelength2xyzTable();
		for (int i = 0; i < WL_MAX; i ++) 
		{
#if 0
			//XYZ系なら輝度はYの値
			luminance_table[i] = wavelength2xyz_table[i * 4 + 2];
#else
			//面倒だけど波長をRGBに直して「輝度」を算出してみる
#ifdef SPECTRUM_USE
			luminance_table[i] = luminance(Spectrum2RGB((i + 360.0) / 1000.0));
#else
			luminance_table[i] = 1.0;
#endif
#endif
			//fprintf(stderr, "%.16f\n", luminance_table[i]);
		}
		
		for (int i = 0; i < WL_MAX; i ++) {
			total += luminance_table[i];
			cdf[i] = total;
			pdf[i] = luminance_table[i];
		}
		for (int i = 0; i < WL_MAX; i ++) {
			cdf[i] /= total;
			pdf[i] /= total;
		}
	}
	inline double luminance(const Color &color) {
		return dot(Vector3d(0.2126, 0.7152, 0.0722), color);
	}


	double wavelength_sampling(Random& rnd, int& wavelength_index, double& pdf_div)
	{
#if 0
		// 波長を重点サンプリング
		double* p = std::lower_bound(cdf, cdf + WL_MAX, rnd.next01());
#ifdef _WIN64
		void* difp = (void*)(p - cdf);
		uint64_t x = reinterpret_cast<uintptr_t>(difp);
		//uint64_t x = reinterpret_cast<uintptr_t>(p) - reinterpret_cast<uintptr_t>(cdf);
		wavelength_index = x;
#else
		wavelength_index = p - cdf;
#endif
#else
		const double a = rnd.next01();
		for ( int i = 0; i < WL_MAX; i++ )
		{
			if ( cdf[i] >= a )
			{
				wavelength_index = i;
				break;
			}
		}
#endif
		if (wavelength_index >= WL_MAX) wavelength_index = WL_MAX - 1;
		const double div = pdf[wavelength_index];

		const double* wavelength2xyz_table = sp.Wavelength2xyzTable();
		//const double wavelength = wavelength2xyz_table[4*wavelength_index]/1000.0;
		double wavelength = (wavelength_index + 360.0)/1000.0;

		pdf_div = div;
		return wavelength;
	}
};

inline Color xyz2rgb(const Color &xyz)
{
#if 0
	Color rgb;
	//Reference White = D65
	// XYZ -> sRGB (D65)
	rgb.x = (xyz.x * ( 3.2404542) + xyz.y * (-1.5371385) + xyz.z * (-0.4985314));
	rgb.y = (xyz.x * (-0.9692660) + xyz.y * ( 1.8760108) + xyz.z * ( 0.0415560));
	rgb.z = (xyz.x * ( 0.0556434) + xyz.y * (-0.2040259) + xyz.z * ( 1.0572252));
	return rgb;
#else
	return Color(dot(Color(2.3655, -0.8971, -0.4683), xyz),
		dot(Color(-0.5151, 1.4264, 0.0887), xyz),
		dot(Color(0.0052, -0.0144, 1.0089), xyz));
#endif
}


//T temperature (Kelvin)
//lambda	wavelength (meter)
inline double BlackBody(  const double &T, const double &lambda) 
{
  static const double h = 6.62606896e-34;   // Plank constant
  static const double c = 2.99792458e+8;    // Speed of light
  static const double k = 1.38064880e-23;   // Boltzmann constant
  const double pl   = h * c / ( 4.97 * k * T );
  const double arg1 = 2 * PS_PI * h * c * c;
  const double arg2 = ( h * c ) / k;
  double value = (arg1 * pow(lambda, -5.0)) / (exp(arg2 / (lambda * T)) - 1.0);
  double peak  = (arg1 * pow(pl    , -5.0)) / (exp(arg2 / (pl     * T)) - 1.0);
  return value / peak;
}

inline double normalizedPlanck(double T, double lambda){ // K, nm
	static const double h = 6.62606896e-34;   // Plank constant
	static const double c = 2.99792458e+8;    // Speed of light
	static const double k = 1.38064880e-23;   // Boltzmann constant
	lambda = lambda / 1e9; // nm to m.
	double peakLambda = h * c / (4.97 * k * T);
	return BlackBody(T, lambda) / BlackBody(T, peakLambda);
}

#ifdef SPECTRUM_USE
double RGB2Spectrum(const Color rgb, const double lambda);
Color Spectrum2RGB(const double lambda);
Color Spectrum2XYZ(const double lambda);
double Color_lambda( int t, const double lambda);
double Color_lambda2( int t, const double lambda);

void TestSpectrum2XYZ();

void TestRedSpectrum();
void TestGreenSpectrum();
void TestBlueSpectrum();
void TestYellowSpectrum();
void TestMagentaSpectrum();
void TestCyanSpectrum();
void TestWhiteSpectrum();

void Test2RedSpectrum();
void Test2GreenSpectrum();
void Test2BlueSpectrum();
void Test2YellowSpectrum();
void Test2MagentaSpectrum();
void Test2CyanSpectrum();
void Test2WhiteSpectrum();

void TestSpectrum();
#endif

};
#endif