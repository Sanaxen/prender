#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "bitmap.h"
#include "Vector3d.h"
#include <memory.h>
#include "hdr.h"
#include "color.h"

#include <vector>
#include <string>

#include "ParticipatingMediaPrm.h"

namespace prender {


//材質
enum ReflectionType {
	REFLECTION_TYPE_DIFFUSE,		// 完全拡散面。いわゆるLambertian面。
	REFLECTION_TYPE_SPECULAR,		// 理想的な鏡面。
	REFLECTION_TYPE_REFRACTION,		// 理想的なガラス的物質。
	REFLECTION_TYPE_WARD_BRDF,		// Wardのグロシー
	REFLECTION_TYPE_PHONG_BRDF,		// Phong
	REFLECTION_TYPE_SSS_REFRACTION,				// 表面下散乱
	REFLECTION_TYPE_SSS_DIFFUSE,
	REFLECTION_TYPE_SSS_WARD_BRDF,
	REFLECTION_TYPE_REFRACTION_FRESNEL		// フレネル反射。
};

#define REFRACTIVE_INDEX_OF_GLASS	1.5			// ガラスの屈折率
#define REFRACTIVE_INDEX_OF_AIR		1.000292	// 空気の屈折率

enum RefractiveType
{
	REFRACTIVE_FORMULA0,	//Refractive index
	REFRACTIVE_FORMULA1,	//sqrt( a[0] - a[1]*pow(x,2) + a[2]*pow(x,-2) + a[3]*pow(x,-4) + a[4]*pow(x,-6) + a[5]*pow(x,-8) )
	REFRACTIVE_FORMULA2,	//sqrt( 1 + a[0]*pow(x,2)/(pow(x,2)-a[1]) + a[2]*pow(x,2)/(pow(x,2)-a[3]) + a[4]*pow(x,2)/(pow(x,2)-a[4]) )
	REFRACTIVE_FORMULA99,	//1.050000 +  9.5 / (wavelength - 156.0);
};

//Ward BRDFパラメータ
struct WardBRDFParameter
{
	double alp_x;
	double alp_y;
};

struct PhongBRDFParameter
{
	double specular_exponent;
};

class Material
{
public:
	std::string name;
	int id;

	bool background;

	inline Material()
	{
		id = -1;
		//name = "";	//初期化されているはずなので無駄
		background = false;
		formulaType = REFRACTIVE_FORMULA0;
		emission = Color();									//発光色
		color    = Color();									//反射率
		specular = Color();

		reflection_type = REFLECTION_TYPE_DIFFUSE;	//材質（表面における反射の種類）
		refractive_index = REFRACTIVE_INDEX_OF_GLASS;	
		//memset(Dispersion_formula1, '\0', sizeof(double)*6);
		Dispersion_formula1[0] = 0.0;
		Dispersion_formula1[1] = 0.0;
		Dispersion_formula1[2] = 0.0;
		Dispersion_formula1[3] = 0.0;
		Dispersion_formula1[4] = 0.0;
		Dispersion_formula1[5] = 0.0;
		r_refractive_index = 1.0;


		repeat = 1;
		texture = 0;
		bump_texture = 0;
		IBL_texture = 0;
		IBL_texture_HDR = 0;
		ibl_texture_coef = 1.0;
		hemisphere_map = 0;
		angular_map = 0;
		panoramic_map = 0;
		disk_map = 0;
		equirectangular_map = 0;

		roughness = 0.0;

	}
	inline ~Material()
	{
		//テクスチャはこのクラスが解放される解きに解放してはいけない(Entityが所持している)
		//if ( texture ) delete texture;
		//texture = 0;
	}

	//媒質のパラメータ
	ParticipatingMedia participatingMediaParam;

	Color emission;					//発光色
	Color color;					//反射率
	Color specular;
	double roughness;

	ReflectionType reflection_type;	//材質（表面における反射の種類）

	union 
	{
		WardBRDFParameter ward_brdf;	//Ward BRDFパラメータ
		PhongBRDFParameter phong_brdf;	//Phong
	} brdfParameter;

	int formulaType;				//屈折率算出方法
	double refractive_index;		//屈折率
	double Dispersion_formula1[6];	//波長に対する屈折率

	double r_refractive_index;		//相対屈折率

	inline double refractive(const double wavelength) const
	{
		//http://refractiveindex.info/legacy/
		const double* a = Dispersion_formula1;
		switch ( formulaType )
		{
		case REFRACTIVE_FORMULA0:
			return refractive_index;
			break;
		case REFRACTIVE_FORMULA1:
			return sqrt( a[0] + a[1]*pow(wavelength,2.0) + a[2]*pow(wavelength,-2.0) + a[3]*pow(wavelength,-4.0) + a[4]*pow(wavelength,-6.0) + a[5]*pow(wavelength,-8.0) );
			break;
		case REFRACTIVE_FORMULA2:
			return sqrt( 1.0 + a[0]*pow(wavelength,2.0)/(pow(wavelength,2.0)-a[1]) + a[2]*pow(wavelength,2.0)/(pow(wavelength,2.0)-a[3]) + a[4]*pow(wavelength,2.0)/(pow(wavelength,2.0)-a[5]) );
			break;
		case REFRACTIVE_FORMULA99:
			return 1.400 + 50.0 / (wavelength*1000.0 - 230.0);
			break;
		default:
			return refractive_index;
			break;
		}
		return refractive_index;
	}

	//正しくないrefractive(const double wavelength)をコピーした関数
	inline double r_refractive(const double wavelength) const
	{
		//http://refractiveindex.info/legacy/
		const double* a = Dispersion_formula1;
		switch ( formulaType )
		{
		case REFRACTIVE_FORMULA0:
			return r_refractive_index;
			break;
		case REFRACTIVE_FORMULA1:
			return sqrt( a[0] + a[1]*pow(wavelength,2.0) + a[2]*pow(wavelength,-2.0) + a[3]*pow(wavelength,-4.0) + a[4]*pow(wavelength,-6.0) + a[5]*pow(wavelength,-8.0) );
			break;
		case REFRACTIVE_FORMULA2:
			return sqrt( 1.0 + a[0]*pow(wavelength,2.0)/(pow(wavelength,2.0)-a[1]) + a[2]*pow(wavelength,2.0)/(pow(wavelength,2.0)-a[3]) + a[4]*pow(wavelength,2.0)/(pow(wavelength,2.0)-a[5]) );
			break;
		case REFRACTIVE_FORMULA99:
			return 1.400 + 50.0 / (wavelength*1000.0 - 230.0);
			break;
		default:
			return r_refractive_index;
			break;
		}
		return r_refractive_index;
	}

	BitMap* texture;			//テクスチャ
	BitMap* bump_texture;		//Bumpテクスチャ
	int repeat;

	// IBL
	BitMap* IBL_texture;
	HDRImage* IBL_texture_HDR;
	double ibl_texture_coef;

	char hemisphere_map;				//テクスチャ（IBL)のマッピング方法
	char angular_map;
	char panoramic_map;
	char disk_map;
	char equirectangular_map;

	//物体をライトとする条件判定
	inline bool isLight() const
	{
		if ( reflection_type == REFLECTION_TYPE_DIFFUSE && color.length() < 0.001 && emission.length() > 0.001 )
		{
			return true;
		}

#if IBL_LIGHT
		//IBLってライト扱いで良いのか?
		if (IBL())
		{
			return true;
		}
#endif
		return false;
	}

	///////////////// IBL関連 //////////////////////////////
	inline bool IBL() const
	{
		return ( IBL_texture || IBL_texture_HDR);
	}
	inline int IBL_W() const
	{
		if ( IBL_texture ) return IBL_texture->W();
		if ( IBL_texture_HDR ) return IBL_texture_HDR->width();
		return 0;
	}
	inline int IBL_H() const
	{
		if ( IBL_texture ) return IBL_texture->H();
		if ( IBL_texture_HDR ) return IBL_texture_HDR->height();
		return 0;
	}

	inline Color IBL_Color(int i, int j)
	{
		if ( IBL_texture )
		{
			Rgb rgb = IBL_texture->cell(i,j);
			return Color(rgb.r, rgb.g, rgb.b) / 255.0;
		}
		if ( IBL_texture_HDR ) return IBL_texture_HDR->image(j,i);
		return Color(0,0,0);
	}
};


class MaterialList
{
public:
	std::vector<Material> list;

	inline Material* getMaterial(const int id) 
	{
		return &(list[id]);
	}
	inline int FindMaterial(std::string& name) 
	{
		if ( name == "" )
		{
			return -1;
		}
		for ( int i = 0; i < list.size(); i++ )
		{
			if ( list[i].name == name )
			{
				return i;
			}
		}
		return -1;
	}

	inline Material* getMaterial(std::string& name) 
	{
		int id = FindMaterial( name );
		if ( id < 0 )
		{
			fprintf(stderr, "Material[%s] not found!!\n", name.c_str());
			return 0;
		}
		return &(list[id]);
	}


	inline int Add( Material mat )
	{
		if ( mat.id >= 0 )
		{
			return mat.id;
		}
		if ( mat.name != "" )
		{
			int id = -1;
			if ( (id = FindMaterial(mat.name)) == -1 )
			{
				mat.id = list.size();
				list.push_back(mat);
				return mat.id;
			}
			return id;
		}
		mat.id = list.size();
		list.push_back( mat );
		return mat.id;
	}

};

};

#endif
