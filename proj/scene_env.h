#ifndef _SCENE_ENV_H_

#define _SCENE_ENV_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "def.h"
#include "Vector3d.h"
#include "matrix.h"
#include "sceneObject.h"

#include <omp.h>
#include <iostream>
#include "random.h"
#include "light.h"
#include "ParticipatingMediaPrm.h"

#include "KerrBlackHole.h"
#include "WormHole.h"

#include "Expr.h"
#define atof(s)	(expr_eval(s))
#define atoi(s)	((int)expr_eval(s))

inline void sscanf1f(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 1; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(buf, val);
}
inline void sscanf2f(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 2; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(buf, val);
}

inline void sscanf3f(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 3; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(buf, val);
}
inline void sscanf4f(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 4; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(buf, val);
}
inline void sscanf1fs(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	val.push_back(va_arg(args, double*));
	va_end(args);
	expr_eval(strchr(buf, ' ')+1, val);
}
inline void sscanf3fs(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 3; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(strchr(buf, ' ')+1, val);
}
inline void sscanf4fs(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 4; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(strchr(buf, ' ')+1, val);
}

inline void sscanf5fs(const char* buf, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

	std::vector<double*> val;
	for (int i = 0; i < 5; i++)
	{
		val.push_back(va_arg(args, double*));
	}
	va_end(args);
	expr_eval(strchr(buf, ' ') + 1, val);
}


#pragma warning(disable:4101)


namespace prender {

void BlackHoleTextureInit(char* drive, char* dir, KerrBlackHole* kb);
void WormHoleTextureInit(char* drive, char* dir, WormHole* kb);

inline char* getLine(char* buf, int size, FILE* fp)
{
	char* p = fgets(buf, size, fp);
	if ( p && p[0] =='#' )
	{
		return fgets(buf, size, fp);
	}
	return p;
}

// シーンの背景
extern BitMap* backgroundTexture[2];
extern SphericalMapping* backgroundMap[2];

class SceneEnv
{
	char drive[_MAX_DRIVE];	// ドライブ名
    char dir[_MAX_DIR];		// ディレクトリ名
    char fname[_MAX_FNAME];	// ファイル名
    char ext[_MAX_EXT];		// 拡張子

public:
	int threads;
	char output[512];

	int image_width;
	int image_height;
	int imageDumpTime;

	double world_screen_width;
	double world_screen_height;
	double world_screen_dist;
	Vector3d camera_position;
	Vector3d camera_dir;
	Vector3d target_position;	//見る位置
	Vector3d camera_up;
	double camera_angle;		//水平方向の視角。 0°～180°の範囲で指定する。
	Matrix4x4	camera_matrix;

	double gamma_offset = 0;
	int samples;
	int supersamples;
	int DepthLimit;
	int Depth;		 // ロシアンルーレットで打ち切らない最大深度

	Color BackgroundColor;
	Color Environment_Light;

	double lensRadius;
	double focalDistance;

	double sensor_response;	//センサ応答

	LightList light_list;				//光源リスト
	MaterialList material_list;			//材質リスト

	int nextEventEstimation;
	int participatingMedia;

	//関与媒質のパラメータ
	ParticipatingMedia participatingMediaParam;

	//メトロポリス光輸送法
	int metropolisTransport;
	double mutation;		//1ピクセル当たり の変異回数
	double pre_sample;		//サンプル数=全ピクセル数*pre_sample
	double mlt_r1;				//摂動パラメータ（MLT の元論文参照）
	double mlt_r2;				//摂動パラメータ（MLT の元論文参照）

	//Energy Redistribution Path Tracing
	int energyRedistributionPathTracing;

	MTRandom rand_;

	double timeLimit;

	//IBL
	int ibl;
	Entity* ibl_ent;

	float luminance_cutoff = 0.025;
	// tent filter
	int useTentFilter;

	//ブラックホール
	KerrBlackHole* blackHole;
	Vector3d blackHolePos;
	bool blackHolePos_set;
	bool blackHole_exist;
	int blackhole_disk;
	double blackhole_dist;
	double geodesics_max_length;
	double accretion_disk;
	double angular_momentum;
	int initial_condition;
	std::vector<std::string> accretion_disk_texture;
	std::string background_texture[2];
	double background_texture_map_coef[2][4];
	double background_texture_coef[2];

	//ワームホール
	WormHole* wormhole_;
	Vector3d wormholePos1;
	Vector3d wormholePos2;
	double wormhole_Rho;
	double wormhole_a;
	double wormhole_W;
	double wormhole_dist;
	double wormhole_sign;
	bool wormHole_exist;
	bool wormHolePos_set;

	ScenObjct EntList;
	void Default()
	{
		initial_condition = 0;
		angular_momentum = 0.0;
		accretion_disk = 0.0;
		geodesics_max_length = 10000.0;
		blackhole_dist = 10000.0;
		blackHolePos_set = false;
		blackHole_exist = false;
		blackhole_disk = false;

		background_texture_map_coef[0][0] = 1.0;
		background_texture_map_coef[0][1] = 0.0;
		background_texture_map_coef[0][2] = 1.0;
		background_texture_map_coef[0][3] = 0.0;

		background_texture_map_coef[1][0] = 1.0;
		background_texture_map_coef[1][1] = 0.0;
		background_texture_map_coef[1][2] = 1.0;
		background_texture_map_coef[1][3] = 0.0;

		background_texture_coef[0] = 1.0;
		background_texture_coef[1] = 1.0;

		wormhole_dist = 10000.0;
		wormHole_exist = false;
		wormHolePos_set = false;
		wormhole_Rho = 1.3;
		wormhole_a = 0.5;
		wormhole_W = 0.05;
		wormhole_sign = 0.0;

		blackHole = 0;
		wormhole_ = 0;

		imageDumpTime = 60;
		rand_.seed(time(NULL));

		sensor_response = 1.0;
		energyRedistributionPathTracing = 0;

		metropolisTransport = 0;
		mutation = 2;
		pre_sample = 0.4;
		mlt_r1 = 0.0;
		mlt_r2 = 0.0;

		useTentFilter = 0;

		ibl = -1;
		ibl_ent = 0;

		timeLimit = -1.0;
		threads = std::min(1, omp_get_max_threads()-1);
		if (threads < 1) threads = 1;
		lensRadius = -1.0;
		focalDistance = 0.0;

		BackgroundColor = Color(0.0, 0.0, 0.0);
		image_width = 320;
		image_height = 240;
		world_screen_width = 30.0;
		world_screen_height = 30.0;
		world_screen_dist = 40.0;
		
		samples = 8;
		supersamples = 2;
		DepthLimit = 64;
		Depth = 8;

		camera_position = Vector3d(50.0, 52.0, 220.0);
		camera_dir      = normalize(Vector3d(0.0, -0.04, -1.0));
		camera_up       = Vector3d(0.0, 1.0, 0.0);
		target_position = Vector3d(PS_INF,PS_INF,PS_INF);
		camera_angle = PS_INF;
		camera_matrix.LoadIdentity();

		nextEventEstimation = 0;
		participatingMedia  = 0;
		Environment_Light = Color(0.0,0.0,0.0);
	}

	SceneEnv()
	{
		Default();
	};

	void Load( char* filename)
	{
		Default();

		std::vector<Matrix4x4> camera_matlist;
		FILE* fp = fopen(filename, "r");

		if ( fp )
		{
			_splitpath( filename, drive, dir, fname, ext );

			std::cout << "Drive=" << drive << std::endl;
			std::cout << "Dir  =" << dir   << std::endl;
			std::cout << "Fname=" << fname << std::endl;
			std::cout << "Ext  =" << ext   << std::endl;

			char* use_brssdf = getenv("USE_BSSRDF");

			strcpy(output, "output.ppm");
			char buf[1024];
			while( getLine(buf, 1024, fp) != NULL )
			{
				//gamma_offset
				if (strcmp(buf, "gamma_offset\n") == 0)
				{
					getLine(buf, 1024, fp);
					gamma_offset = atof(buf);
					continue;
				}
				//LUMINANCE_CUTOFF
				if (strcmp(buf, "luminance_cutoff\n") == 0)
				{
					getLine(buf, 1024, fp);
					luminance_cutoff = atof(buf);
					continue;
				}
				if (strcmp(buf, "tentfilter\n") == 0)
				{
					getLine(buf, 1024, fp);
					useTentFilter = atoi(buf);
					continue;
				}
				if (strcmp(buf, "IMAGEDUMP\n") == 0)
				{
					getLine(buf, 1024, fp);
					imageDumpTime = atoi(buf);
					continue;
				}
				
				if (strncmp(buf, "VAR ", 4) == 0)
				{
					char* p = strchr(buf, '\n');
					if (p) *p = '\0';
					std::string name = buf + 4;
					getLine(buf, 1024, fp);
					double val = atof(buf);
					set_Variable(name, val);
					continue;
				}
				if (strcmp(buf, "TIMELIMIT\n") == 0)
				{
					getLine(buf, 1024, fp);
					timeLimit = atof(buf);
					continue;
				}
				if (strcmp(buf, "THREAD\n") == 0)
				{
					getLine(buf, 1024, fp);
					threads = atoi(buf);
					continue;
				}
				if ( strcmp(buf, "DEPTH_LIMIT\n") == 0 )
				{
					getLine(buf, 1024, fp);
					DepthLimit = atoi(buf);
					continue;
				}
				if ( strcmp(buf, "DEPTH\n") == 0 )
				{
					getLine(buf, 1024, fp);
					Depth = atoi(buf);
					continue;
				}
				if ( strcmp(buf, "OUTPUT\n") == 0 )
				{
					getLine(buf, 1024, fp);
					char* p = strchr(buf, '\n');
					if ( p ) *p = '\0';
					strcpy(output, buf);
					continue;
				}
				if ( strcmp(buf, "IMAGE\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d %d", &image_width, &image_height);
					continue;
				}
				if ( strcmp(buf, "SCREEN\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &world_screen_width, &world_screen_height, &world_screen_dist);
					continue;
				}
				if ( strcmp(buf, "SAMPLING\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d", &samples);
					continue;
				}
				if ( strcmp(buf, "SUPERSAMPLING\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d", &supersamples);
					continue;
				}
				if ( strcmp(buf, "CAMERA_POS\n") == 0 )
				{
					double x, y, z;
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &x, &y, &z);
					camera_position = Vector3d(x, y, z);
					continue;
				}
				if ( strcmp(buf, "CAMERA_DIR\n") == 0 )
				{
					double x, y, z;
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &x, &y, &z);
					camera_dir = Vector3d(x, y, z);
					continue;
				}
				if ( strcmp(buf, "TARGET_POS\n") == 0 )
				{
					double x, y, z;
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &x, &y, &z);
					target_position = Vector3d(x, y, z);
					continue;
				}
				if ( strcmp(buf, "CAMERA_UPVEC\n") == 0 )
				{
					double x, y, z;
					getLine(buf, 1024, fp);
					sscanf(buf, "%lf %lf %lf", &x, &y, &z);
					camera_up = Vector3d(x, y, z);
					continue;
				}
				if ( strcmp(buf, "CAMERA_ANGLE\n") == 0 )
				{
					getLine(buf, 1024, fp);
					camera_angle = atof(buf);
					continue;
				}
				if ( strcmp(buf, "CAMERA_MATRIX\n") == 0 )
				{
					double x, y, z;
					do{
						getLine(buf, 1024, fp);
						if ( buf[0] == '\n'||buf[0] == '\0') break;

						if ( strncmp(buf, "rotation ", 9) == 0 )
						{
							sscanf3fs(buf, "rotation %lf %lf %lf", &x, &y, &z);
							Matrix4x4 mat;
							camera_matlist.push_back(mat.RotationX(x*PS_RAD));
							camera_matlist.push_back(mat.RotationY(y*PS_RAD));
							camera_matlist.push_back(mat.RotationZ(z*PS_RAD));
							continue;
						}
						if ( strncmp(buf, "translate ", 10) == 0 )
						{
							sscanf3fs(buf, "translate %lf %lf %lf", &x, &y, &z);
							Matrix4x4 mat;
							camera_matlist.push_back(mat.Translation(x, y, z));
							continue;
						}
					}while(1);
					continue;
				}
				if ( strcmp(buf, "LENS\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf2f(buf, "%lf %lf", &lensRadius, &focalDistance);
					continue;
				}

				if ( strcmp(buf, "ENVLIGHT\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &Environment_Light.x, &Environment_Light.y, &Environment_Light.z);
					continue;
				}

				if ( strcmp(buf, "SCATTERING\n") == 0 )
				{
					getLine(buf, 1024, fp);
					participatingMediaParam.scattering = atof(buf);
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &participatingMediaParam.scatteringColor.x, &participatingMediaParam.scatteringColor.y, &participatingMediaParam.scatteringColor.z);
					continue;
				}
				if ( strcmp(buf, "ABSORBING\n") == 0 )
				{
					getLine(buf, 1024, fp);
					participatingMediaParam.absorbing = atof(buf);
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &participatingMediaParam.absorbingColor.x, &participatingMediaParam.absorbingColor.y, &participatingMediaParam.absorbingColor.z);
					continue;
				}
				if ( strcmp(buf, "PAHASE\n") == 0 )
				{
					getLine(buf, 1024, fp);
					participatingMediaParam.phase_prm[0] = atof(buf);
					participatingMediaParam.phase_prm[1] = participatingMediaParam.phase_prm[0];
					participatingMediaParam.phase_prm[2] = participatingMediaParam.phase_prm[0];
					continue;
				}
				if (strncmp(buf, "PAHASE ", 7) == 0)
				{
					double x, y, z;
					sscanf3fs(buf, "PAHASE %lf %lf %lf", &x, &y, &z);
					participatingMediaParam.phase_prm[0] = x;
					participatingMediaParam.phase_prm[1] = y;
					participatingMediaParam.phase_prm[2] = z;
					continue;
				}

				if ( strcmp(buf, "BACKGROUND\n") == 0 )
				{
					double x, y, z;
					getLine(buf, 1024, fp);
					sscanf3f(buf, "%lf %lf %lf", &BackgroundColor.x, &BackgroundColor.y, &BackgroundColor.z);
					continue;
				}
				if ( strcmp(buf, "nextEventEstimation\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d", &nextEventEstimation);
					continue;
				}
				if ( strcmp(buf, "participatingMedia\n") == 0 )
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d", &participatingMedia);
					continue;
				}
				if (strcmp(buf, "metropolisTransport\n") == 0)
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d", &metropolisTransport);
					continue;
				}
				if (strcmp(buf, "mutation\n") == 0)
				{
					getLine(buf, 1024, fp);
					sscanf1f(buf, "%lf", &mutation);
					continue;
				}
				if (strcmp(buf, "pre_sample\n") == 0)
				{
					getLine(buf, 1024, fp);
					sscanf1f(buf, "%lf", &pre_sample);
					continue;
				}
				if (strcmp(buf, "mlt_r1r2\n") == 0)
				{
					getLine(buf, 1024, fp);
					sscanf2f(buf, "%lf %lf", &mlt_r1, &mlt_r2);
					continue;
				}

				if (strcmp(buf, "energyRedistributionPathTracing\n") == 0)
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%d", &energyRedistributionPathTracing);
					continue;
				}
				if (strcmp(buf, "sensor_response\n") == 0)
				{
					getLine(buf, 1024, fp);
					sscanf(buf, "%lf", &sensor_response);
					continue;
				}

				if (strncmp(buf, "background_texture ", 19) == 0)
				{
					char* p = strchr(buf, '\n');
					if (p) *p = '\0';
					background_texture[0] = buf + 19;
					printf("background_texture=[%s]\n", background_texture[0].c_str());
					continue;
				}
				if (strncmp(buf, "map_coef ", 9) == 0)
				{
					sscanf4f(buf + 9, "%lf %lf %lf %lf", &(background_texture_map_coef[0][0]), &(background_texture_map_coef[0][1]), &(background_texture_map_coef[0][2]), &(background_texture_map_coef[0][3]));
					continue;
				}
				if (strncmp(buf, "coef ", 5) == 0)
				{
					sscanf1f(buf + 5, "%lf", &(background_texture_coef[0]));
					continue;
				}



				if (strncmp(buf, "BLACKHOLE", 9) == 0 )
				{
					blackHole_exist = true;
					if (strncmp(buf, "BLACKHOLE ", 10) == 0)
					{
						sscanf3fs(buf, "BLACKHOLE %lf %lf %lf", &blackHolePos.x, &blackHolePos.y, &blackHolePos.z);
						blackHolePos_set = true;
					}
					//属性定義
					while (fgets(buf, 1024, fp) != NULL)
					{
						if (buf[0] == '\n')
						{
							break;
						}
						if (buf[0] == '#')
						{
							continue;
						}
						if ( strncmp(buf, "dist ", 5) == 0 )
						{
							sscanf1fs(buf, "dist %lf", &blackhole_dist);
						}							
						if (strncmp(buf, "geodesics_max_length ", 21) == 0)
						{
							sscanf1fs(buf, "geodesics_max_length %lf", &geodesics_max_length);
						}
						if (strncmp(buf, "accretion_disk ", 15) == 0)
						{
							sscanf1fs(buf, "accretion_disk %lf", &accretion_disk);
						}
						if (strncmp(buf, "angular_momentum ", 17) == 0)
						{
							sscanf1fs(buf, "angular_momentum %lf", &angular_momentum );
						}
						if (strncmp(buf, "initial_condition ", 18) == 0)
						{
							double dmy;
							sscanf1fs(buf, "initial_condition %lf", &dmy);
							initial_condition = (int)dmy;
						}
						if (strncmp(buf, "accretion_disk_texture ", 23) == 0)
						{
							char* p = strchr(buf, '\n');
							if (p) *p = '\0';
							accretion_disk_texture.push_back( buf + 23 );
							printf("accretion_disk_texture[%d]=[%s]\n", (int)accretion_disk_texture.size() - 1, accretion_disk_texture[accretion_disk_texture.size() - 1].c_str());
						}
						if (strncmp(buf, "background_texture ", 19) == 0)
						{
							char* p = strchr(buf, '\n');
							if (p) *p = '\0';
							background_texture[0] = buf + 19;
							printf("background_texture=[%s]\n", background_texture[0].c_str());
						}
						if (strncmp(buf, "map_coef ", 9) == 0)
						{
							sscanf4f(buf+9, "%lf %lf %lf %lf", &(background_texture_map_coef[0][0]), &(background_texture_map_coef[0][1]), &(background_texture_map_coef[0][2]), &(background_texture_map_coef[0][3]));
						}
						if (strncmp(buf, "coef ", 5) == 0)
						{
							sscanf1f(buf+5, "%lf", &(background_texture_coef[0]));
						}
					}
					continue;
				}

				if (strncmp(buf, "WORMHOLE", 8) == 0)
				{
					wormHole_exist = true;
					//属性定義
					while (fgets(buf, 1024, fp) != NULL)
					{
						if (buf[0] == '\n')
						{
							break;
						}
						if (buf[0] == '#')
						{
							continue;
						}
						if (strncmp(buf, "in ", 3) == 0)
						{
							sscanf3fs(buf, "in %lf %lf %lf", &wormholePos1.x, &wormholePos1.y, &wormholePos1.z);
							wormHolePos_set = true;
						}
						if (strncmp(buf, "out ", 4) == 0)
						{
							sscanf3fs(buf, "out %lf %lf %lf", &wormholePos2.x, &wormholePos2.y, &wormholePos2.z);
							wormHolePos_set = true;
						}
						if (strncmp(buf, "Rho ", 4) == 0)
						{
							sscanf1fs(buf, "Rho %lf %lf %lf", &wormhole_Rho);
						}
						if (strncmp(buf, "a ", 2) == 0)
						{
							sscanf1fs(buf, "a %lf %lf %lf", &wormhole_a);
						}
						if (strncmp(buf, "W ", 2) == 0)
						{
							sscanf1fs(buf, "W %lf %lf %lf", &wormhole_W);
						}
						if (strncmp(buf, "sign ", 5) == 0)
						{
							sscanf1fs(buf, "sign %lf", &wormhole_sign);
						}


						if (strncmp(buf, "dist ", 5) == 0)
						{
							sscanf1fs(buf, "dist %lf", &wormhole_dist);
						}
						if (strncmp(buf, "geodesics_max_length ", 21) == 0)
						{
							sscanf1fs(buf, "geodesics_max_length %lf", &geodesics_max_length);
						}
						if (strncmp(buf, "background_texture1 ", 20) == 0)
						{
							char* p = strchr(buf, '\n');
							if (p) *p = '\0';
							background_texture[0] = buf + 20;
							printf("background_texture1=[%s]\n", background_texture[0].c_str());
						}
						if (strncmp(buf, "background_texture2 ", 20) == 0)
						{
							char* p = strchr(buf, '\n');
							if (p) *p = '\0';
							background_texture[1] = buf + 20;
							printf("background_texture2=[%s]\n", background_texture[1].c_str());
						}
						if (strncmp(buf, "map_coef1 ", 10) == 0)
						{
							sscanf4f(buf+10, "%lf %lf %lf %lf", &(background_texture_map_coef[0][0]), &(background_texture_map_coef[0][1]), &(background_texture_map_coef[0][2]), &(background_texture_map_coef[0][3]));
						}
						if (strncmp(buf, "map_coef2 ", 10) == 0)
						{
							sscanf4f(buf+10, "%%lf %lf %lf %lf", &(background_texture_map_coef[1][0]), &(background_texture_map_coef[1][1]), &(background_texture_map_coef[1][2]), &(background_texture_map_coef[1][3]));
						}
						if (strncmp(buf, "coef1 ", 6) == 0)
						{
							sscanf1f(buf+6, "%lf", &(background_texture_coef[0]));
						}
						if (strncmp(buf, "coef2 ", 6) == 0)
						{
							sscanf1f(buf+6, "%lf", &(background_texture_coef[1]));
						}

					}
					continue;
				}

				//形状定義
				if ( strcmp(buf, "OBJECT\n") == 0 )
				{
					std::vector<Matrix4x4> matrixList;

					EntityType type;
					char file[512];			//OBJファイル（メッシュ）

					double x, y, z, r, rr;		//球
					int hemisphere = 0;		//半球フラグ
					Vector3d org, normal;	//平面
					Vector3d U, V;	//UV平面

					double dx, dy, dz;		//平行移動
					double sx, sy, sz;		//スケール
					Material material;		//材質
					double material_scale = 1.0;

					int smooth = 0;
					bool background = false;

					int normal_vector_inverse = 1;
					int back = 1;
					int material_id = -1;

					int parallel_light = 0;
					Vector3d light_dir;

					int spot_light = 0;
					double spot_angle[2] = { 0.0, 0.0 };
					double falloff = 0.0;
					double attenuation[3] = { 0.0, 0.0, 1.0 };

					int shadow = 0;
					int media_ignore_front_back = '\0';

					getLine(buf, 1024, fp);
					
					if ( strncmp(buf, "infinity_light ", 15) == 0 )
					{
						sscanf3fs(buf, "infinity_light %lf %lf %lf", &normal.x, &normal.y, &normal.z);
						type = ENTITY_TYPE_INF_LIGHT;
					}

					//無限平面
					if ( strncmp(buf, "plane ", 6) == 0 )
					{
						sscanf3fs(buf, "plane %lf %lf %lf", &org.x, &org.y, &org.z);
						getLine(buf, 1024, fp);
						sscanf3fs(buf, "normal %lf %lf %lf", &normal.x, &normal.y, &normal.z);
						type = ENTITY_TYPE_PLANE;
					}

					//有限平面
					if ( strncmp(buf, "uvplane ", 8) == 0 )
					{
						sscanf3fs(buf, "uvplane %lf %lf %lf", &org.x, &org.y, &org.z);
						getLine(buf, 1024, fp);
						if ( buf[0] == 'U' )
						{
							sscanf3fs(buf, "U %lf %lf %lf", &U.x, &U.y, &U.z);
						}else
						{
							sscanf3fs(buf, "V %lf %lf %lf", &V.x, &V.y, &V.z);
						}
						getLine(buf, 1024, fp);
						if ( buf[0] == 'U' )
						{
							sscanf3fs(buf, "U %lf %lf %lf", &U.x, &U.y, &U.z);
						}else
						{
							sscanf3fs(buf, "V %lf %lf %lf", &V.x, &V.y, &V.z);
						}
						type = ENTITY_TYPE_UVPLANE;
					}

					//球
					if ( strncmp(buf, "sphere ", 7) == 0 )
					{
						sscanf4fs(buf, "sphere %lf %lf %lf %lf", &x, &y, &z, &r);
						type = ENTITY_TYPE_SPHERE;
					}
					if ( strncmp(buf, "hemisphere ", 11) == 0 )
					{
						sscanf4fs(buf, "hemisphere %lf %lf %lf %lf", &x, &y, &z, &r);
						hemisphere = 1;
						type = ENTITY_TYPE_SPHERE;
					}

					//円盤
					if (strncmp(buf, "disk ", 5) == 0)
					{
						rr = 0.0;
						sscanf4fs(buf, "disk %lf %lf %lf %lf", &x, &y, &z, &r);
						getLine(buf, 1024, fp);
						sscanf3fs(buf, "normal %lf %lf %lf", &normal.x, &normal.y, &normal.z);
						type = ENTITY_TYPE_CIRCLE;
					}
					if (strncmp(buf, "disk2 ", 6) == 0)
					{
						sscanf5fs(buf, "disk2 %lf %lf %lf %lf %lf", &x, &y, &z, &r, &rr);
						getLine(buf, 1024, fp);
						sscanf3fs(buf, "normal %lf %lf %lf", &normal.x, &normal.y, &normal.z);
						type = ENTITY_TYPE_CIRCLE;
					}


					//OBJメッシュ
					if ( strncmp(buf, "objfile ", 8) == 0 )
					{
						strcpy(file, buf+8);
						char* p = strchr(file, '\n');
						if ( p ) *p = '\0';
						type = ENTITY_TYPE_POLYGON;
					}

					//形状属性定義
					while(	fgets(buf, 1024, fp) != NULL )
					{
						if ( buf[0] == '\n' )
						{
							break;
						}

						if ( strncmp(buf, "media_ignore_front_back ", 24) == 0 )
						{
							double z;
							sscanf1fs(buf, "media_ignore_front_back %lf", &z);
							media_ignore_front_back = (int)z;
							if ( media_ignore_front_back ) fprintf(stderr, "表裏を無視しました.\n");
						}
						if ( strncmp(buf, "texture ", 8) == 0 )
						{
							char* p = strchr(buf, '\n');
							if ( p ) *p = '\0';
							char* filename =  buf+8;

							char filename2[1024];
							FILE* fp = fopen(filename, "r");
							if ( fp )
							{
								printf("texture=[%s]\n", filename);
								fclose(fp);
							}else
							{
								sprintf(  filename2, "%s%s%s", drive, dir, filename);
								filename = filename2;
								printf("texture=[%s]\n", filename);
								fp = fopen(filename, "r");
								if ( fp )
								{
									fclose(fp);
								}else
								{
									printf("ファイルが開けません[%s]\n", filename);
									return ;
								}
							}

							BitMap* bmp = new BitMap();
							bmp->Read(filename);
							if ( bmp->GetImage())
							{
								material.texture = bmp;
							}else
							{
								material.texture = 0;
								delete bmp;
							}
							continue;
						}		
						if ( strncmp(buf, "bump_texture ", 13) == 0 )
						{
							char* p = strchr(buf, '\n');
							if ( p ) *p = '\0';
							char* filename =  buf+13;

							char filename2[1024];
							FILE* fp = fopen(filename, "r");
							if ( fp )
							{
								printf("bump_texture=[%s]\n", filename);
								fclose(fp);
							}else
							{
								sprintf(  filename2, "%s%s%s", drive, dir, filename);
								filename = filename2;
								printf("bump_texture=[%s]\n", filename);
								fp = fopen(filename, "r");
								if ( fp )
								{
									fclose(fp);
								}else
								{
									printf("ファイルが開けません[%s]\n", filename);
									return ;
								}
							}

							BitMap* bmp = new BitMap();
							bmp->Read(filename);
							if ( bmp->GetImage() )
							{
								material.bump_texture = bmp;
							}else
							{
								material.bump_texture = 0;
								delete bmp;
							}
							continue;
						}		
						
						if ( strncmp(buf, "ibl_texture ", 12) == 0 )
						{
							char* p = strchr(buf, '\n');
							if ( p ) *p = '\0';
							char* filename =  buf+12;

							char filename2[1024];
							FILE* fp = fopen(filename, "r");
							if ( fp )
							{
								printf("ibl_texture=[%s]\n", filename);
								fclose(fp);
							}else
							{
								sprintf(  filename2, "%s%s%s", drive, dir, filename);
								filename = filename2;
								printf("ibl_texture=[%s]\n", filename);
								fp = fopen(filename, "r");
								if ( fp )
								{
									fclose(fp);
								}else
								{
									printf("ファイルが開けません[%s]\n", filename);
									return ;
								}
							}
							char drive[_MAX_DRIVE];	// ドライブ名
							char dir[_MAX_DIR];		// ディレクトリ名
							char fname[_MAX_FNAME];	// ファイル名
							char ext[_MAX_EXT];		// 拡張子
							_splitpath( filename, drive, dir, fname, ext );

							if ( strcmpi(ext, ".bmp") == 0 )
							{
								BitMap* bmp = new BitMap();
								bmp->Read(filename);
								if ( bmp->GetImage())
								{
									material.IBL_texture = bmp;
								}else
								{
									material.IBL_texture = 0;
									delete bmp;
								}
							}else
							if ( strcmpi(ext, ".hdr") == 0 )
							{
								HDRImage* hdr = new HDRImage();
								hdr->load_unsafe( std::string(filename));
								if ( hdr->width() > 0 && hdr->height() >= 0 )
								{
									material.IBL_texture_HDR = hdr;
								}else
								{
									material.IBL_texture_HDR = 0;
									delete hdr;
								}
							}else
							{
								printf("no support image file format\n");
							}
							continue;
						}	
						if ( strncmp(buf, "ibl_texture_coef ", 17) == 0 )
						{
							sscanf1fs(buf, "ibl_texture_coef %lf", &material.ibl_texture_coef);
							continue;
						}
						if ( strncmp(buf, "hemisphere_map ", 15) == 0 )
						{
							int dmy;
							sscanf(buf, "hemisphere_map %d", &dmy);
							material.hemisphere_map = dmy;
							continue;
						}
						if ( strncmp(buf, "angular_map ", 12) == 0 )
						{
							int dmy;
							sscanf(buf, "angular_map %d", &dmy);
							material.angular_map = dmy;
							continue;
						}
						if ( strncmp(buf, "panoramic_map ", 14) == 0 )
						{
							int dmy;
							sscanf(buf, "panoramic_map %d", &dmy);
							material.panoramic_map = dmy;
							continue;
						}
						if (strncmp(buf, "disk_map ", 9) == 0)
						{
							int dmy;
							sscanf(buf, "disk_map %d", &dmy);
							material.disk_map = dmy;
							continue;
						}
						if (strncmp(buf, "repeat ", 7) == 0)
						{
							sscanf(buf, "repeat %d", &material.repeat);
							continue;
						}


						if ( strncmp(buf, "blackhole_disk ", 15 ) == 0 )
						{
							sscanf(buf, "blackhole_disk %d", &blackhole_disk);
							continue;
						}
						if ( strncmp(buf, "normal ", 7 ) == 0 )
						{
							sscanf(buf, "normal %d", &normal_vector_inverse);
							continue;
						}
						if ( strncmp(buf, "background ", 11 ) == 0 )
						{
							int s;
							sscanf(buf, "background %d", &s);
							if ( s ) background = true;
							else background = false;
							continue;
						}
						if ( strncmp(buf, "smooth ", 7 ) == 0 )
						{
							sscanf(buf, "smooth %d", &smooth);
							continue;
						}
						if ( strncmp(buf, "shadow ", 7 ) == 0 )
						{
							sscanf(buf, "shadow %d", &shadow);
							continue;
						}
						if ( strncmp(buf, "back ", 5) == 0 )
						{
							sscanf(buf, "back %d", &back);
							continue;
						}
						if (strncmp(buf, "parallel_light ", 15) == 0)
						{
							sscanf(buf, "parallel_light %d", &parallel_light);
							printf("parallel_light %d\n", parallel_light);
							continue;
						}
						if (strncmp(buf, "parallel_light_dir ", 19) == 0)
						{
							sscanf3f(buf + 19, "%lf %lf %lf", &light_dir.x, &light_dir.y, &light_dir.z);
							printf("parallel_light_dir %f %f %f\n", light_dir.x, light_dir.y, light_dir.z);
							light_dir = normalize(light_dir);
							continue;
						}
						if (strncmp(buf, "spot_light ", 11) == 0)
						{
							sscanf(buf, "spot_light %d", &spot_light);
							printf("spot_light %d\n", spot_light);
							continue;
						}
						if (strncmp(buf, "spot_light_dir ", 15) == 0)
						{
							sscanf3f(buf + 15, "%lf %lf %lf", &light_dir.x, &light_dir.y, &light_dir.z);
							printf("spot_light_dir %f %f %f\n", light_dir.x, light_dir.y, light_dir.z);
							light_dir = normalize(light_dir);
							continue;
						}
						if (strncmp(buf, "spot_light_angle ", 17) == 0)
						{
							sscanf2f(buf + 17, "%lf %lf", spot_angle, spot_angle+1);
							printf("spot_light_angle %f %f\n", spot_angle[0], spot_angle[1]);
							continue;
						}
						if (strncmp(buf, "spot_light_falloff ", 19) == 0)
						{
							sscanf1f(buf + 19, "%lf", &falloff);
							printf("spot_light_falloff %f\n", falloff);
							continue;
						}
						if (strncmp(buf, "spot_attenuation ", 17) == 0)
						{
							sscanf3f(buf + 17, "%lf %lf %lf", attenuation, attenuation + 1, attenuation+2);
							printf("spot_attenuation %f %f %f\n", attenuation[0], attenuation[1], attenuation[2]);
							continue;
						}
						

						
						if ( strcmp(buf, "MATERIAL_SYMBOL\n") == 0 )
						{
							getLine(buf, 1024, fp);
							char* p = strchr(buf, '\n');
							if ( p ) *p = '\0';
							bool exist = false;
							for ( int ii = 0; true; ii++ )
							{
								if ( strcmp(buf, material_symbol_[ii].name) == 0 )
								{
									material.participatingMediaParam.scatteringColor = Color(material_symbol_[ii].scattering);
									material.participatingMediaParam.absorbingColor = Color(material_symbol_[ii].absorbing);
									material.participatingMediaParam.eta = material_symbol_[ii].eta;
									material.refractive_index = material_symbol_[ii].eta;
									material.r_refractive_index = material_symbol_[ii].eta;

									material.participatingMediaParam.phase_prm[0] = material_symbol_[ii].g[0];
									material.participatingMediaParam.phase_prm[1] = material_symbol_[ii].g[1];
									material.participatingMediaParam.phase_prm[2] = material_symbol_[ii].g[2];

									material.participatingMediaParam.scattering = (material.participatingMediaParam.scatteringColor.x + material.participatingMediaParam.scatteringColor.y + material.participatingMediaParam.scatteringColor.z) / 3.0;
									material.participatingMediaParam.absorbing = (material.participatingMediaParam.absorbingColor.x + material.participatingMediaParam.absorbingColor.y + material.participatingMediaParam.absorbingColor.z)/3.0;
									exist = true;
									break;
								}
							}
							if (!exist )
							{
								fprintf(stderr, "MATERIAL_SYMBOL [%s] ? not found!!\n", buf);
							}
							if ( use_brssdf )
							{
								material.participatingMediaParam.use_brssdf = (int)atof(use_brssdf);
								fprintf(stderr, "material.participatingMediaParam.use_brssdf %d\n", material.participatingMediaParam.use_brssdf);
							}
							continue;
						}
						if (strcmp(buf, "MATERIAL_R_MAX\n") == 0)
						{
							getLine(buf, 1024, fp);
							material.participatingMediaParam.r_max = atof(buf);
							continue;
						}
						if (strcmp(buf, "MATERIAL_SCALE\n") == 0)
						{
							getLine(buf, 1024, fp);
							material_scale = atof(buf);
							continue;
						}
						if ( strcmp(buf, "SCATTERING\n") == 0 )
						{
							getLine(buf, 1024, fp);
							material.participatingMediaParam.scattering = atof(buf);
							getLine(buf, 1024, fp);
							sscanf3f(buf, "%lf %lf %lf", &material.participatingMediaParam.scatteringColor.x, &material.participatingMediaParam.scatteringColor.y, &material.participatingMediaParam.scatteringColor.z);
							
							if ( use_brssdf )
							{
								material.participatingMediaParam.use_brssdf = (int)atof(use_brssdf);
								fprintf(stderr, "material.participatingMediaParam.use_brssdf %d\n", material.participatingMediaParam.use_brssdf);
							}
							continue;
						}
						if ( strcmp(buf, "ABSORBING\n") == 0 )
						{
							getLine(buf, 1024, fp);
							material.participatingMediaParam.absorbing = atof(buf);
							getLine(buf, 1024, fp);
							sscanf3f(buf, "%lf %lf %lf", &material.participatingMediaParam.absorbingColor.x, &material.participatingMediaParam.absorbingColor.y, &material.participatingMediaParam.absorbingColor.z);
							continue;
						}
						if ( strcmp(buf, "BSSRDF\n") == 0 )
						{
							getLine(buf, 1024, fp);
							material.participatingMediaParam.use_brssdf = (int)atof(buf);
							if ( use_brssdf )
							{
								material.participatingMediaParam.use_brssdf = (int)atof(use_brssdf);
								fprintf(stderr, "material.participatingMediaParam.use_brssdf %d\n", material.participatingMediaParam.use_brssdf);
							}
							getLine(buf, 1024, fp);
							continue;
						}

						if ( strcmp(buf, "PAHASE\n") == 0 )
						{
							getLine(buf, 1024, fp);
							material.participatingMediaParam.phase_prm[0] = atof(buf);
							material.participatingMediaParam.phase_prm[1] = material.participatingMediaParam.phase_prm[0];
							material.participatingMediaParam.phase_prm[2] = material.participatingMediaParam.phase_prm[0];
							continue;
						}
						if (strncmp(buf, "PAHASE ", 7) == 0)
						{
							sscanf3fs(buf, "PAHASE %lf %lf %lf", &x, &y, &z);
							material.participatingMediaParam.phase_prm[0] = x;
							material.participatingMediaParam.phase_prm[1] = x;
							material.participatingMediaParam.phase_prm[2] = z;
							continue;
						}

						if ( strcmp(buf, "level\n") == 0 )
						{
							getLine(buf, 1024, fp);
							material.participatingMediaParam.level = (int)atof(buf);
							
							if ( material.participatingMediaParam.level == -2 )  material.participatingMediaParam.use_brssdf = 1;
							if ( material.participatingMediaParam.level >= 0 )  material.participatingMediaParam.use_brssdf = 0;
							fprintf(stderr, "=>material.participatingMediaParam.use_brssdf %d\n", material.participatingMediaParam.use_brssdf);
							continue;
						}

						if ( strncmp(buf, "rotation ", 9) == 0 )
						{
							sscanf3fs(buf, "rotation %lf %lf %lf", &x, &y, &z);
							Matrix4x4 mat;
							matrixList.push_back(mat.RotationX(x*PS_RAD));
							matrixList.push_back(mat.RotationY(y*PS_RAD));
							matrixList.push_back(mat.RotationZ(z*PS_RAD));
							continue;
						}
						if ( strncmp(buf, "translate ", 10 ) == 0 )
						{
							sscanf3fs(buf, "translate %lf %lf %lf", &dx, &dy, &dz);
							Matrix4x4 mat;
							matrixList.push_back(mat.Translation(dx, dy, dz));
							continue;
						}
						if ( strncmp(buf, "scale ", 6 ) == 0 )
						{
							sscanf3fs(buf, "scale %lf %lf %lf", &sx, &sy, &sz);
							Matrix4x4 mat;
							matrixList.push_back(mat.Scale(sx, sy, sz));
							continue;
						}
						if ( strncmp(buf, "emission ", 9 ) == 0 )
						{
							sscanf3fs(buf, "emission %lf %lf %lf", &material.emission.x, &material.emission.y, &material.emission.z);
							continue;
						}
						if ( strncmp(buf, "color ", 6 ) == 0 )
						{
							sscanf3fs(buf, "color %lf %lf %lf", &material.color.x, &material.color.y, &material.color.z);
							continue;
						}
						if ( strncmp(buf, "specular ", 6 ) == 0 )
						{
							sscanf3fs(buf, "specular %lf %lf %lf", &material.specular.x, &material.specular.y, &material.specular.z);
							continue;
						}
						if ( strncmp(buf, "refractive_index ", 17) == 0 )
						{
							sscanf1fs(buf, "refractive_index %lf", &material.refractive_index);
							continue;
						}
						if ( strncmp(buf, "r_refractive_index ", 19) == 0 )
						{
							sscanf1fs(buf, "r_refractive_index %lf", &material.r_refractive_index);
							continue;
						}
						if ( strncmp(buf, "ward_brdf ", 10 ) == 0 )
						{
							double lo_s;
							double lo_d;

							sscanf4fs(buf, "ward_brdf %lf %lf %lf %lf", &lo_d, &lo_s, &material.brdfParameter.ward_brdf.alp_x, &material.brdfParameter.ward_brdf.alp_y);
							material.color = material.color * lo_d;
							material.specular = material.specular * lo_s;
							continue;
						}
						if (strncmp(buf, "phong_brdf ", 11) == 0)
						{
							double lo_s;
							double lo_d;

							sscanf3fs(buf, "phong_brdf %lf %lf %lf", &lo_d, &lo_s, &material.brdfParameter.phong_brdf.specular_exponent);
							material.color = material.color * lo_d;
							material.specular = material.specular * lo_s;
							continue;
						}
						if (strncmp(buf, "roughness ", 10) == 0)
						{
							sscanf1fs(buf, "roughness %lf", &material.roughness);
							continue;
						}
						

						if ( strcmp(buf, "reflection diffuse\n") == 0 )
						{
							material.reflection_type = REFLECTION_TYPE_DIFFUSE;
							material.roughness = 1.0;
							continue;
						}
						if ( strcmp(buf, "reflection reflection\n") == 0 )
						{
							material.reflection_type = REFLECTION_TYPE_REFRACTION;
							if ( material.refractive_index == 0.0 )
							{
								material.refractive_index = REFRACTIVE_INDEX_OF_GLASS;
								material.r_refractive_index = REFRACTIVE_INDEX_OF_GLASS/REFRACTIVE_INDEX_OF_AIR;
							}
							continue;
						}
						if ( strcmp(buf, "reflection reflection1\n") == 0 )
						{
							material.reflection_type = REFLECTION_TYPE_REFRACTION;
							material.formulaType = REFRACTIVE_FORMULA1;
							getLine(buf, 512, fp);
							sscanf(buf, "%lf %lf %lf %lf %lf %lf", 
								material.Dispersion_formula1, material.Dispersion_formula1+1, material.Dispersion_formula1+2,
								material.Dispersion_formula1+3, material.Dispersion_formula1+4, material.Dispersion_formula1+5);
							continue;
						}
						if ( strcmp(buf, "reflection reflection2\n") == 0 )
						{
							material.reflection_type = REFLECTION_TYPE_REFRACTION;
							material.formulaType = REFRACTIVE_FORMULA2;
							getLine(buf, 512, fp);
							sscanf(buf, "%lf %lf %lf %lf %lf %lf", 
								material.Dispersion_formula1, material.Dispersion_formula1+1, material.Dispersion_formula1+2,
								material.Dispersion_formula1+3, material.Dispersion_formula1+4, material.Dispersion_formula1+5);
							continue;
						}
						if ( strcmp(buf, "reflection reflection99\n") == 0 )
						{
							material.reflection_type = REFLECTION_TYPE_REFRACTION;
							material.formulaType = REFRACTIVE_FORMULA99;
							continue;
						}
						if (strcmp(buf, "reflection specular\n") == 0)
						{
							material.reflection_type = REFLECTION_TYPE_SPECULAR;
							material.roughness = 0.0;
							continue;
						}						
						if ( strcmp(buf, "reflection ward_brdf\n") == 0 )
						{
							material.reflection_type = REFLECTION_TYPE_WARD_BRDF;
							continue;
						}
						if (strcmp(buf, "reflection phong_brdf\n") == 0)
						{
							material.reflection_type = REFLECTION_TYPE_PHONG_BRDF;
							continue;
						}
						if (strcmp(buf, "reflection Subsurface_Scattering\n") == 0)
						{
							if (material.reflection_type == REFLECTION_TYPE_REFRACTION)
							{
								material.reflection_type = REFLECTION_TYPE_SSS_REFRACTION;
							}else
							if (material.reflection_type == REFLECTION_TYPE_WARD_BRDF)
							{
								material.reflection_type = REFLECTION_TYPE_SSS_WARD_BRDF;
							}else
							if ( material.reflection_type == REFLECTION_TYPE_DIFFUSE )
							{
								material.reflection_type = REFLECTION_TYPE_SSS_DIFFUSE;
							}else
							{
								material.reflection_type = REFLECTION_TYPE_SSS_DIFFUSE;
							}
							continue;
						}
						
					}

					//形状登録
					if ( type == ENTITY_TYPE_INF_LIGHT )
					{

						if ( material.isLight() )
						{
							Light lit;
							lit.light = 0;
							lit.light_id = light_list.list.size();
							lit.dir = normal;
							lit.infinity_light_emitssin = material.emission;
							lit.infinity_light = 1;
							lit.power = material.emission.length();
							light_list.list.push_back(lit);
							printf("infinity_light dir %f %f %f\n", lit.dir.x, lit.dir.y, lit.dir.z);
							fprintf(stderr, "infinity_light dir %f %f %f\n", lit.dir.x, lit.dir.y, lit.dir.z);
							//fprintf(stderr, "--------------\n");
						}
					}

					if ( type == ENTITY_TYPE_UVPLANE )
					{
						material.background = background;
						material_id = material_list.Add(material);

						UVPlane* plane = new UVPlane(org, U, V, material_id);
						plane->materialList = &material_list;
						plane->base_plane->materialList = plane->materialList;

						plane->normal_vector_inverse = normal_vector_inverse;
						plane->MatrixTransformation(matrixList);
						plane->back = back;
						plane->CalcArea();
						plane->rnd = &rand_;
						plane->shadow = (shadow)? true: false;
						plane->media_ignore_front_back = (unsigned char)media_ignore_front_back;
						EntList.Add((Entity*)plane);

						if ( material.participatingMediaParam.scattering > PS_EPS )
						{
							plane->material()->participatingMediaParam.setup_Scattering(material_scale);
							plane->scattering = true;
						}
						if ( material.isLight() )
						{
							plane->light = true;
							Light lit;
							lit.light = (Entity*)plane;
							lit.light_id = light_list.list.size();
							lit.parallel_light = parallel_light;
							lit.dir = light_dir;
							if (lit.dir.length() < PS_EPS16)
							{
								lit.dir = plane->normal*plane->normal_vector_inverse;
							}
							plane->light_id = lit.light_id;
							lit.power = plane->area *  material.emission.length();
							if ( material.IBL())
							{
								ibl = lit.light_id;
								lit.power = 1.0;
							}
							light_list.list.push_back(lit);
							
							if (parallel_light ) printf("parallel_light dir %f %f %f\n", lit.dir.x, lit.dir.y, lit.dir.z);
							if (parallel_light) fprintf(stderr, "parallel_light dir %f %f %f\n", lit.dir.x, lit.dir.y, lit.dir.z);
						}
						if (material.IBL())
						{
							ibl_ent = plane;
						}
					}
					if ( type == ENTITY_TYPE_PLANE )
					{
						material.background = background;
						material_id = material_list.Add(material);

						Plane* plane = new Plane(org, normal, material_id);
						plane->materialList = &material_list;
						plane->normal_vector_inverse = normal_vector_inverse;
						plane->MatrixTransformation(matrixList);
						plane->back = back;
						plane->shadow = (shadow)? true: false;
						plane->CalcArea();
						if ( material.participatingMediaParam.scattering > PS_EPS )
						{
							plane->material()->participatingMediaParam.setup_Scattering(material_scale);
							plane->scattering = true;
						}

						plane->media_ignore_front_back = (unsigned char)media_ignore_front_back;
						plane->rnd = &rand_;
						EntList.Add((Entity*)plane);
					}

					if (type == ENTITY_TYPE_CIRCLE)
					{
						material.background = background;
						material_id = material_list.Add(material);

						Circle* disk = new Circle(Vector3d(x, y, z), normal, r, material_id);
						disk->radius_inner = rr;
						disk->materialList = &material_list;
						disk->blackhole_disk = blackhole_disk;

						disk->base->materialList = disk->materialList;
						disk->base->base_plane->materialList = disk->materialList;
						disk->normal_vector_inverse = normal_vector_inverse;
						disk->MatrixTransformation(matrixList);
						disk->back = back;
						disk->CalcArea();
						if ( material.participatingMediaParam.scattering > PS_EPS )
						{
							disk->material()->participatingMediaParam.setup_Scattering(material_scale);
							disk->scattering = true;
						}

						disk->media_ignore_front_back = (unsigned char)media_ignore_front_back;
						disk->rnd = &rand_;
						EntList.Add((Entity*)disk);

						if ( material.isLight() )
						{
							disk->light = true;
							Light lit;
							lit.light = (Entity*)disk;
							lit.light_id = light_list.list.size();
							lit.parallel_light = parallel_light;
							lit.spot_light = spot_light;
							lit.dir = light_dir;
							disk->light_id = lit.light_id;
							lit.power = disk->area *  material.emission.length();

							if ( material.IBL())
							{
								ibl = lit.light_id;
								lit.power = 1.0;
							}
							light_list.list.push_back(lit);

						}
						if (material.IBL())
						{
							ibl_ent = disk;
						}
					}

					if (type == ENTITY_TYPE_SPHERE)
					{
						material_id = material_list.Add(material);

						Sphere* sph = new Sphere(r, Vector3d(x, y, z), material_id);
						sph->materialList = &material_list;
						sph->normal_vector_inverse = normal_vector_inverse;
						sph->hemisphere = hemisphere;
						sph->MatrixTransformation(matrixList);
						sph->back = back;
						sph->CalcArea();
						if (material.participatingMediaParam.scattering > PS_EPS)
						{
							sph->material()->participatingMediaParam.setup_Scattering(material_scale);
							sph->scattering = true;
						}

						sph->media_ignore_front_back = (unsigned char)media_ignore_front_back;
						sph->rnd = &rand_;
						EntList.Add((Entity*)sph);

						if (material.isLight())
						{
							sph->light = true;
							Light lit;
							lit.light = (Entity*)sph;
							lit.light_id = light_list.list.size();
							lit.parallel_light = parallel_light;
							lit.spot_light = spot_light;
							lit.dir = light_dir;
							sph->light_id = lit.light_id;
							lit.power = sph->area *  material.emission.length();

							if (lit.spot_light)
							{
								lit.power = lit.SpotLightPower(material.emission.length());
								lit.spot_light_cos_angle[0] = cos(spot_angle[0] * PS_RAD);
								lit.spot_light_cos_angle[1] = cos((spot_angle[0] - spot_angle[1]) * PS_RAD);
								lit.falloff = falloff;
								lit.attenuation[0] = attenuation[0];
								lit.attenuation[1] = attenuation[1];
								lit.attenuation[2] = attenuation[2];
							}
							if (material.IBL())
							{
								ibl = lit.light_id;
								lit.power = 1.0;
							}
							light_list.list.push_back(lit);

						}
						if (material.IBL())
						{
							ibl_ent = sph;
						}
					}
					if (type == ENTITY_TYPE_POLYGON)
					{
						char objfile[512];

						strcpy(objfile, file);
						FILE* fp = fopen(objfile, "r");
						if ( fp )
						{
							printf("OBJFILE=[%s]\n", objfile);
							fclose(fp);
						}else
						{
							sprintf(  objfile, "%s%s%s", drive, dir, file);
							printf("OBJFILE=[%s]\n", objfile);
							fp = fopen(objfile, "r");
							if ( fp )
							{
								fclose(fp);
							}else
							{
								printf("ファイルが開けません[%s]\n", objfile);
								return ;
							}
						}
						material_id = material_list.Add(material);

						Obj *obj = new Obj( objfile);
						Polygon* poly = new Polygon(obj, material_id);
						poly->materialList = &material_list;

						poly->normal_vector_inverse = normal_vector_inverse;
						
						poly->MatrixTransformation(matrixList);
						poly->LoadTexture();
						poly->LoadBumpTexture();
						obj->dump();

						poly->mesh->smooth = smooth;
						poly->CalcArea();//MakeTrianglesEntityで計算するので何もしない
						if ( material.participatingMediaParam.scattering > PS_EPS )
						{
							poly->material()->participatingMediaParam.setup_Scattering(material_scale);
							poly->scattering = true;
						}
						poly->back = back;

						poly->media_ignore_front_back = (unsigned char)media_ignore_front_back;
						poly->rnd = &rand_;


						int lit_idx = -1;
						if (material.isLight())
						{
							poly->light = true;
							Light lit;
							lit.light = (Entity*)poly;
							lit.light_id = light_list.list.size();
							lit.parallel_light = parallel_light;
							lit.spot_light = spot_light;
							lit.dir = light_dir;
							poly->light_id = lit.light_id;

							//areaはまだ計算されていない。
							//lit.power = poly->area *  material.emission.length();
							//fprintf(stderr, "poly->area %f\n", poly->area);

							if (lit.spot_light)
							{
								lit.power = lit.SpotLightPower(material.emission.length());
								lit.spot_light_cos_angle[0] = cos(spot_angle[0] * PS_RAD);
								lit.spot_light_cos_angle[1] = cos((spot_angle[0] - spot_angle[1]) * PS_RAD);
								lit.falloff = falloff;
								lit.attenuation[0] = attenuation[0];
								lit.attenuation[1] = attenuation[1];
								lit.attenuation[2] = attenuation[2];
							}
							if (material.IBL())
							{
								ibl = lit.light_id;
								lit.power = 1.0;
							}

							//powerを設定していないので覚えておく（areaが未計算なので)
							light_list.list.push_back(lit);
							lit_idx = light_list.list.size() - 1;

						}

						if (BVH_USE)
						{
							try{
								printf("Build BVH start\n");
								EntList.Add((Entity*)poly->ConstructBVH());
								printf("end.\n");
								continue;
							}catch(...)
							{}
						}
						if ( QBVH_USE )
						{
							try{
								printf("Build QBVH start\n");
								EntList.Add((Entity*)poly->ConstructQBVH());
								printf("end.\n");
							}catch(...)
							{
								printf("re tary ==> Build BVH start\n");
								EntList.Add((Entity*)poly->ConstructBVH());
								printf("end.\n");
							}
							continue;
						}

						
						{
							poly->MakeTrianglesEntity();
							EntList.Add((Entity*)poly);

							if (lit_idx >= 0)
							{
								light_list.list[lit_idx].power = poly->area *  material.emission.length();
							}
						}
					}
					continue;
				}
			
				if (buf[0] != '#' && buf[0] != '\n' && buf[0] != ' ' && buf[0] != '\t')
				{
					fprintf(stderr, "?? %s\n", buf);
					printf("?? %s\n", buf);
				}
			 }
			fclose(fp);
		}

		if ( target_position.x < PS_INF )
		{
			printf("=>注視点からカメラまでの距離と方向を算出\n");
			camera_dir = normalize( target_position - camera_position);
		}
		if ( camera_matlist.size() )
		{
			Vector3d org(0,0,0);

			printf("=>カメラの座標変換\n");
			for ( int i = 0; i < camera_matlist.size() ; i++ )
			{
				org = camera_matlist[i]*org;
				camera_dir = camera_matlist[i]*(camera_dir - org);
				camera_position = camera_matlist[i]*camera_position;
			}
			camera_dir = normalize(camera_dir);
		}
		if ( camera_angle < PS_INF )
		{
			printf("=>水平方向の視角からスクリーン(H,W)算出\n");
			world_screen_width = 2.0*world_screen_dist*tan( camera_angle*PS_RAD / 2.0 );
			world_screen_height = world_screen_width;
		}

		//複数の光源から確率的に一つを選択する準備
		light_list.random_light_setup();

		//σt =σs + σa 減衰係数(attenuation coefficient)
		participatingMediaParam.setup_Scattering();
		if (participatingMediaParam.transmittance < PS_EPS12)
		{
			participatingMedia = 0;
		}

		char* env = 0;
		if ( (env =getenv("SAMPLING")) ) samples = atoi(env);
		if ( (env =getenv("SUPERSAMPLING")) ) supersamples = atoi(env);
		if ( (env =getenv("NEXTEVENTESTIMATION")) ) nextEventEstimation = atoi(env);
		if ( (env =getenv("IMAGE_X")) ) image_width = atoi(env);
		if ( (env =getenv("IMAGE_Y")) ) image_height = atoi(env);
		if ( (env = getenv("TIMELIMIT"))) timeLimit= atof(env);
		if ( (env = getenv("IMAGEDUMP"))) imageDumpTime= atoi(env);
		if ((env = getenv("LUMINANCE_CUTOFF"))) luminance_cutoff = atof(env);

		if ((env = getenv("METROPOLIS"))) metropolisTransport = atoi(env);
		if ((env = getenv("ERPT"))) energyRedistributionPathTracing = atoi(env);
		if ((env = getenv("MUTATION")))	mutation = atof(env);
		if ((env = getenv("PRESAMPLE"))) pre_sample = atof(env);

		if ((env = getenv("ENVLIGHT_R"))) Environment_Light.x = atof(env);
		if ((env = getenv("ENVLIGHT_G"))) Environment_Light.y = atof(env);
		if ((env = getenv("ENVLIGHT_B"))) Environment_Light.z = atof(env);

		if ((env = getenv("OMP_NUM_THREADS"))) threads = (int)atof(env);
		printf("use OMP_NUM_THREADS->threads=%d\n", threads);
		
		camera_dir = normalize(camera_dir);
		printf("\nOUTPUT [%s]\n", output);
		printf("IMAGE %d x %d\n", image_width, image_height);
		printf("SAMPLING %d\n", samples);
		printf("SUPERSAMPLING %d\n", supersamples);
		printf("SCREEN %f x %f dist %f\n", world_screen_width, world_screen_height, world_screen_dist);
		printf("CAMERA_POS %f %f %f\n", camera_position.x, camera_position.y, camera_position.z);
		printf("CAMERA_DIR %f %f %f\n", camera_dir.x, camera_dir.y, camera_dir.z);

		printf("  CAMERA_UPVEC %f %f %f\n", camera_up.x, camera_up.y, camera_up.z);
		printf("gamma_offset:%f\n", gamma_offset);
		//const Vector3d screen_x = normalize(cross(camera_dir, camera_up));
		//camera_up = normalize(cross(screen_x, camera_dir));

		//printf("=>CAMERA_UPVEC %f %f %f\n", camera_up.x, camera_up.y, camera_up.z);

		if ( camera_angle < PS_INF )printf("CAMERA_ANGLE %f\n", camera_angle);
		else printf("*CAMERA_ANGLE %f\n", atan((world_screen_width / 2.0) / world_screen_dist)*2.0 / PS_RAD);

		if ( target_position.x < PS_INF )printf("TARGET_POS %f %f %f\n", target_position.x, target_position.y, target_position.z);

		printf("BACKGROUND %f %f %f\n", BackgroundColor.x, BackgroundColor.y, BackgroundColor.z);

		if (metropolisTransport)
		{
			if (pre_sample <= 0.0) pre_sample = 0.4;
			if (mutation <= 0.0) mutation = 32;
			if (mlt_r1 == 0.0 || mlt_r2 == 0.0)
			{
				mlt_r1 = 0.1;
				mlt_r2 = 0.2;
			}
			//nextEventEstimation = 1;
			useTentFilter = 1;
			sensor_response = 0.25*sensor_response;
		}
		if (energyRedistributionPathTracing)
		{
			if (pre_sample <= 0.0) pre_sample = 0.4;
			if (mutation <= 0.0) mutation = 32;
			if (mlt_r1 == 0.0 || mlt_r2 == 0.0)
			{
				mlt_r1 = 1.0;
				mlt_r2 = 0.25;
			}
			//nextEventEstimation = 1;
			useTentFilter = 1;
			sensor_response = 0.25*sensor_response;
		}
		if (metropolisTransport || energyRedistributionPathTracing)
		{
			printf("pre_sample %f\n", pre_sample);
			printf("mutation %f\n", mutation);
			printf("mlt_r1 %f  mlt_r2 %f\n", mlt_r1, mlt_r2);
			if (mlt_r1 >= mlt_r2)
			{
				printf("parameter error, mlt_r1(%f) < mlt_r2(%f)\n", mlt_r1, mlt_r2);
			}
		}

		printf("metropolisTransport %d\n", metropolisTransport);
		if (metropolisTransport)
		{
			printf("metropolisTransport mutation %.3f\n", mutation);
		}

		printf("energyRedistributionPathTracing %d\n", energyRedistributionPathTracing);
		if (energyRedistributionPathTracing)
		{
			printf("metropolisTransport mutation %.3f\n", mutation);
		}
		if (energyRedistributionPathTracing && metropolisTransport)
		{
			fprintf(stderr, "option error:\"energyRedistributionPathTracing\" and \"energyRedistributionPathTracing\" !!\n");
		}

		printf("nextEventEstimation %d\n", nextEventEstimation);
		printf("participatingMedia %d\n", participatingMedia);
		printf("THREAD %d\n", threads);

		printf("DEPTH_LIMIT %d\n", DepthLimit);
		printf("DEPTH %d\n", Depth);

		printf("LUMINANCE_CUTOFF:%f\n", luminance_cutoff);

		if (!blackHole_exist && !backgroundMap[0])
		{
			bool texture_exist = false;
			if (!backgroundTexture[0])
			{
				char filename2[1024];
				FILE* fp = fopen(background_texture[0].c_str(), "r");
				if (fp)
				{
					printf("texture=[%s]\n", background_texture[0].c_str());
					fclose(fp);
					texture_exist = true;
				}
				else
				{
					sprintf(filename2, "%s%s%s", drive, dir, background_texture[0].c_str());
					printf("texture=[%s]\n", filename2);
					background_texture[0] = filename2;
					fp = fopen(filename2, "r");
					if (fp)
					{
						fclose(fp);
						texture_exist = true;
					}
					else
					{
						printf("background_textureファイルが開けません[%s]\n", background_texture[0].c_str());
					}
				}

				if (texture_exist)
				{
					backgroundTexture[0] = new BitMap;
					backgroundTexture[0]->Read((char*)background_texture[0].c_str());
					if (!backgroundTexture[0]->GetImage())
					{
						fprintf(stderr, "backgroundTexture error.\n");
					}
					else
					{
						printf("backgroundTexture %d x %d\n", backgroundTexture[0]->W(), backgroundTexture[0]->H());
						fprintf(stderr, "backgroundTexture %d x %d\n", backgroundTexture[0]->W(), backgroundTexture[0]->H());
						backgroundMap[0] = new SphericalMapping(backgroundTexture[0]->W(), backgroundTexture[0]->H());
					}
				}
			}
		}

		if (blackHole_exist)
		{
			if (!blackHolePos_set)
			{
				double Z = blackhole_dist;
				if (getenv("BLACK_HOLE_Z"))
				{
					Z = atof(getenv("BLACK_HOLE_Z"));
				}
				blackHolePos = camera_position + camera_dir * Z;
			}
			blackHole = new KerrBlackHole(blackHolePos.x, blackHolePos.y, blackHolePos.z, accretion_disk, 1.0, angular_momentum, camera_position, camera_dir);
			
			blackHole->initial_condition = initial_condition;
			blackHole->geodesics_max_length = geodesics_max_length;
			blackHole->accretion_disk_texture = accretion_disk_texture;
			blackHole->background_texture = background_texture[0];
			for (int ii = 0; ii < blackHole->accretion_disk_texture.size(); ii++)
			{
				printf("accretion_disk_texture[%d][%s]\n", ii, blackHole->accretion_disk_texture[ii].c_str());
			}
			printf("background_texture[%s]\n", blackHole->background_texture.c_str());
			blackHole->background_texture_map_coef[0] = background_texture_map_coef[0][0];
			blackHole->background_texture_map_coef[1] = background_texture_map_coef[0][1];
			blackHole->background_texture_map_coef[2] = background_texture_map_coef[0][2];
			blackHole->background_texture_map_coef[3] = background_texture_map_coef[0][3];
			blackHole->background_texture_coef = background_texture_coef[0];

			BlackHoleTextureInit(drive, dir, blackHole);
			fprintf(stderr, "BlackHole Position:%f %f %f\n", blackHolePos.x, blackHolePos.y, blackHolePos.z);
			fprintf(stderr, "事象の地平面半径:%f\n", blackHole->Rhor);
			fprintf(stderr, "内部 %f 膠着円盤 %f\n", blackHole->Rmstable, blackHole->Rdisk);
			fprintf(stderr, "angular_momentum %f\n", blackHole->a);
			fprintf(stderr, "initial_condition %d\n", blackHole->initial_condition);

			printf("BlackHole Position:%f %f %f\n", blackHolePos.x, blackHolePos.y, blackHolePos.z);
			printf("事象の地平面半径:%f\n", blackHole->Rhor);
			printf("内部 %f 膠着円盤 %f\n", blackHole->Rmstable, blackHole->Rdisk);
			printf("angular_momentum %f\n", blackHole->a);
			printf("initial_condition %d\n", blackHole->initial_condition);
		}

		if (wormHole_exist)
		{
			if (!wormHolePos_set)
			{
				double Z =wormhole_dist;
				if (getenv("WORM_HOLE_Z"))
				{
					Z = atof(getenv("WORM_HOLE_Z"));
				}
				wormholePos1 = camera_position + camera_dir * Z;
			}
			wormhole_ = new WormHole(
				wormholePos1.x, wormholePos1.y, wormholePos1.z, 
				wormholePos2.x, wormholePos2.y, wormholePos2.z,
				wormhole_Rho, wormhole_a, wormhole_W,
				camera_position, camera_dir);

			wormhole_->usr_set_sign_ = (int)wormhole_sign;
			wormhole_->camera_up = camera_up;
			wormhole_->screen_x = normalize(cross(camera_dir, camera_up));
			wormhole_->screen_y = normalize(cross(wormhole_->screen_x, camera_dir));

			wormhole_->geodesics_max_length = geodesics_max_length;
			wormhole_->background_texture[0] = background_texture[0];
			wormhole_->background_texture[1] = background_texture[1];
			for (int ii = 0; ii < 2; ii++)
			{
				printf("background_texture[%d][%s]\n", ii, wormhole_->background_texture[ii].c_str());
			}
			wormhole_->background_texture_map_coef[0][0] = background_texture_map_coef[0][0];
			wormhole_->background_texture_map_coef[0][1] = background_texture_map_coef[0][1];
			wormhole_->background_texture_map_coef[0][2] = background_texture_map_coef[0][2];
			wormhole_->background_texture_map_coef[0][3] = background_texture_map_coef[0][3];

			wormhole_->background_texture_map_coef[1][0] = background_texture_map_coef[1][0];
			wormhole_->background_texture_map_coef[1][1] = background_texture_map_coef[1][1];
			wormhole_->background_texture_map_coef[1][2] = background_texture_map_coef[1][2];
			wormhole_->background_texture_map_coef[1][3] = background_texture_map_coef[1][3];

			wormhole_->background_texture_coef[0] = background_texture_coef[0];
			wormhole_->background_texture_coef[1] = background_texture_coef[1];

			WormHoleTextureInit(drive, dir, wormhole_);
			fprintf(stderr, "WormkHole Position:%f %f %f -> %f %f %f\n", wormholePos1.x, wormholePos1.y, wormholePos1.z, wormholePos2.x, wormholePos2.y, wormholePos2.z);
			fprintf(stderr, "事象の地平面半径:%f\n", wormhole_->Rho);
			fprintf(stderr, "内部長さ %f\n", wormhole_->a);
			fprintf(stderr, "内部と接続する湾曲幅 %f\n", wormhole_->W);

			printf("WormkHole Position:%f %f %f -> %f %f %f\n", wormholePos1.x, wormholePos1.y, wormholePos1.z, wormholePos2.x, wormholePos2.y, wormholePos2.z);
			printf("事象の地平面半径:%f\n", wormhole_->Rho);
			printf("内部長さ %f\n", wormhole_->a);
			printf("内部と接続する湾曲幅 %f\n", wormhole_->W);
		}

#if USE_PARTICIPATING_MEDIA
		printf("散乱(scattering %f) %f %f %f\n", participatingMediaParam.scattering, participatingMediaParam.scatteringColor.x, participatingMediaParam.scatteringColor.y, participatingMediaParam.scatteringColor.z);
		printf("吸収(absorption %f) %f %f %f\n", participatingMediaParam.absorbing, participatingMediaParam.absorbingColor.x, participatingMediaParam.absorbingColor.y, participatingMediaParam.absorbingColor.z);
		printf("減衰係数(attenuation coefficient %f) %f %f %f\n", participatingMediaParam.transmittance, participatingMediaParam.transmittanceColor.x, participatingMediaParam.transmittanceColor.y, participatingMediaParam.transmittanceColor.z);
		printf("位相関数パラメータ[Henyey-Greenstein関数] %f %f %f\n", participatingMediaParam.phase_prm[0], participatingMediaParam.phase_prm[1], participatingMediaParam.phase_prm[2]);
#endif
	}
};

};

#endif
