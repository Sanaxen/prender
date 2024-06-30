#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include <math.h>
#include "ObjLoader.h"

#define BUFFER_OFFSET(bytes) ((GLubyte *)NULL + (bytes))

void Obj::realloc_face(int f)
{
	int cur_f = nf;
	nf++;

	try {
		fnorm.resize(nf);
		face.resize(nf);
		norm_id.resize(nf);
		texture_id.resize(nf);
		texture_material_id.resize(nf);

	}
	catch (std::bad_alloc e) {
		std::cerr << "メモリが足りません" << std::endl;
		abort();
	}

	for (int i = cur_f; i < nf; i++)
	{
		texture_material_id[i] = -1;
	}
}

/*
** コンストラクタ
*/
Obj::Obj(char *name, int* rgb)
{
	max[0] = -9999999.0f;
	max[1] = -9999999.0f;
	max[2] = -9999999.0f;
	min[0] = 99999999.0f;
	min[1] = 99999999.0f;
	min[2] = 99999999.0f;
	daiag = 1.0;

  std::ifstream file(name, std::ios::in | std::ios::binary);
  char buf[1024];
  char buf0[1024];
  int i, v, f;
  int uvNum;
  int vnNum;

  if (!file) {
    std::cerr << name << " が開けません" << std::endl;
    return;
  }
  std::cout << "Success : 指定されたOBJファイル\n";
  std::cout << "File Name : " << materialName.c_str() << std::endl;

  printf("Load [%s]..", name);
  v = f = nvn = ntn = 0;
 int color_id = 0;
  /* データの数を調べる */
  while (file.getline(buf0, sizeof buf)) {
	{
		char *p = buf0;
		while(*p == ' ' || *p == '\t' ) p++;
		strcpy(buf, buf0);
	}
    if (buf[0] == 'v' && buf[1] == ' ') {
      ++v;
    }
    if (buf[0] == 'v' && buf[1] == 'n') {
      ++nvn;
    }
    if (buf[0] == 'v' && buf[1] == 't') {
      ++ntn;
    }
    else if (buf[0] == 'f' && buf[1] == ' ') {
      ++f;
    }
  }

  nv = v;
  nf = f;

  try {
    vert.resize(v);
	norm.resize(std::max(v,nvn));
	color.resize(v);
	fnorm.resize(f);
	face.resize(f);
	norm_id.resize(f);
	uv.resize(ntn);
	texture_id.resize(f);
	texture_material_id.resize(f);

  }
  catch (std::bad_alloc e) {
    std::cerr << "メモリが足りません" << std::endl;
    abort();
  }

  for ( int i = 0; i < f; i++ )
  {
	  texture_material_id[i] = -1;
  }

  /* ファイルの巻き戻し */
  file.clear();
  file.seekg(0L, std::ios::beg);

  double umin, umax, vmin,vmax;

  umin = 99999.0;
  umax = -99999.0;
  vmin = 99999.0;
  vmax = -99999.0;

  int current_texture_material_id = -1;

#if 10
	/* データの読み込み */
	v = f = 0;
	uvNum = 0;
	vnNum = 0;
	while (file.getline(buf0, sizeof buf)) 
	{
		{
			char *p = buf0;
			while(*p == ' ' || *p == '\t' ) p++;
			strcpy(buf, p);
		}
		//マテリアル定義ファイル名
		if (strncmp(buf, "mtllib ", 7) == 0) {
			char tmp[128];
			strcpy(tmp, buf+7);
			char* p = strchr(tmp, '\n');
			if ( p ) *p = '\0';
			p = strchr(tmp, '\r');
			if ( p ) *p = '\0';
			materialName = tmp;
			if ( !LoadMTLFile())
			{
				char drive[_MAX_DRIVE];	// ドライブ名
				char dir[_MAX_DIR];		// ディレクトリ名
				char fname[_MAX_FNAME];	// ファイル名
				char ext[_MAX_EXT];		// 拡張子
				_splitpath( name, drive, dir, fname, ext );
				materialName = drive;
				materialName += dir;
				materialName += tmp;
				LoadMTLFile();
			}
			continue;
		}else
		if (strncmp(buf, "usemtl ", 7) == 0) {
			//マテリアル定義ファイルに書かれている識別名
			char tmp[128];
			strcpy(tmp, buf+7);
			char* p = strchr(tmp, '\n');
			if ( p ) *p = '\0';
			p = strchr(tmp, '\r');
			if ( p ) *p = '\0';

			//　マテリアル名から検索
			current_texture_material_id = -1;
			for ( int i=0; i<texture_material.size(); i++ )
			{
				//　名前が一致したらマテリアル番号を格納
				if ( strcmpi(texture_material[i].name.c_str(), tmp) == 0 )
				{
					printf("texture_material[%d]=>%s\n", i, tmp);
					current_texture_material_id = i;
				}
			}
			if ( current_texture_material_id < 0 )
			{
				printf("名前が一致したマテリアルが見つからない.\n");
			}
			continue;
		}else
		if (buf[0] == 'v' && buf[1] == 'n') {
			//fprintf(stderr, "%d %s\n\n", vnNum, buf);
			sscanf(buf, "%*s %f %f %f", norm[vnNum].g, norm[vnNum].g + 1, norm[vnNum].g + 2);
			++vnNum;
		}else
		if (buf[0] == 'v' && buf[1] == 't') {
			sscanf(buf, "%*s %f %f", uv[uvNum].g, uv[uvNum].g + 1);
			if (umin > uv[uvNum].g[0]) umin = uv[uvNum].g[0];
			if (vmin > uv[uvNum].g[1]) vmin = uv[uvNum].g[1];
			if (umax < uv[uvNum].g[0]) umax = uv[uvNum].g[0];
			if (vmax < uv[uvNum].g[1]) vmax = uv[uvNum].g[1];
			++uvNum;
		}else
		if (buf[0] == 'v' && buf[1] == ' ') {
			int n = 0;
			float c[4];
			c[3] = 1;
			set_color = true;
			n = sscanf(buf, "%*s %f %f %f %f %f %f %d", vert[v].g, vert[v].g + 1, vert[v].g + 2, c, c + 1, c + 2, &color_id);
			if ( n != 7 )
			{
				color_id = -1;
				n = sscanf(buf, "%*s %f %f %f %f %f %f", vert[v].g, vert[v].g + 1, vert[v].g + 2, c, c + 1, c + 2);
				if ( n != 6 )
				{
					set_color = false;
					sscanf(buf, "%*s %f %f %f", vert[v].g, vert[v].g + 1, vert[v].g + 2);
				}
			}

			if ( set_color )
			{
				color[v].g[0] = c[0];
				color[v].g[1] = c[1];
				color[v].g[2] = c[2];
			}else
			{
				if ( rgb )
				{
					color[v].g[0] = (double)rgb[0];
					color[v].g[1] = (double)rgb[1];
					color[v].g[2] = (double)rgb[2];
				}
			}
			++v;
		} else if (buf[0] == 'f' && buf[1] == ' ')
		{
			int dmy[12];
#if 0
			if (sscanf(buf + 2, "%d/%d/%d %d/%d/%d %d/%d/%d", 
				face[f].g + 0, texture_id[f].g+0, norm_id[f]+0,
				face[f].g + 1, texture_id[f].g+1, norm_id[f]+1,
				face[f].g + 2, texture_id[f].g+2, norm_id[f]+2) == 9) 
			{
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--texture_id[f].g[0];
				--texture_id[f].g[1];
				--texture_id[f].g[2];
				--norm_id[f][0];
				--norm_id[f][1];
				--norm_id[f][2];
				++f;
				continue;
			}
			if (sscanf(buf + 2, "%d//%d %d//%d %d//%d", 
				face[f].g + 0, norm_id[f]+0,
				face[f].g + 1, norm_id[f]+1,
				face[f].g + 2, norm_id[f]+2) == 6) 
			{
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--norm_id[f][0];
				--norm_id[f][1];
				--norm_id[f][2];
				++f;
				continue;
			}
			if (sscanf(buf + 2, "%d/%d %d/%d %d/%d", 
				face[f].g + 0, texture_id[f].g+0,
				face[f].g + 1, texture_id[f].g+1,
				face[f].g + 2, texture_id[f].g+2) == 6) 
			{
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--texture_id[f].g[0];
				--texture_id[f].g[1];
				--texture_id[f].g[2];
				norm_id[f][0] = face[f].g[0];
				norm_id[f][1] = face[f].g[1];
				norm_id[f][2] = face[f].g[2];
				++f;
				continue;
			}
			if (sscanf(buf + 2, "%d %d %d", 
				face[f].g + 0,
				face[f].g + 1,
				face[f].g + 2) == 3) 
			{
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				norm_id[f][0] = face[f].g[0];
				norm_id[f][1] = face[f].g[1];
				norm_id[f][2] = face[f].g[2];
				++f;
				continue;
			}
#else
			if (sscanf(buf + 2, "%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",
				face[f].g + 0, texture_id[f].g + 0, norm_id[f].g + 0,
				face[f].g + 1, texture_id[f].g + 1, norm_id[f].g + 1,
				face[f].g + 2, texture_id[f].g + 2, norm_id[f].g + 2, dmy, dmy+1, dmy+2) == 12)
			{
				texture_material_id[f] = current_texture_material_id;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--texture_id[f].g[0];
				--texture_id[f].g[1];
				--texture_id[f].g[2];
				--norm_id[f].g[0];
				--norm_id[f].g[1];
				--norm_id[f].g[2];
				++f;

				realloc_face(f);
				texture_material_id[f] = current_texture_material_id;
				face[f].g[0] = face[f - 1].g[2];
				face[f].g[1] = dmy[0] - 1;
				face[f].g[2] = face[f - 1].g[0];
				texture_id[f].g[0] = texture_id[f - 1].g[2];
				texture_id[f].g[1] = dmy[1] - 1;
				texture_id[f].g[2] = texture_id[f - 1].g[0];
				norm_id[f].g[0] = norm_id[f - 1].g[2];
				norm_id[f].g[1] = dmy[2] - 1;
				norm_id[f].g[2] = norm_id[f - 1].g[0];
				++f;

				continue;
			}
			if (sscanf(buf + 2, "%d/%d/%d %d/%d/%d %d/%d/%d",
				face[f].g + 0, texture_id[f].g + 0, norm_id[f].g+0,
				face[f].g + 1, texture_id[f].g + 1, norm_id[f].g+1,
				face[f].g + 2, texture_id[f].g + 2, norm_id[f].g+2) == 9)
			{
				texture_material_id[f] = current_texture_material_id;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--texture_id[f].g[0];
				--texture_id[f].g[1];
				--texture_id[f].g[2];
				--norm_id[f].g[0];
				--norm_id[f].g[1];
				--norm_id[f].g[2];
				++f;
				continue;
			}

			if (sscanf(buf + 2, "%d//%d %d//%d %d//%d %d//%d",
				face[f].g + 0, norm_id[f].g+0,
				face[f].g + 1, norm_id[f].g+1,
				face[f].g + 2, norm_id[f].g + 2, dmy, dmy+1) == 8)
			{
				texture_material_id[f] = current_texture_material_id;
				//texture_material_id[f] = -1;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--norm_id[f].g[0];
				--norm_id[f].g[1];
				--norm_id[f].g[2];
				++f;

				realloc_face(f);
				texture_material_id[f] = current_texture_material_id;
				//texture_material_id[f] = -1;
				face[f].g[0] = face[f - 1].g[2];
				face[f].g[1] = dmy[0] - 1;
				face[f].g[2] = face[f - 1].g[0];
				norm_id[f].g[0] = norm_id[f - 1].g[2];
				norm_id[f].g[1] = dmy[1] - 1;
				norm_id[f].g[2] = norm_id[f - 1].g[0];
				++f;
				continue;
			}
			if (sscanf(buf + 2, "%d//%d %d//%d %d//%d",
				face[f].g + 0, norm_id[f].g+0,
				face[f].g + 1, norm_id[f].g+1,
				face[f].g + 2, norm_id[f].g+2) == 6)
			{
				texture_material_id[f] = current_texture_material_id;
				//texture_material_id[f] = -1;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--norm_id[f].g[0];
				--norm_id[f].g[1];
				--norm_id[f].g[2];
				++f;
				continue;
			}

			if (sscanf(buf + 2, "%d/%d %d/%d %d/%d %d/%d",
				face[f].g + 0, texture_id[f].g + 0,
				face[f].g + 1, texture_id[f].g + 1,
				face[f].g + 2, texture_id[f].g + 2,
				dmy, dmy+1) == 8)
			{
				texture_material_id[f] = current_texture_material_id;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--texture_id[f].g[0];
				--texture_id[f].g[1];
				--texture_id[f].g[2];
				++f;

				realloc_face(f);
				texture_material_id[f] = current_texture_material_id;
				face[f].g[0] = face[f - 1].g[2];
				face[f].g[1] = dmy[0] - 1;
				face[f].g[2] = face[f - 1].g[0];
				texture_id[f].g[0] = texture_id[f - 1].g[2];
				texture_id[f].g[1] = dmy[1] - 1;
				texture_id[f].g[2] = texture_id[f - 1].g[0];
				++f;
				continue;
			}
			if (sscanf(buf + 2, "%d/%d %d/%d %d/%d",
				face[f].g + 0, texture_id[f].g+0,
				face[f].g + 1, texture_id[f].g+1,
				face[f].g + 2, texture_id[f].g+2) == 6) 
			{
				texture_material_id[f] = current_texture_material_id;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				--texture_id[f].g[0];
				--texture_id[f].g[1];
				--texture_id[f].g[2];
				++f;
				continue;
			}


			if (sscanf(buf + 2, "%d %d %d %d",
				face[f].g + 0,
				face[f].g + 1,
				face[f].g + 2, dmy) == 4)
			{
				texture_material_id[f] = current_texture_material_id;
				//texture_material_id[f] = -1;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				++f;

				realloc_face(f);
				texture_material_id[f] = current_texture_material_id;
				//texture_material_id[f] = -1;
				face[f].g[0] = face[f - 1].g[2];
				face[f].g[1] = dmy[0] - 1;
				face[f].g[2] = face[f - 1].g[0];
				++f;
				continue;
			}
			if (sscanf(buf + 2, "%d %d %d",
				face[f].g + 0,
				face[f].g + 1,
				face[f].g + 2) == 3) 
			{
				texture_material_id[f] = current_texture_material_id;
				//texture_material_id[f] = -1;
				--face[f].g[0];
				--face[f].g[1];
				--face[f].g[2];
				++f;
				continue;
			}
#endif
		}
	}
#endif

	if ( uvNum == 0 )
	{
		uv.clear();
		texture_id.clear();
	}


	if ( vnNum )
	{
	  if(vnNum == nvn ) printf("vertex normal...OK\n");
	  else  printf("vertex normal...ERROR,\n");
	}
	
	if ( uvNum )
	{
	  if (uvNum == ntn) printf("vertex UV...OK\n");
	  else  printf("vertex UV...ERROR,\n");
	}
	
	CalcNormalVector();
  printf("...done\n");



  	printf("Create Obj(%x)\n", (int)((void*)this));
	fflush(stdout);
}


void Obj::CalcNormalVectorChk()
{
	if (nvn == 0) return;

	int i;

	Vector3df* fnorm_org = new Vector3df[nf];

  /* 面法線ベクトルの算出 */
#pragma omp parallel for
	for (i = 0; i < nf; ++i) 
	{
		float dx1 = vert[face[i].g[1]].g[0] - vert[face[i].g[0]].g[0];
		float dy1 = vert[face[i].g[1]].g[1] - vert[face[i].g[0]].g[1];
		float dz1 = vert[face[i].g[1]].g[2] - vert[face[i].g[0]].g[2];
		float dx2 = vert[face[i].g[2]].g[0] - vert[face[i].g[0]].g[0];
		float dy2 = vert[face[i].g[2]].g[1] - vert[face[i].g[0]].g[1];
		float dz2 = vert[face[i].g[2]].g[2] - vert[face[i].g[0]].g[2];

		fnorm_org[i].g[0] = dy1 * dz2 - dz1 * dy2;
		fnorm_org[i].g[1] = dz1 * dx2 - dx1 * dz2;
		fnorm_org[i].g[2] = dx1 * dy2 - dy1 * dx2;
	}
#pragma omp parallel for
	for (i = 0; i < nf; ++i)
	{
		float a = sqrt(fnorm_org[i].g[0] * fnorm_org[i].g[0] + fnorm_org[i].g[1] * fnorm_org[i].g[1] + fnorm_org[i].g[2] * fnorm_org[i].g[2]);

		if (a != 0.0) {
			fnorm_org[i].g[0] /= a;
			fnorm_org[i].g[1] /= a;
			fnorm_org[i].g[2] /= a;
		}
		//printf("%f %f %f\n", fnorm[i][0],fnorm[i][1],fnorm[i][2]);
	}

	if ( 1/*nvn == 0*/ )
	{
		nvn = nv;

#pragma omp parallel for
		for (i = 0; i < nf; ++i) 
		{
			norm_id[i].g[0] = face[i].g[0];
			norm_id[i].g[1] = face[i].g[1];
			norm_id[i].g[2] = face[i].g[2];

			float x = (norm[face[i].g[0]].g[0] + norm[face[i].g[1]].g[0] + norm[face[i].g[2]].g[0]) / 3.0;
			float y = (norm[face[i].g[0]].g[1] + norm[face[i].g[1]].g[1] + norm[face[i].g[2]].g[1]) / 3.0;
			float z = (norm[face[i].g[0]].g[2] + norm[face[i].g[1]].g[2] + norm[face[i].g[2]].g[2]) / 3.0;

			double ln = x*x + y*y + z*z;
			if (ln > 1.0e-16)
			{
				ln = sqrt(ln);
				x /= ln;
				y /= ln;
				z /= ln;

				if (fnorm_org[i].g[0] * x + fnorm_org[i].g[1] * y + fnorm_org[i].g[2] * z < -0.99)
				{
					unsigned int k = face[i].g[0];
					face[i].g[0] = face[i].g[2];
					face[i].g[2] = k;
					fnorm[i].g[0] = x;
					fnorm[i].g[1] = y;
					fnorm[i].g[2] = z;
				}
			}
		}
	}
	delete[] fnorm_org;
}

void Obj::CalcNormalVector()
{
	int i;

	/* 面法線ベクトルの算出 */
	if (nvn == 0)
	{
#pragma omp parallel for
		for (i = 0; i < nf; ++i)
		{
			norm_id[i].g[0] = face[i].g[0];
			norm_id[i].g[1] = face[i].g[1];
			norm_id[i].g[2] = face[i].g[2];


			float dx1 = vert[face[i].g[1]].g[0] - vert[face[i].g[0]].g[0];
			float dy1 = vert[face[i].g[1]].g[1] - vert[face[i].g[0]].g[1];
			float dz1 = vert[face[i].g[1]].g[2] - vert[face[i].g[0]].g[2];
			float dx2 = vert[face[i].g[2]].g[0] - vert[face[i].g[0]].g[0];
			float dy2 = vert[face[i].g[2]].g[1] - vert[face[i].g[0]].g[1];
			float dz2 = vert[face[i].g[2]].g[2] - vert[face[i].g[0]].g[2];

			fnorm[i].g[0] = dy1 * dz2 - dz1 * dy2;
			fnorm[i].g[1] = dz1 * dx2 - dx1 * dz2;
			fnorm[i].g[2] = dx1 * dy2 - dy1 * dx2;
		}
	}
	else
	{
#pragma omp parallel for
		for (i = 0; i < nf; ++i)
		{
			float x = (norm[norm_id[i].g[0]].g[0] + norm[norm_id[i].g[1]].g[0] + norm[norm_id[i].g[2]].g[0]) / 3.0;
			float y = (norm[norm_id[i].g[0]].g[1] + norm[norm_id[i].g[1]].g[1] + norm[norm_id[i].g[2]].g[1]) / 3.0;
			float z = (norm[norm_id[i].g[0]].g[2] + norm[norm_id[i].g[1]].g[2] + norm[norm_id[i].g[2]].g[2]) / 3.0;

			fnorm[i].g[0] = x;
			fnorm[i].g[1] = y;
			fnorm[i].g[2] = z;

			double ln = x*x + y*y + z*z;
			if (ln > 1.0e-16)
			{
				ln = sqrt(ln);
				x /= ln;
				y /= ln;
				z /= ln;
				fnorm[i].g[0] = x;
				fnorm[i].g[1] = y;
				fnorm[i].g[2] = z;
			}
		}
	}

	if (nvn == 0)
	{
		//delete [] norm;
		//delete [] norm_id;
		//norm = new Vector3df[nv];
		//norm_id = new idx[nv];
		nvn = nv;

		/* 頂点の仮想法線ベクトルの算出 */
		for (i = 0; i < nv; ++i)
		{
			norm[i].g[0] = norm[i].g[1] = norm[i].g[2] = 0.0;
		}

		int *count = new int[nv];
		memset(count, '\0', nv*sizeof(int));

#pragma omp parallel for
		for (i = 0; i < nf; ++i)
		{
			norm_id[i].g[0] = face[i].g[0];
			norm_id[i].g[1] = face[i].g[1];
			norm_id[i].g[2] = face[i].g[2];

			norm[face[i].g[0]].g[0] += fnorm[i].g[0];
			norm[face[i].g[0]].g[1] += fnorm[i].g[1];
			norm[face[i].g[0]].g[2] += fnorm[i].g[2];
			count[face[i].g[0]]++;

			norm[face[i].g[1]].g[0] += fnorm[i].g[0];
			norm[face[i].g[1]].g[1] += fnorm[i].g[1];
			norm[face[i].g[1]].g[2] += fnorm[i].g[2];
			count[face[i].g[1]]++;

			norm[face[i].g[2]].g[0] += fnorm[i].g[0];
			norm[face[i].g[2]].g[1] += fnorm[i].g[1];
			norm[face[i].g[2]].g[2] += fnorm[i].g[2];
			count[face[i].g[2]]++;
		}

		/* 法線ベクトルの正規化 */
#pragma omp parallel for
		for (i = 0; i < nv; ++i)
		{
			if (count[i] == 0) continue;
			norm[i].g[0] /= (double)count[i];
			norm[i].g[1] /= (double)count[i];
			norm[i].g[2] /= (double)count[i];
			float a = sqrt(norm[i].g[0] * norm[i].g[0] + norm[i].g[1] * norm[i].g[1] + norm[i].g[2] * norm[i].g[2]);

			if (a != 0.0) {
				norm[i].g[0] /= a;
				norm[i].g[1] /= a;
				norm[i].g[2] /= a;
			}
		}
		delete[] count;

#pragma omp parallel for
		for (i = 0; i < nf; ++i)
		{
			float a = sqrt(fnorm[i].g[0] * fnorm[i].g[0] + fnorm[i].g[1] * fnorm[i].g[1] + fnorm[i].g[2] * fnorm[i].g[2]);

			if (a != 0.0) {
				fnorm[i].g[0] /= a;
				fnorm[i].g[1] /= a;
				fnorm[i].g[2] /= a;
			}
			//printf("%f %f %f\n", fnorm[i][0],fnorm[i][1],fnorm[i][2]);
		}
	}
}

bool Obj::LoadMTLFile()
{
	std::ifstream file;
	char buf[128];
	char buf0[128];
	float tmp_float=0.0f;

	//　ファイルを開く
	file.open(materialName.c_str(), std::ios::in);
	if ( !file.is_open() )
	{
		std::cout << "Error : 指定されたMTLファイルが開けませんでした\n";
		std::cout << "File Name : " << materialName.c_str() << std::endl;
		return false;
	}
	std::cout << "Success : 指定されたMTLファイル\n";
	std::cout << "File Name : " << materialName.c_str() << std::endl;

	ObjMaterial* tmp_mat = 0;
	//　ファイルの末端までループ
	while ( !file.eof() )
	{

		//　1行読み取り
		file.getline(buf0, sizeof(buf));
		{
			char *p = buf0;
			while(*p == ' ' || *p == '\t' ) p++;
			strcpy(buf, p);
		}

		//　バッファの1文字目で判断
		switch ( buf[0] )
		{
		//　newmtl
		case 'n':
			texture_material.push_back(ObjMaterial());
			char tmp[128];
			{
				sscanf(buf, "newmtl %s", tmp);
				char* p = strchr(tmp, '\n');
				if ( p ) *p = '\0';
				p = strchr(tmp, '\r');
				if ( p ) *p = '\0';
			}
			tmp_mat = &texture_material[texture_material.size()-1];
			tmp_mat->name = tmp;
			tmp_mat->user_ReflectionType = -1;
			break;

		//　Ka, Kd, Ks
		case 'K':
			switch ( buf[1] )
			{
			//　Ambient
			case 'a':
				sscanf(buf, "Ka %f %f %f", tmp_mat->ambient.g, tmp_mat->ambient.g + 1, tmp_mat->ambient.g + 2);
				break;

			//　Diffuse
			case 'd':
				sscanf(buf, "Kd %f %f %f", tmp_mat->diffuse.g, tmp_mat->diffuse.g + 1, tmp_mat->diffuse.g + 2);
				break;

			//　Specular
			case 's':
				sscanf(buf, "Ks %f %f %f", tmp_mat->specular.g, tmp_mat->specular.g + 1, tmp_mat->specular.g + 2);
				break;
			}
			break;

		//　d
		case 'd':
			if ( sscanf(buf, "d %f", &tmp_float) == 1 )
			{
				tmp_mat->ambient.g[3] = tmp_float;
				tmp_mat->diffuse.g[3] = tmp_float;
				tmp_mat->specular.g[3] = tmp_float;
				tmp_mat->emission.g[3] = tmp_float;
			}
			break;

		//　Tr
		case 'T':
			if ( buf[1] == 'r' )
			{
				if ( sscanf(buf, "Tr %f", &tmp_float) == 1 )
				{
					tmp_mat->ambient.g[3] = std::max(0.0,1.0 - tmp_float);
					tmp_mat->diffuse.g[3] = std::max(0.0,1.0 - tmp_float);
					tmp_mat->specular.g[3] = std::max(0.0,1.0 - tmp_float);
					tmp_mat->emission.g[3] = std::max(0.0,1.0 - tmp_float);
				}
			}
			break;

		//　Ni
		case 'N':
			if ( buf[1] == 'i' )	//Refraction index
			{
				sscanf(buf, "Ni %f", &tmp_mat->refraction_index);
			}
			if (buf[1] == 's')	//optical_density 
			{
				sscanf(buf, "Ns %f", &tmp_mat->shininess);
			}
			break;
		case 'm':
			if ( strncmp(buf, "map_Kd ", 6) == 0 )
			{
				char tmp[128];
				strcpy(tmp, buf + 7);
				char* p = strchr(tmp, '\n');
				if ( p ) *p = '\0';
				p = strchr(tmp, '\r');
				if ( p ) *p = '\0';
				tmp_mat->texture_name = tmp;
				printf("texture=[%s]\n", tmp);
			}
			if ( strncmp(buf, "map_Bump ", 9) == 0 )
			{
				char tmp[128];
				strcpy(tmp, buf + 10);
				char* p = strchr(tmp, '\n');
				if ( p ) *p = '\0';
				p = strchr(tmp, '\r');
				if ( p ) *p = '\0';
				tmp_mat->bump_texture_name = tmp;
			}
			break;
		case '#':
			if (strncmp(buf, "#SMOOTH ", 8) == 0)
			{
				int dmy = 0;
				sscanf(buf, "#SMOOTH %d", &dmy);
				if ( dmy )
				{
					tmp_mat->smooth = 1;
				}else
				{
					tmp_mat->smooth = 2;
				}
			}
			if (strncmp(buf, "#REFLECTION ", 12) == 0)
			{
				sscanf(buf, "#REFLECTION %d", &tmp_mat->user_ReflectionType);
				printf("特殊処理(自前の反射属性)=>%s =>[%d]　name[%s]\n", buf, tmp_mat->user_ReflectionType, tmp_mat->name.c_str());
				if (tmp_mat->user_ReflectionType == 0)	tmp_mat->roughness = 0.0;
				if (tmp_mat->user_ReflectionType == 1)	tmp_mat->roughness = 0.0;
				if (tmp_mat->user_ReflectionType == 2 && tmp_mat->refraction_index == 0.0f)
				{
					printf("=== ERROR? 屈折率(Ni)未設定 => 1.5に設定\n");
					tmp_mat->refraction_index = 1.5f;
				}

			}
			if (strncmp(buf, "#EMISSION ", 10) == 0)
			{
				sscanf(buf, "#EMISSION %f %f %f", tmp_mat->emission.g, tmp_mat->emission.g + 1, tmp_mat->emission.g + 2);
				printf("特殊処理(自前の反射属性)=>%s\n", buf);
			}
			if (strncmp(buf, "#COLOR ", 7) == 0)
			{
				sscanf(buf, "#COLOR %f %f %f", tmp_mat->diffuse.g, tmp_mat->diffuse.g + 1, tmp_mat->diffuse.g + 2);
				printf("特殊処理(自前の反射属性)=>%s\n", buf);
			}
			if (strncmp(buf, "#WARD ", 6) == 0)
			{
				sscanf(buf, "#WARD %f %f", tmp_mat->ward, tmp_mat->ward + 1);
				printf("特殊処理(自前の反射属性)=>%s\n", buf);
			}
			if (strncmp(buf, "#SPECULAR ", 10) == 0)
			{
				sscanf(buf, "#SPECULAR %f %f %f", tmp_mat->specular.g, tmp_mat->specular.g + 1, tmp_mat->specular.g + 2);
				printf("特殊処理(自前の反射属性)=>%s\n", buf);
			}
			if (strncmp(buf, "#ROUGHNESS ", 11) == 0)
			{
				sscanf(buf, "#ROUGHNESS %f", &tmp_mat->roughness);
				printf("特殊処理(自前の反射属性)=>%s\n", buf);
			}
			if (strncmp(buf, "#REPEAT ", 8) == 0)
			{
				sscanf(buf, "#REPEAT %d", &tmp_mat->texture_repeat);
			}

			break;

		//　該当なし
		default:
			break;
		}
	}
	
	//　ファイルを閉じる
	file.close();

	for ( int i = 0; i < texture_material.size(); i++ )
	{
		if ( texture_material[i].user_ReflectionType < 0 )
		{
			continue;
		}
		if ( texture_material[i].roughness )
		{
			continue;
		}
		float l1 = vectorlength_g(texture_material[i].diffuse.g);
		float l2 = vectorlength_g(texture_material[i].specular.g);
		if ( l1 > 0.001 && l2 > 0.001 )
		{
			texture_material[i].roughness = l1/(l1+l2);
		}
	}
	return true;
}

/*
** デストラクタ
*/
Obj::~Obj()
{
	vert.clear();
	norm.clear();
	norm_id.clear();
	fnorm.clear();
	face.clear();
	color.clear();
	uv.clear();
	texture_id.clear();
	printf("Delete Obj(%x)\n", (int)((void*)this));
}


void Obj::dump()
{
	FILE* fp = fopen("dump.obj", "w");

	for ( int i = 0; i < nv; i++ )
	{
		fprintf(fp, "v %f %f %f\n", vert[i].g[0], vert[i].g[1], vert[i].g[2]);
	}
	for ( int i = 0; i < nvn; i++ )
	{
		fprintf(fp, "vn %f %f %f\n", norm[i].g[0], norm[i].g[1], norm[i].g[2]);
	}
	for ( int i = 0; i < nf; i++ )
	{
		fprintf(fp, "f %d//%d %d//%d %d//%d\n", 
			face[i].g[0]+1, norm_id[i].g[0]+1, 
			face[i].g[1] + 1, norm_id[i].g[1] + 1,
			face[i].g[2] + 1, norm_id[i].g[2] + 1);
	}
	fclose(fp);
}

void Obj::dump2()
{
	FILE* fp = fopen("dump.obj", "w");

	for ( int i = 0; i < nv; i++ )
	{
		fprintf(fp, "v %f %f %f\n", vert[i].g[0], vert[i].g[1], vert[i].g[2]);
	}
	for ( int i = 0; i < nf; i++ )
	{
		fprintf(fp, "vn %f %f %f\n", fnorm[i].g[0], fnorm[i].g[1], fnorm[i].g[2]);
	}
	for ( int i = 0; i < nf; i++ )
	{
		fprintf(fp, "f %d//%d %d//%d %d//%d\n", 
			face[i].g[0]+1, face[i].g[0]+1, 
			face[i].g[1]+1, face[i].g[1]+1, 
			face[i].g[2]+1, face[i].g[2]+1);
	}
	fclose(fp);
}
