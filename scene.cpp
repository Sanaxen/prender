#include "scene.h"

#include "KerrBlackHole.h"
#include "WormHole.h"

namespace prender {

	// レンダリングするシーンデータ
	ScenObjct* EntList;

	KerrBlackHole* BlackHole;
	WormHole* wormHole;

	std::vector<BitMap*> discTexture;
	std::vector<DiscMapping*> discMap;
	BitMap* backgroundTexture[2];
	SphericalMapping* backgroundMap[2];

	SceneEnv* Env;
	class prenderInit
	{
	public:
		prenderInit()
		{
			BlackHole = 0;
			backgroundTexture[0] = 0;
			backgroundTexture[1] = 0;
			backgroundMap[0] = 0;
			backgroundMap[1] = 0;
			wormHole = 0;
		}
	};

	prenderInit _prenderInit;

	void CreateEntList(SceneEnv& env)
	{

		Env = &env;
#if 10
		BlackHole = env.blackHole;
		wormHole = env.wormhole_;
		EntList = &env.EntList;
#else
		Sphere* spheres[10];
		spheres[0] = new Sphere(1e5, Vector3d( 1e5+1, 40.8, 81.6), Color(),      Color(0.75, 0.25, 0.25), REFLECTION_TYPE_DIFFUSE); // 左
		spheres[1] = new Sphere(1e5, Vector3d(-1e5+99, 40.8, 81.6),Color(),      Color(0.25, 0.25, 0.75), REFLECTION_TYPE_DIFFUSE); // 右
		spheres[2] = new Sphere(1e5, Vector3d(50, 40.8, 1e5),      Color(),      Color(0.75, 0.75, 0.75), REFLECTION_TYPE_DIFFUSE); // 奥
		spheres[3] = new Sphere(1e5, Vector3d(50, 40.8, -1e5+250), Color(),      Color(),                 REFLECTION_TYPE_DIFFUSE); // 手前
		spheres[4] = new Sphere(1e5, Vector3d(50, 1e5, 81.6),      Color(),      Color(0.75, 0.75, 0.75), REFLECTION_TYPE_DIFFUSE); // 床
		spheres[5] = new Sphere(1e5, Vector3d(50, -1e5+81.6, 81.6),Color(),      Color(0.75, 0.75, 0.75), REFLECTION_TYPE_DIFFUSE); // 天井
		spheres[6] = new Sphere(20,Vector3d(65, 20, 20),           Color(),      Color(0.25, 0.75, 0.25), REFLECTION_TYPE_DIFFUSE); // 緑球
		spheres[7] = new Sphere(16.5,Vector3d(27, 16.5, 47),       Color(),      Color(0.99, 0.99, 0.99), REFLECTION_TYPE_SPECULAR); // 鏡
		spheres[8] = new Sphere(15.0,Vector3d(50.0, 90.0, 81.6),   Color(36,36,36), Color(),  REFLECTION_TYPE_DIFFUSE); //照明
#if 0		
		spheres[9] = new Sphere(16.5,Vector3d(77, 16.5, 78),       Color(),      Color(0.99, 0.99, 0.99), REFLECTION_TYPE_REFRACTION); //ガラス

		for ( int i =0; i < 10;i++ ) EntList.Add( (Entity*)spheres[i] );
#else
		for ( int i =0; i < 9;i++ ) EntList.Add( (Entity*)spheres[i] );
		Obj *obj = new Obj("C:\\Users\\neutral\\Desktop\\レンダラ\\パストレーシング\\test\\sphere.obj");
		Polygon* poly = new Polygon(obj);
		poly->emission = Color();

#if 0
		poly->color =  Color(0.25, 0.75, 0.25);
		poly->reflection_type = REFLECTION_TYPE_DIFFUSE;
#else
		poly->color =  Color(0.99, 0.99, 0.99);
		poly->reflection_type = REFLECTION_TYPE_REFRACTION;
#endif
		poly->translate(77.0, 16.5, 78.0);
		EntList.Add( (Entity*)poly);
#endif

#endif
	}


	// シーンとの交差判定関数
	bool intersect_scene_(const Ray &ray, Intersection *intersection, const int depth,  int target_id)
	{
		// 初期化
		intersection->hitpoint.distance = PS_INF;
		intersection->object_id = -1;

		const size_t n = EntList->List.size();

		//自己ヒット防止
		const Ray ray2(ray.org + ray.dir*0.000001, ray.dir);

		for (int i = 0; i < n; i++)
		{
			if ( target_id >= 0 ) i = target_id;
			const Entity* entity = EntList->List[i];

			Intersection hit;
			if (entity->bvh)
			{
				if (entity->bvh->CheckIntersection(ray2, hit))
				{
					//printf("---------------START-------------\n");
					//printf("#BVH dist:%f %f %f %f  %f %f %f\n", 
					//	hit.hitpoint.distance,
					//		hit.hitpoint.position.x,
					//		hit.hitpoint.position.y,
					//		hit.hitpoint.position.z,
					//		hit.hitpoint.normal.x,
					//		hit.hitpoint.normal.y,
					//		hit.hitpoint.normal.z);

					Intersection hit3 = *intersection;

					if (hit.hitpoint.distance < intersection->hitpoint.distance)
					{
						*intersection = hit;
						intersection->object_id = i;
					}

					//総あたり計算（チェック用）
					//Intersection hit2;
					//if (entity->intersect(ray2, &hit2.hitpoint))
					//{
					//	if (hit2.hitpoint.distance < hit3.hitpoint.distance)
					//	{
					//		printf("*    dist %f %f %f %f  %f %f %f\n", 
					//			hit2.hitpoint.distance,
					//			hit2.hitpoint.position.x,
					//			hit2.hitpoint.position.y,
					//			hit2.hitpoint.position.z,
					//			hit2.hitpoint.normal.x,
					//			hit2.hitpoint.normal.y,
					//			hit2.hitpoint.normal.z);
					//		fflush(stdout);
					//		intersection->hitpoint = hit.hitpoint;
					//		intersection->object_id = i;
					//	}
					//}
					//printf("---------------END-------------\n\n");
				}
			}
			else
				if (entity->qbvh)
				{
					if (entity->qbvh->CheckIntersection(ray2, hit))
					{
						if (hit.hitpoint.distance < intersection->hitpoint.distance)
						{
							*intersection = hit;
							intersection->object_id = i;
						}
					}
				}
				else
				{
					//if (!entity->boundingBox.CheckIntersection(ray2))
					//{
					//	continue;
					//}
					Intersection hit;
					if (entity->intersect(ray2, &hit.hitpoint))
					{
						if (hit.hitpoint.distance < intersection->hitpoint.distance)
						{
							intersection->hitpoint = hit.hitpoint;
							intersection->object_id = i;
						}
					}
				}
				if ( target_id >= 0 ) break;
		}

		//背景に届いた
		if (intersection->object_id < 0 && BlackHole == 0 && wormHole == 0)
		{
			Spherical s = Cartesian(ray.dir).ToSpherical();
			if (backgroundMap[0] && backgroundTexture[0] && backgroundTexture[0]->GetImage())
			{
				s.th = s.th * Env->background_texture_map_coef[0][0] + Env->background_texture_map_coef[0][1];
				s.ph = s.ph * Env->background_texture_map_coef[0][2] + Env->background_texture_map_coef[0][3];
				s.th = fmod(s.th, M_PI);
				s.ph = fmod(s.ph, 2.0*M_PI);

				intersection->object_id = BACK_GLOUND;
				int xPos, yPos;

				backgroundMap[0]->Map(PS_INF, s.th, s.ph, xPos, yPos);
				Rgb t = backgroundTexture[0]->cell(yPos, xPos);
				intersection->hitpoint.material.color = Color(Env->background_texture_coef[0] * t.r / 255.0, Env->background_texture_coef[0] * t.g / 255.0, Env->background_texture_coef[0]*t.b / 255.0);
				//fprintf(stderr, "Hit background!!\n");
			}
		}

#if 0
		//真っ黒なオブジェクトの場合の特殊処理
		if (intersection->object_id >= 0)
		{
			if (intersection->hitpoint.material.color.sqr() < 1.0e-14 && intersection->hitpoint.material.emission.sqr() < 1.0e-14)
			{
				if (intersection->hitpoint.material.bump_texture || intersection->hitpoint.material.texture)
				{
					/* empty */
				}else
				{
					intersection->hitpoint.material.color = ONE();
					intersection->hitpoint.material.formulaType = REFLECTION_TYPE_REFRACTION;
				}
			}
		}
#endif
		return intersection->object_id >= 0;
	}



	inline float GetBrightness(Color& color)
	{
		//return luminance(color);
		float r = color.x;
		float g = color.y;
		float b = color.z;

		float max, min;

		max = r; min = r;

		if (g > max) max = g;
		if (b > max) max = b;

		if (g < min) min = g;
		if (b < min) min = b;

		return (max + min) / 2;
	}

	inline int Cap(float x, float max)
	{
		return x > max ? max : x;
	}

	inline int CapMin(float x, float min)
	{
		return x < min ? min : x;
	}

	inline Color AddColor(Color hitColor, Color tintColor)
	{

		float brightness = GetBrightness(tintColor);
		Color result = Color(
			Cap((int)((1 - brightness) * hitColor.x) + CapMin(tintColor.x - 20 / 255.0, 0) * (205 / 255.0), 1.0),
			Cap((int)((1 - brightness) * hitColor.y) + CapMin(tintColor.y - 20 / 255.0, 0) * (205 / 255.0), 1.0),
			Cap((int)((1 - brightness) * hitColor.z) + CapMin(tintColor.z - 20 / 255.0, 0) * (205 / 255.0), 1.0)
			);
		return result;
	}


	std::vector<Vector3d> pointList;
	void PointsDump()
	{
		if (pointList.size() < 1000) return;

		FILE* fp = fopen("points.obj", "w");

		for (int i = 0; i < pointList.size(); i++)
		{
			Vector3d p = pointList[i];
			fprintf(fp, "v %f %f %f\n", p.x, p.y, p.z);
		}
		fclose(fp);
		exit(0);
	}


	void BlackHoleTextureInit(char* drive, char* dir, KerrBlackHole* kb)
	{
		fprintf(stderr, "=============\n");

		discTexture.resize(kb->accretion_disk_texture.size());
		discMap.resize(kb->accretion_disk_texture.size());

		for (int i = 0; i < kb->accretion_disk_texture.size(); i++)
		{
			if (!discMap[i])
			{
				bool texture_exist = false;

				if (!discTexture[i])
				{
					char filename2[1024];
					FILE* fp = fopen(kb->accretion_disk_texture[i].c_str(), "r");
					if (fp)
					{
						printf("texture=[%s]\n", kb->accretion_disk_texture[i].c_str());
						fclose(fp);
						texture_exist = true;
					}
					else
					{
						sprintf(filename2, "%s%s%s", drive, dir, kb->accretion_disk_texture[i].c_str());
						kb->accretion_disk_texture[i] = filename2;
						printf("texture=[%s]\n", filename2);
						fp = fopen(filename2, "r");
						if (fp)
						{
							fclose(fp);
							texture_exist = true;
						}
						else
						{
							printf("accretion_disk_textureファイルが開けません[%s]\n", kb->accretion_disk_texture[i].c_str());
						}
					}

					if (texture_exist)
					{
						discTexture[i] = new BitMap;
						discTexture[i]->Read((char*)kb->accretion_disk_texture[i].c_str());
						if (!discTexture[i]->GetImage())
						{
							fprintf(stderr, "discTexture error.\n");
						}
						else
						{
							printf("discTexture %d x %d\n", discTexture[i]->W(), discTexture[i]->H());
							fprintf(stderr, "discTexture %d x %d\n", discTexture[i]->W(), discTexture[i]->H());
							fprintf(stderr, "Rmstable %f Rdisk %f\n", kb->Rmstable, kb->Rdisk);
							discMap[i] = new DiscMapping(kb->Rmstable, kb->Rdisk, discTexture[i]->W(), discTexture[i]->H());
						}
					}
				}
			}
		}

		if (!backgroundMap[0])
		{
			bool texture_exist = false;
			if (!backgroundTexture[0])
			{
				char filename2[1024];
				FILE* fp = fopen(kb->background_texture.c_str(), "r");
				if (fp)
				{
					printf("texture=[%s]\n", kb->background_texture.c_str());
					fclose(fp);
					texture_exist = true;
				}
				else
				{
					sprintf(filename2, "%s%s%s", drive, dir, kb->background_texture.c_str());
					kb->background_texture = filename2;
					printf("texture=[%s]\n", filename2);
					fp = fopen(filename2, "r");
					if (fp)
					{
						fclose(fp);
						texture_exist = true;
					}
					else
					{
						printf("background_textureファイルが開けません[%s]\n", kb->background_texture.c_str());
					}
				}

				if (texture_exist)
				{
					backgroundTexture[0] = new BitMap;
					backgroundTexture[0]->Read((char*)kb->background_texture.c_str());
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
	}

	void WormHoleTextureInit(char* drive, char* dir, WormHole* kb)
	{
		fprintf(stderr, "=============\n");


		for (int ii = 0; ii < 2; ii++)
		{
			if (!backgroundMap[ii])
			{
				bool texture_exist = false;
				if (!backgroundTexture[ii])
				{
					char filename2[1024];
					FILE* fp = fopen(kb->background_texture[ii].c_str(), "r");
					if (fp)
					{
						printf("texture=[%s]\n", kb->background_texture[ii].c_str());
						fclose(fp);
						texture_exist = true;
					}
					else
					{
						sprintf(filename2, "%s%s%s", drive, dir, kb->background_texture[ii].c_str());
						kb->background_texture[ii] = filename2;
						printf("texture[%d]=[%s]\n", ii, filename2);
						fp = fopen(filename2, "r");
						if (fp)
						{
							fclose(fp);
							texture_exist = true;
						}
						else
						{
							printf("background_textureファイルが開けません[%s]\n", kb->background_texture[ii].c_str());
						}
					}

					if (texture_exist)
					{
						backgroundTexture[ii] = new BitMap;
						backgroundTexture[ii]->Read((char*)kb->background_texture[ii].c_str());
						if (!backgroundTexture[ii]->GetImage())
						{
							fprintf(stderr, "backgroundTexture error.\n");
						}
						else
						{
							printf("backgroundTexture[%d] %d x %d\n", ii, backgroundTexture[ii]->W(), backgroundTexture[ii]->H());
							fprintf(stderr, "backgroundTexture[%d] %d x %d\n", ii, backgroundTexture[ii]->W(), backgroundTexture[ii]->H());
							backgroundMap[ii] = new SphericalMapping(backgroundTexture[ii]->W(), backgroundTexture[ii]->H());
						}
					}
				}
			}
		}
	}

	int Ok = 0;
	int Ng = 0;

	// シーンとの交差判定関数
	bool intersect_scene(const Ray &ray, Intersection *intersection, const int depth, int target_id)
	{
		if (BlackHole == 0 && wormHole == 0) return intersect_scene_(ray, intersection, depth);
		if (wormHole) return intersect_scene__(ray, intersection, depth);

		//if (depth > 0) return intersect_scene_(ray, intersection, depth);

		KerrBlackHole kb = *BlackHole;


		//ブラックホールの中心を(0,0,0)にする座標系に変更
		if (!kb.Setup(ray.org.x, ray.org.y, ray.org.z))
		{
			fprintf(stderr, "座標系エラー\n");
			return false;
		}


		double y[KERRBLACKHOLE_ODE_N];
		double dydx[KERRBLACKHOLE_ODE_N];
		double y_next[KERRBLACKHOLE_ODE_N];
		double dydx_next[KERRBLACKHOLE_ODE_N];
		double ylaststep[KERRBLACKHOLE_ODE_N];

		KerrBlackHoleWrkParam prm;

		//レイ・ベクトルをブラックホールの中心を(0,0,0)にする座標系に変換
		Vector3d dir = ray.dir;

		//事象の地平面に到達
		if (kb.r0 < kb.Rhor)
		{
			//fprintf(stderr, "############事象の地平面に到達\n");
			intersection->object_id = EVENT_HORIZON;
			return false;
		}


		y[0] = kb.r0;
		y[1] = kb.theta0;
		y[2] = kb.phi0;
		y[3] = 0;
		y[4] = 0;
		y[5] = 0;


		kb.traceDir = 1.0;
		//RayからKerrBlackHoleクラスの内部座標でのrayにしたときの微分方程式初期条件を設定する
		int stat = kb.RayToOrdinaryDifferentialEquationInitial(depth, prm, 1.0*dir, y, dydx);
		if (stat < 0)
		{
			y[0] = kb.r0;
			y[1] = kb.theta0;
			y[2] = kb.phi0;
			 kb.traceDir = -1.0;
			 stat = kb.RayToOrdinaryDifferentialEquationInitial(depth, prm, -1.0*dir, y, dydx);
		}
		if (stat == 0)
		{
			Ok++;
			if ((Ok + Ng) % 10000 == 0) fprintf(stderr, "\t\tTracking of geodesic success rate %d / %d  %.3f%%\n", Ok, (Ok + Ng), 100.0*(double)Ok / (double)(Ok + Ng));
		}
		else
		{
			Ng++;
		}
		if (stat < 0)
		{
			fprintf(stderr, "\t\tsuccess rate %d / %d  %.3f%%\n", Ok, (Ok + Ng), 100.0*(double)Ok / (double)(Ok + Ng));
			return false;
			bool s = intersect_scene_(ray, intersection, depth);
			//fprintf(stderr, "---------------------------\n");
			if (s)
			{
				//if (fabs(intersection->hitpoint.material.color.x - 1.0) < 1.0e-10 && fabs(intersection->hitpoint.material.color.y) < 1.0e-10 && fabs(intersection->hitpoint.material.color.z) < 1.0e-10)
				//{
				//	//fprintf(stderr, "############事象の地平面に到達 %f %f\n", kb.Rhor, y[0]);
				//	return false;
				//}
				return s;
			}
			return s;
		}

		//fprintf(stderr, "---------------------------\n");
		int side;

		int cnt = 0;
		int cntMax = kb.geodesics_max_length;
		Vector3d tnv;
		double hdid;

		Intersection pre_Intersection;
		pre_Intersection.object_id = -1;
		int pre_from_outside = -1;

		Vector3d org2 = ray.org;
		Vector3d dir2 = ray.dir;

		Vector3d org1 = ray.org;

		//Color pixel(-1, 0, 0);
		//Color hitPixel(0, 0, 0);

		double h_dist = 2.0;
		double h_step = h_dist;
		while ( (cntMax < 0) ? 1 : (cnt < cntMax) )
		{
			double delta = (org2 -org1).length();

			//世界座標系の現在のレイ（測地線の接ベクトル）で交点計算
			if (intersect_scene_(Ray(org1, dir2), intersection, depth))
			{
				double d = h_dist;

				//交点があって対象にぶつかる直前なら交点とする
				if (intersection->hitpoint.distance < d )
				{
					Cartesian cc(kb.Coordinate_transformation(intersection->hitpoint.position));
					Spherical ss = cc.ToBoyerLindquist(kb.a);

					if (ss.r < kb.Rhor)
					{
						//fprintf(stderr, "############事象の地平面に到達 %d\n", __LINE__);
						intersection->object_id = EVENT_HORIZON;
						return false;
					}
					return true;
				}
			}/*else
			{
				fprintf(stderr, "no intersection\n");
				return false;
			}*/

			kb.geodesic(prm, y, dydx);

#if 10
			if ( fabs(h_step) < 1.0e-10 )
			{
				//fprintf(stderr, "%.16f\n", h_step);
				return false;
				h_step = 1.0e-10;
			}
			h_step = kb.NextPosition(prm, y, dydx, y_next, dydx_next, h_step, 0/*&tnv*/, &hdid);
			if (h_step > h_dist) h_step = h_dist;
#else
			Spherical p1(y[0], y[1], y[2]);
			const int chkloopMax = 300;
			for (int kk = 0; kk < chkloopMax; kk++)
			{
				hdid = 0;
				//次の位置を求める
				if (h_step < 1.0e-6) h_step = 1.0e-6;
				h_step = kb.NextPosition(prm, y, dydx, y_next, dydx_next, h_step, 0/*&tnv*/, &hdid);
				if (h_step > h_dist) h_step = h_dist;

				//事象の地平面に到達
				if (y_next[0] < kb.Rhor)
				{
					//fprintf(stderr, "############事象の地平面に到達 %d\n", __LINE__);
					intersection->object_id = EVENT_HORIZON;
					return false;
				}


				Spherical p2(y_next[0], y_next[1], y_next[2]);

				//２点間距離を求める
				double ds = (kb.WorldAxisSystem(p2.ToBoyerLindquist(kb.a)) - kb.WorldAxisSystem(p1.ToBoyerLindquist(kb.a))).length();

				if (ds < 1.0e-6)
				{
					//fprintf(stderr, "--%d----- %.16f-----------------------------\n", kk, ds);
				}
				else
				{
					//fprintf(stderr, "OK!!--%d----- %.16f-----------------------------\n\n", kk, ds);
					break;
				}

				if (kk == chkloopMax-1)
				{
					fprintf(stderr, "ERROR!!--%d----- %.16f--(%f)---------------------------\n\n", kk, ds, h_step);
					intersection->object_id = -1;
					return false;
				}
				memcpy(y, y_next, KERRBLACKHOLE_ODE_N*sizeof(double));
				memcpy(dydx, dydx_next, KERRBLACKHOLE_ODE_N*sizeof(double));

				kb.geodesic(prm, y, dydx);
			}
#endif
			memcpy(y, y_next, KERRBLACKHOLE_ODE_N*sizeof(double));
			memcpy(dydx, dydx_next, KERRBLACKHOLE_ODE_N*sizeof(double));


			//事象の地平面に到達
			if (y[0] < kb.Rhor)
			{
				//fprintf(stderr, "############事象の地平面に到達 %d\n", __LINE__);
				intersection->object_id = EVENT_HORIZON;
				return false;
			}


			org1 = org2;	//現在の位置を覚えておく

			//次の位置に更新する
			//微小移動した位置を世界座標系に戻す
			Spherical sp(y[0], y[1], y[2]);
			// pointList.push_back(sp.ToVector());

			//org2 = kb.WorldAxisSystem(sp.ToCartesian());
			org2 = kb.WorldAxisSystem(sp.ToBoyerLindquist(kb.a));
			//fprintf(stderr, "%05d %f %f %f\n", cnt, org2.x, org2.y, org2.z);

			//dir2 = tnv;	//本当の接ベクトル
			dir2 = normalize(org2 - org1);	//ポリラインとして計算してるので折れ線方向を接ベクトルとする。
			cnt++;
		}
		//fprintf(stderr, "%d\n", cnt);
		// PointsDump();


		//fprintf(stderr, "Tracking of geodesic Limit!![%d]\n", cnt);
		//Spherical sp(y[0], y[1], y[2]);
		//org2 = kb.WorldAxisSystem(sp);

		if (intersect_scene_(Ray(org2, dir2), intersection, depth))
		{
			//交点があるから背景には届いていない
			return false;

			Cartesian cc(intersection->hitpoint.position);
			Spherical ss = cc.ToBoyerLindquist(kb.a);
			if (ss.r < kb.Rhor)
			{
				//fprintf(stderr, "############事象の地平面に到達 %d\n", __LINE__);
				intersection->object_id = EVENT_HORIZON;
				return false;
			}

			//hitPixel = intersection->hitpoint.material.color;
			////if (pixel.x >= 0)
			//{
			// intersection->hitpoint.material.color = AddColor(hitPixel, pixel);
			//}
			return true;
		}
		else
		{
#if 10
			//背景に届いた
			//fprintf(stderr, "%f %f %f\n", org2.x, org2.y, org2.z);
			if (backgroundMap[0] && backgroundTexture[0] && backgroundTexture[0]->GetImage())
			{
				y[1] = y[1]*kb.background_texture_map_coef[0] + kb.background_texture_map_coef[1];
				y[2] = y[2]*kb.background_texture_map_coef[2] + kb.background_texture_map_coef[3];
				y[1] = fmod(y[1], M_PI);
				y[2] = fmod(y[2], 2.0*M_PI);

				intersection->object_id = BACK_GLOUND;
				int xPos, yPos;

				backgroundMap[0]->Map(y[0], y[1], y[2], xPos, yPos);
				Rgb t = backgroundTexture[0]->cell(yPos, xPos);
				intersection->hitpoint.material.color = Color(kb.background_texture_coef*t.r / 255.0, kb.background_texture_coef*t.g / 255.0, kb.background_texture_coef*t.b / 255.0);
				//fprintf(stderr, "Hit background!!\n");
			}
#endif
		}
		return false;
	}


	// シーンとの交差判定関数
	bool intersect_scene__(const Ray &ray, Intersection *intersection, const int depth,  int target_id)
	{
		WormHole kb = *wormHole;

		Spherical rayS = Cartesian(ray.dir.x, ray.dir.y, ray.dir.z).ToSpherical();

		double sign = 1.0;
		if ( ray.org.length() > 1000000)
		{
			sign = -1.0;
		}
		kb.sign = sign;

		//ブラックホールの中心を(0,0,0)にする座標系に変更
		if (!kb.Setup(ray.org.x, ray.org.y, ray.org.z, rayS, sign))
		{
			fprintf(stderr, "座標系エラー\n");
			return false;
		}


		double y[WORMHOLE_ODE_N];
		double dydx[WORMHOLE_ODE_N];
		double y_next[WORMHOLE_ODE_N];
		double dydx_next[WORMHOLE_ODE_N];

		//レイ・ベクトルをワームホールの中心を(0,0,0)にする座標系に変換
		Vector3d dir = ray.dir;
		Vector3d tnv = dir;

		kb.initial(y, dydx, rayS, tnv);


		int cnt = 0;
		int cntMax = kb.geodesics_max_length;
		double hdid;

		Intersection pre_Intersection;
		pre_Intersection.object_id = -1;
		int pre_from_outside = -1;

		Vector3d org2 = ray.org;
		Vector3d dir2 = ray.dir;

		Vector3d org1 = ray.org;


		double h_dist = 0.5;
		double h_step = h_dist;
		while ((cntMax < 0) ? 1 : (cnt < cntMax))
		{
			if (fabs(y[0]) > 1.0e10 || fabs(kb.r_wh(y[0])) > 1.0e10)
			{
				//fprintf(stderr, " l=%f cnt %d\n",  y[0], cnt);
				break;
			}

			//double delta = (org2 - org1).length();

			//世界座標系の現在のレイ（測地線の接ベクトル）で交点計算
			if (intersect_scene_(Ray(org1, dir2), intersection, depth))
			{
				//fprintf(stderr, "hit object!!\n");
				double d = h_dist;

				//交点があって対象にぶつかる直前なら交点とする
				if (intersection->hitpoint.distance < d)
				{
					return true;
				}
			}/*else
			 {
			 fprintf(stderr, "no intersection\n");
			 return false;
			 }*/

			kb.geodesic( y, dydx);

			if (fabs(h_step) < 1.0e-14)
			{
				h_step = 1.0e-14;
			}
			h_step = kb.NextPosition( y, dydx, y_next, dydx_next, h_step, 0/*&tnv*/, &hdid);
			//fprintf(stderr, "%d %f->%f\n", cnt, y[0], y_next[0]);

			if (h_step > h_dist) h_step = h_dist;
			memcpy(y, y_next, WORMHOLE_ODE_N*sizeof(double));
			memcpy(dydx, dydx_next, WORMHOLE_ODE_N*sizeof(double));


			org1 = org2;	//現在の位置を覚えておく

			//次の位置に更新する
			//微小移動した位置を世界座標系に戻す
			Spherical sp(kb.r_wh(y[0]), y[1], y[2]);

			org2 = kb.WorldAxisSystem(sp.ToCartesian(), y[0]);
			//fprintf(stderr, "%05d %f %f %f\n", cnt, org2.x, org2.y, org2.z);

			//dir2 = tnv;	//本当の接ベクトル
			dir2 = normalize(org2 - org1);	//ポリラインとして計算してるので折れ線方向を接ベクトルとする。
			cnt++;

		}
		//fprintf(stderr, "%d\n", cnt);
		// PointsDump();


		//fprintf(stderr, "Tracking of geodesic Limit!![%d]\n", cnt);
		//Spherical sp(y[0], y[1], y[2]);
		//org2 = kb.WorldAxisSystem(sp);

		//if (intersect_scene_(Ray(org2, dir2), intersection, depth))
		//{
		//	//交点があるから背景には届いていない
		//	return false;
		//}
		//else
		{
#if 10
			Spherical s = Cartesian(dir2).ToSpherical();
			y[1] = s.th;
			y[2] = s.ph;

			//背景に届いた
			if (y[0] >= 0)
			{
				y[1] = y[1]*kb.background_texture_map_coef[0][0] + kb.background_texture_map_coef[0][1];
				y[2] = y[2]*kb.background_texture_map_coef[0][2] + kb.background_texture_map_coef[0][3];
				y[1] = fmod(y[1], M_PI);
				y[2] = fmod(y[2], 2.0*M_PI);

				//fprintf(stderr, "背景に届いた l=%f\n", y[0]);
				//fprintf(stderr, "%f %f %f\n", org2.x, org2.y, org2.z);
				if (backgroundMap[0] && backgroundTexture[0] && backgroundTexture[0]->GetImage())
				{
					intersection->object_id = BACK_GLOUND;
					int xPos, yPos;

					backgroundMap[0]->Map(kb.r_wh(y[0]), y[1], y[2], xPos, yPos);
					Rgb t = backgroundTexture[0]->cell(yPos, xPos);
					intersection->hitpoint.material.color = Color(kb.background_texture_coef[0]*t.r / 255.0, kb.background_texture_coef[0]*t.g / 255.0, kb.background_texture_coef[0]*t.b / 255.0);
					//fprintf(stderr, "Hit background!!\n");
				}
			}
			if (y[0] < 0)
			{
				y[1] = y[1]*kb.background_texture_map_coef[1][0] + kb.background_texture_map_coef[1][1];
				y[2] = y[2]*kb.background_texture_map_coef[1][2] + kb.background_texture_map_coef[1][3];
				y[1] = fmod(y[1], M_PI);
				y[2] = fmod(y[2], 2.0*M_PI);

				//fprintf(stderr, "通過完了 l=%f\n", y[0]);
				//背景に届いた
				//fprintf(stderr, "%f %f %f\n", org2.x, org2.y, org2.z);
				if (backgroundMap[1] && backgroundTexture[1] && backgroundTexture[1]->GetImage())
				{
					intersection->object_id = BACK_GLOUND2;
					int xPos, yPos;

					backgroundMap[1]->Map(kb.r_wh(y[0]), y[1], y[2], xPos, yPos);
					Rgb t = backgroundTexture[1]->cell(yPos, xPos);
					intersection->hitpoint.material.color = Color(kb.background_texture_coef[1]*t.r / 255.0, kb.background_texture_coef[1]*t.g / 255.0, kb.background_texture_coef[1]*t.b / 255.0);
					//fprintf(stderr, "Hit background!!\n");
				}
			}
#endif
		}
		return false;
	}

};
