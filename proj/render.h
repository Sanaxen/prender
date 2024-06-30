#ifndef _RENDER_H_
#define _RENDER_H_

#include <iostream>
#include <algorithm>

#include "radiance.h"
#include "ppm.h"
#include "random.h"

#include "scene_env.h"
#include "def.h"

#include "bitmap.h"
#include "rgbe.h"
#include "hdr.h"
#include "mpiutil.h"

#include "mlt.h"
#include "erpt.h"

namespace prender {

// concentricにサンプリング
void concentric_sample_disk(double u1, double u2, double *dx, double *dy) 
{
    double r, theta;
    // [0, 1]の一様乱数u1,u2を[-1, 1]の一様乱数sx,syに写像
    const double sx = 2 * u1 - 1;
    const double sy = 2 * u2 - 1;


    // sx, syが0,0だった場合は特別に処理
    if (sx == 0.0 && sy == 0.0) {
        *dx = 0.0;
        *dy = 0.0;
        return;
    }
	// 四つに分割した円の各部位で別々の処理になる
    if (sx >= -sy) {
        if (sx > sy) {
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0f + sy/r;
        }
        else {
            r = sy;
            theta = 2.0f - sx/r;
        }
    }
    else {
        if (sx <= sy) {
            r = -sx;
            theta = 4.0f - sy/r;
        }
        else {
            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= PS_PI_4;
    *dx = r * cosf(theta);
    *dy = r * sinf(theta);
}


class CameraScreenPram
{
public:
	// カメラ位置
	Vector3d camera_position;
	Vector3d camera_dir;
	Vector3d camera_up;

	// ワールド座標系でのスクリーンの大きさ
	double screen_width;
	double screen_height;

	// スクリーンまでの距離
	double screen_dist;

	// スクリーンを張るベクトル
	Vector3d screen_x;
	Vector3d screen_y;
	Vector3d screen_center;


	int width;
	int height;
	int size;
	double iw;
	double ih;

	double focalDistance;
	double lensRadius;

	CameraScreenPram(SceneEnv* env_p)
	{
		// カメラ位置
		camera_position = env_p->camera_position;
		camera_dir = env_p->camera_dir;
		camera_up = env_p->camera_up;

		// ワールド座標系でのスクリーンの大きさ
		screen_width = env_p->world_screen_width * env_p->image_width / env_p->image_height;
		screen_height = env_p->world_screen_height;

		// スクリーンまでの距離
		screen_dist = env_p->world_screen_dist;

		// スクリーンを張るベクトル
		screen_x = normalize(cross(camera_dir, camera_up)) * screen_width;
		screen_y = normalize(cross(screen_x, camera_dir)) * screen_height;
		screen_center = camera_position + camera_dir * screen_dist;
		printf("%f,%f,%f\n", screen_x.x, screen_x.y, screen_x.z);
		printf("%f,%f,%f\n", screen_y.x, screen_y.y, screen_y.z);

		width = env_p->image_width;
		height = env_p->image_height;
		iw = 1.0 / width;
		ih = 1.0 / height;

		size = width*height;
		focalDistance = env_p->focalDistance;
		lensRadius = env_p->lensRadius;
	}
};

std::vector<MTRandom*> PrimarySampleRondomVector;
class Render
{
	SceneEnv* env_p;
	Color *image;
public:

#ifdef SPECTRUM_USE
	std::vector<int> samplineg_wavelength;
	int samplineg_wavelength_sz;
#endif
	double renderingTime;

	Render(SceneEnv& env)
	{
		env_p = &env;

#ifdef SPECTRUM_USE
		const int step = 5;
		for ( int i = 380; i<= 720; i += step)
		{
			samplineg_wavelength.push_back(i);
		}
		samplineg_wavelength_sz = samplineg_wavelength.size();
		printf("SPECTRUM %d-%d step %d sample %d\n", samplineg_wavelength[0], samplineg_wavelength[samplineg_wavelength.size() - 1], step, samplineg_wavelength.size());
#endif
		renderingTime = 0;
		logvalue = 0.0;
	}


private:
	double logvalue;
	// MLTのために新しいパスをサンプリングする関数。
	// 今回はradiance()（パストレ）を使ったが何でもいい。
	inline PathSample generate_new_path(const SceneEnv* env_p, const CameraScreenPram& cameraScreen, spectrum_df& spectrum_fnc, RandomMC *mlt, int xx, int yy, int local = 0)
	{
		double pdf = 1.0;
		float x = xx;
		float y = yy;

		if (0 && (xx != -1 && yy != -1))
		{
			const float x0 = x;
			const float y0 = y;
#if 0
			/*
			r1, r2 は摂動量を調整するパラメータ。
			MLT の元論文では、r1 は 0.1、r2 は画像面の面積(横ピクセル数 x 縦ピクセル数)となっている。
			A practical 〜 では r1 は 0.1、r2 は横ピクセル数の 10 % に設定している

			ちなみに、MLT のレンズ部分経路変異(lens subpath mutation) に相当する処理を行う場合は、
			A practical 〜 では pixel_offset() により大きな値を設定、r1 は 1.0、r2 は横ピクセル数の 25 %、することで代用している
			*/
			//A Practical Introduction to Metropolis Light	Transport
			//Appendix B Calculating a Pixel Offset
			//経路変異(lens subpath mutation) に相当する処理を行う場合は
			const float r1 = env_p->mlt_r1;
			const float r2 = cameraScreen.width*env_p->mlt_r2;
			//const float r2 = sqrt(cameraScreen.width*env_p->mlt_r2/PS_PI);

				if (logvalue == 0.0) logvalue = logf(r2 / r1);
			do{
				const float phi = mlt->NextSample()*PS_TWOPI;
				const float r = r2 * expf(-logvalue *mlt->NextSample());
				x = x0 + r*cosf(phi);
				y = y0 + r*sinf(phi);
			}while(x < 0 || cameraScreen.width <= x || y < 0 || cameraScreen.height <= y);
			pdf *= cameraScreen.size;

#else
			const int image_plane_mutation_value = cameraScreen.width*0.25;
			//double u, v;
			//
			//concentric_sample_disk(mlt->Next01(), mlt->Next01(), &u, &v);
			//
			//x = x + image_plane_mutation_value*u;
			//y = y + image_plane_mutation_value*v;
			//pdf *= image_plane_mutation_value*image_plane_mutation_value*PS_PI;
			do{
				const float r3 = image_plane_mutation_value;
				const float phi = mlt->NextSample()*PS_TWOPI;
				const float r = r3 * mlt->NextSample();
				x = round_(x0 + r*cosf(phi));
				y = round_(y0 + r*sinf(phi));
			}while(x < 0 || cameraScreen.width <= x || y < 0 || cameraScreen.height <= y);
			pdf *= cameraScreen.size;
#endif

			////0.5は小数点以下を四捨五入する。そうしないpixel_offset()が1未満だ摂動の度にとどんどんゼロになって居てしまう
			//x += 0.5f;
			//y += 0.5f;

			///* Immediately reject if we went off the image plane */
			//if (x < 0 || cameraScreen.width <= x || y < 0 || cameraScreen.height <= y)
			//{
			//	return PathSample(x, y, ZERO(), -1);
			//}
			xx = x;
			yy = y;
		}
		else
		{
			pdf *= cameraScreen.width;
			pdf *= cameraScreen.height;
			do{
				x = (mlt->NextSample() * cameraScreen.width);
				y = (mlt->NextSample() * cameraScreen.height);
			}while(x < 0 || cameraScreen.width <= x || y < 0 || cameraScreen.height <= y);
			xx = x;
			yy = y;
		}

		int sx = mlt->NextSample() < 0.5 ? 0 : 1;
		int sy = mlt->NextSample() < 0.5 ? 0 : 1;

		// テントフィルターによってサンプリング
		// ピクセル範囲で一様にサンプリングするのではなく、ピクセル中央付近にサンプルがたくさん集まるように偏りを生じさせる
		double dx = 0;
		double dy = 0;
		if (env_p->useTentFilter)
		{
			const double rr1 = 2.0 * mlt->NextSample();
			const double rr2 = 2.0 * mlt->NextSample();
			dx = rr1 < 1.0 ? sqrt(rr1) - 1.0 : 1.0 - sqrt(2.0 - rr1);
			dy = rr2 < 1.0 ? sqrt(rr2) - 1.0 : 1.0 - sqrt(2.0 - rr2);
			pdf *= 4.0;
		}
		const double r1 = (sx + 0.5 + dx)*0.5;
		const double r2 = (sy + 0.5 + dy)*0.5;

		////スクリーン座標に正規化
		//xx = (int)(r1 + x);
		//yy = (int)(r2 + y);
		//if (xx < 0 )
		//{
		//	xx = 0;
		//}
		//if (cameraScreen.width <= xx)
		//{
		//	xx = cameraScreen.width-1;
		//}
		//if ( yy < 0)
		//{
		//	yy = 0;
		//}
		//if ( cameraScreen.height <= yy)
		//{
		//	yy = cameraScreen.height-1;
		//}

		// スクリーン上の位置
		const Vector3d screen_position =
			cameraScreen.screen_center +
			cameraScreen.screen_x * ((r1 + xx) * cameraScreen.iw - 0.5) +
			cameraScreen.screen_y * ((r2 + yy) * cameraScreen.ih - 0.5);

		// レイを飛ばす方向
		const Vector3d dir = normalize(screen_position - cameraScreen.camera_position);


		Ray ray(cameraScreen.camera_position, dir);

#if 10
		// レイがレンズにあたって屈折することでDOFが発生する。
		if (cameraScreen.lensRadius > 0.0)
		{
			double lensU, lensV;
			concentric_sample_disk(mlt->NextSample(), mlt->NextSample(), &lensU, &lensV);
			lensU *= cameraScreen.lensRadius;
			lensV *= cameraScreen.lensRadius;
			double ft = fabs(cameraScreen.focalDistance / ray.dir.z);
			Vector3d Pfocus = ray.org + ray.dir * ft;

			ray.org = ray.org + Vector3d(lensU, lensV, 0.0);
			ray.dir = normalize(Pfocus - ray.org);
		}
#endif

		double wavelength_pdf = 1.0;
#ifdef SPECTRUM_USE
		double wavelength/* = rnd2.next01()*0.720+0.380*/;

		wavelength = (double)samplineg_wavelength[(int)(mlt->Next01()*(double)(samplineg_wavelength_sz - 1))];
		wavelength *= 0.001;

		//int wavelength_index;
		//wavelength = spectrum_fnc.wavelength_sampling(*mlt.rnd, wavelength_index, wavelength_pdf);
#else
		//屈折率を表す代表スペクトル
		const double wavelength = D_LINE_WAVELENGTH_SODIUM;
		double russian_roulette_probability = 1.0;
#endif
		int direct_light_hit = 0;
		Spectrum c = radiance(&direct_light_hit, env_p, ray, (Random*)mlt, 0, wavelength, env_p->nextEventEstimation, env_p->participatingMedia) / wavelength_pdf;
		return PathSample(xx, yy, Spectrum2RGB(wavelength)*c, pdf /*1.0/(1.0/pdf)*/, (direct_light_hit != 0));
	}


public:
	// MLTする
	// 画像平面上で大域的に突然変異(mutation)させる
	// Primary sample space MLT (PSSMLT) [Kelemen et al. 2002] 
	int RunMLT()
	{
		const CameraScreenPram cameraScreen(env_p);

		image = new Color[cameraScreen.size];

		spectrum_df spectrum_fnc;

		int num = 0;
		// OpenMP
		omp_set_num_threads(env_p->threads);

		clock_t startTime = clock();

		const clock_t startTimeD = startTime;

		const int nextEventEstimation = env_p->nextEventEstimation;
		const int participatingMedia = env_p->participatingMedia;

		int output_count = 0;

		int thred_max = 1;
#ifdef _OPENMP
		thred_max = omp_get_max_threads();
#endif
		std::vector<MTRandom*> rnd;
		for (int i = 0; i < thred_max; i++)
		{
			MTRandom* r = new MTRandom;
			r->seed(prime_numbers[i+32]);
			rnd.push_back(r);
		}
		PrimarySampleRondomVector = rnd;


		//1ピクセル当たりのサンプル数
		// このパスからMLTで使う最初のパスを得る。(Markov Chain Monte Carlo）
		int SeedPathMax = (int)(cameraScreen.size*env_p->pre_sample);


		//MLTの実施回数
		const int samlesNum = env_p->samples;
		printf("Markov Chain Monte Carlo sample(%d) %.2f/pixel\n", SeedPathMax, (double)samlesNum*SeedPathMax / (double)cameraScreen.size);
		fprintf(stderr, "Markov Chain Monte Carlo sample(%d) %.2f/pixel\n", SeedPathMax, (double)samlesNum*SeedPathMax / (double)cameraScreen.size);

		const bool luminance_sort = false;	//輝度値でソート

		bool timeOver = false;

		int mlt_last = samlesNum;
		int mlt_init = 0;
		int myrank = 0;
		char processor_name[128] = "0";
		char processor_name2[128] = "0";
#ifdef USE_MPI
		int procnum;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		MPI_Comm_size(MPI_COMM_WORLD, &procnum);
		int namelen;
		MPI_Get_processor_name(processor_name, &namelen);
		sprintf(processor_name2, "%s(%d)", processor_name, myrank);

		int dh = mlt_last / procnum;
		if (dh == 0) dh = 1;

		mlt_init = myrank*dh;
		mlt_last = mlt_init + dh;
		if (mlt_last + dh > samlesNum)
		{
			mlt_last = samlesNum;
		}
		output_count = mlt_init;
#endif
		const double mutation = cameraScreen.size*env_p->mutation;
		printf("mutation(%d) %.2f/pixel\n", (int)mutation, (double)samlesNum*mutation / (double)cameraScreen.size);
		fprintf(stderr, "mutation(%d) %.2f/pixel\n", (int)mutation, (double)samlesNum*mutation / (double)cameraScreen.size);

		
		// MLTする
#pragma omp parallel for schedule(dynamic, 1)
		for (int s = mlt_init; s < mlt_last; s++)
		{
			int progress = ((int)(10000.0 * ((s - mlt_init) + 1) / (mlt_last - mlt_init))) / 100.0;
			std::cerr << "[" << processor_name2 << "] " << "Rendering ( " << progress << " )% " << ((double)(clock() - startTimeD) / CLOCKS_PER_SEC) << "sec" << std::endl;

#ifdef _OPENMP
			if (mlt_init == 0 && mlt_last == samlesNum && myrank != 0)
			{
				continue;
			}
			const int thread_id = omp_get_thread_num();
#else
			const int thread_id = 0;
#endif
			if (env_p->timeLimit > 0)
			{
				clock_t endTime = clock();
				if (!timeOver && (double)(endTime - startTimeD) / (double)CLOCKS_PER_SEC > env_p->timeLimit*60.0*60.0)
				{
#pragma omp critical
					{
						timeOver = true;
						std::cerr << "Time Over!!" << std::endl;
						printf("Time Over!! [%d%%]\n", ((int)(10000.0 * s / (samlesNum))) / 100.0);
					}
					continue;
				}
			}

			std::vector<Color> tmp_image;
			tmp_image.resize(cameraScreen.size, Color());

			KelemenMLT mlt;
			mlt.large_step = 1;


			//沢山サンプリングして見る
			std::vector<PathSample> seed_paths(SeedPathMax);
			double sumI = 0.0;

			double luminance_max = 0.0;
			clock_t time_st = clock();
			//サンプル生成
			if (1)
			{
				for (int i = 0; i < SeedPathMax; i++)
				{
					mlt.InitUsedRandCoords();
					if ((i + 1) % (SeedPathMax / 10) == 0)
					{
						std::cerr << "\t[" << processor_name2 << "] (" << progress << "%) Sampling ( " << ((int)(10000.0 * (i + 1) / (SeedPathMax))) / 100.0 << " )%" << std::endl;
					}

					PathSample sample = generate_new_path(env_p, cameraScreen, spectrum_fnc, &mlt, -1, -1, 0);
					mlt.global_time++;
					mlt.StackClear();

					double luminance1 = luminance(sample.F);
					sumI += luminance1;
					if ( luminance1 > luminance_max) luminance_max = luminance1;
					seed_paths[i] = sample;
				}
#if 0
#if 10
				//昇順ソート
				if (luminance_sort) std::sort(seed_paths.begin(), seed_paths.end());
#else
				//降順ソート
				if (luminance_sort) std::sort(seed_paths.begin(), seed_paths.end(), std::greater<PathSample>());
#endif
#endif
			}
			if (sumI < PS_EPS14 )
			{
				continue;
			}
			std::cerr << "\t[" << processor_name2 << "] " << "Initial Sampling ( " << SeedPathMax << ") " << std::endl;


			//double total = 0.0;
			//std::vector<double> cdf_(SeedPathMax, 0.0);
			//std::vector<double> pdf_(SeedPathMax, 0.0);
			//for (int i = 0; i < SeedPathMax; i ++)
			//{
			//	total += luminance(seed_paths[i].F);
			//	cdf_[i] = total;
			//	pdf_[i] = luminance(seed_paths[i].F);
			//}
			//for (int i = 0; i < SeedPathMax; i ++)
			//{
			//	cdf_[i] /= total;
			//	pdf_[i] /= total;
			//}

			// 最初のパスを求める。輝度値に基づく重点サンプリング。
			int selecetd_path = SPPMLTSamplingPath(mlt, seed_paths, sumI, SeedPathMax);
			{
				int k = 0;
				while (k < SeedPathMax && seed_paths[selecetd_path].cancel)
				{
					selecetd_path = SPPMLTSamplingPath(mlt, seed_paths, sumI, SeedPathMax);
					k++;
				}
			}
			//見つけたパスを変異させる。

			// 論文参照
			const double b = sumI / SeedPathMax;
			const double p_large = 0.5;
			const int M = mutation;
			int accept = 0, reject = 0;

			PathSample old_path = seed_paths[selecetd_path];	//輝度値に基づく重点サンプリングで最初のパス
			//old_path.weight = 1.0/pdf_[selecetd_path];


			time_st = clock();
			//仕組みは[A Practical Introduction to Metropolis Light Transport]の2.2 Color Imagesを参照
			//Summarizing, the pseudo-code of the Metropolis light	transport algorithm is as follows を参照

			std::cerr << "\tmutation " << M << " start" << std::endl;

			const double iM = 1.0 / M;
			for (int i = 0; i < M; i++)
			{
				if ((i + 1) % (M / 10) == 0)
				{
					std::cerr << "\t\t[" << processor_name2 << "] (" << progress << "%) Mutation ( " << ((int)(10000.0 * (i + 1) / (M))) / 100.0 << " )% " << ((double)(clock() - time_st) / CLOCKS_PER_SEC) << "sec" << std::endl;
				}

				mlt.InitUsedRandCoords();
				mlt.large_step = (rnd[thread_id]->next01() < p_large) ? 1 : 0;


				//見つけたパスを変異させる。
				PathSample new_path;
				new_path = generate_new_path(env_p, cameraScreen, spectrum_fnc, &mlt, old_path.x, old_path.y, 1);

				// 以下全部論文と同じ（Next()）
				const double Lnew = luminance(new_path.F);
				const double Lold = luminance(old_path.F);
				double a = 1.0;
				
				if ( Lold > 1.0e-16 ) a = std::min(1.0, Lnew / Lold);				// accept prob.
				else a = 1.0;

				const double new_path_weight = (a + mlt.large_step) / (Lnew / b + p_large) * iM;	// accumulate weight
				const double old_path_weight = (1.0 - a) / (Lold / b + p_large) * iM;


				int new_image_index = (cameraScreen.height - new_path.y - 1) * cameraScreen.width + new_path.x;
				int old_image_index = (cameraScreen.height - old_path.y - 1) * cameraScreen.width + old_path.x;


				tmp_image[old_image_index] = tmp_image[old_image_index] + old_path_weight * old_path.weight * old_path.F;
				tmp_image[new_image_index] = tmp_image[new_image_index] + new_path_weight * new_path.weight * new_path.F;

				if (random()->next01() < a)
				{
					// 受理(accept)
					accept++;
					old_path = new_path;	//受理したので次はこのパスを変異させる
					if (mlt.large_step)
					{
						mlt.large_step_time = mlt.global_time;
					}
					mlt.global_time++;
					mlt.StackClear();
				}
				else
				{
					// 棄却
					reject++;
					mlt.RestoreState();
					//棄却したのでもう一回変異させてみる
				}



				//現時点の中間結果を出力する
#pragma omp critical
				{
					//現時点の中間結果を出力する
					clock_t endTime = clock();
					if ((double)(endTime - startTime) / (double)CLOCKS_PER_SEC > env_p->imageDumpTime)
					{
						const int sz = cameraScreen.size;
						const double w = env_p->sensor_response / (double)(mlt_last - mlt_init);

						//オリジナルはバックアップしておく
						Color* imgsv = new Color[sz];
						memcpy(imgsv, image, sizeof(Color)*sz);

						for (int i = 0; i < sz; i++)
						{
							image[i] = (image[i]  + tmp_image[i]) * w;
						}

						
						//輝度値に基づく重点サンプリングによって選んだ位置に印を付ける
						int x = seed_paths[selecetd_path].x;
						int y = seed_paths[selecetd_path].y;
						for (int i = -2; i <= 2; i++)
						{
							for (int j = -2; j <= 2; j++)
							{
								int index = (cameraScreen.height - (y + j) - 1) * cameraScreen.width + (x + i);
								if (index >= 0 && index < cameraScreen.size)
								{
									image[index] = Color(255 * 10, 0, 0);
								}
							}
						}
						x = old_path.x;
						y = old_path.y;
						for (int i = -2; i <= 2; i++)
						{
							for (int j = -2; j <= 2; j++)
							{
								int index = (cameraScreen.height - (y + j) - 1) * cameraScreen.width + (x + i);
								if (index >= 0 && index < cameraScreen.size)
								{
									image[index] = Color(0, 255 * 10, 0);
								}
							}
						}
						x = new_path.x;
						y = new_path.y;
						for (int i = -2; i <= 2; i++)
						{
							for (int j = -2; j <= 2; j++)
							{
								int index = (cameraScreen.height - (y + j) - 1) * cameraScreen.width + (x + i);
								if (index >= 0 && index < cameraScreen.size)
								{
									image[index] = Color(0, 0, 255 * 10);
								}
							}
						}

						OutputBmp(output_count);
						output_count++;
						startTime = endTime;
						std::cerr << "Image put!! [ " << output_count << " ]" << std::endl;

						memcpy(image, imgsv, sizeof(Color)*sz);
						delete[] imgsv;
					}
				}
			}
			std::cerr << "\t[" << processor_name2 << "] " << "Mutation " << ((double)(clock() - time_st) / CLOCKS_PER_SEC) << "sec " << std::endl;
			std::cerr << "\t[" << processor_name2 << "] Accept: " << accept << " Reject: " << reject << " Rate: " << (100.0 * accept / (accept + reject)) << "% " << std::endl;


			//現時点の中間結果を出力する
#pragma omp critical
			{
				const int sz = cameraScreen.size;
				for (int i = 0; i < sz; i++)
				{
					image[i] = image[i] + tmp_image[i];
				}

				//オリジナルはバックアップしておく
				Color* imgsv = new Color[sz];
				memcpy(imgsv, image, sizeof(Color)*sz);

				const double w = env_p->sensor_response / (double)(mlt_last - mlt_init);
				for (int i = 0; i < sz; i++)
				{
					image[i] = image[i] * w;
				}

				clock_t endTime = clock();
				if ((double)(endTime - startTime) / (double)CLOCKS_PER_SEC > env_p->imageDumpTime)
				{
					OutputBmp(output_count);
					output_count++;
					startTime = endTime;
					std::cerr << "Image put[" << processor_name2 << "]!! [ " << output_count << " ]" << std::endl;
				}
				memcpy(image, imgsv, sizeof(Color)*sz);
				delete [] imgsv;

			}
		}

		{
			const int sz = cameraScreen.size;
			const double w = env_p->sensor_response / (double)samlesNum;

			for (int i = 0; i < sz; i++)
			{
				image[i] = image[i] * w;
			}
		}
#ifdef USE_MPI
		if (env_p->timeLimit > 0 && timeOver)
		{
			/* empty */
		}
		else
		{
			fprintf(stderr, "MPI[%d]%s\n", myrank, processor_name);
			if (myrank == 0)
			{
				Color *image2 = new Color[cameraScreen.size];
				for (int i = 1; i < procnum; i++)
				{
					fprintf(stderr, "MPI_Recv[%d]%s\n", i, processor_name);
					MPI_Status status;
					memset(image2, '\0', sizeof(Color)*cameraScreen.size);
					MPI_Recv((void*)image2, sizeof(Color)*cameraScreen.size, MPI_BYTE, i, 99, MPI_COMM_WORLD, &status);

					const int sz = cameraScreen.size;
					for (int j = 0; j < sz; j++)
					{
						image[j] = image[j] + image2[j];
					}
				}
				delete[] image2;
			}
			else
			{
				fprintf(stderr, "MPI_Send[%d]%s\n", myrank, processor_name);
				MPI_Send((void*)image, sizeof(Color)*cameraScreen.size, MPI_BYTE, 0, 99, MPI_COMM_WORLD);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif

		for (int i = 0; i < PrimarySampleRondomVector.size(); i++)
		{
			delete  PrimarySampleRondomVector[i];
		}
		PrimarySampleRondomVector.clear();

		OutputBmp(output_count);
		OutputFinish(output_count);
		return 0;
	}


	//Energy redistribution PT (ERPT) [Cline et al. 2005] 
	int RunERPT()
	{
		const CameraScreenPram cameraScreen(env_p);

		const int samples = env_p->samples;
		const int supersamples = env_p->supersamples;

		image = new Color[cameraScreen.size];

		std::cout << cameraScreen.width << "x" << cameraScreen.height << " " << samples * (supersamples * supersamples) << " spp" << std::endl;

		spectrum_df spectrum_fnc;

		int num = 0;
		// OpenMP
		omp_set_num_threads(env_p->threads);

		clock_t startTime = clock();

		const clock_t startTimeD = startTime;

		const int nextEventEstimation = env_p->nextEventEstimation;
		const int participatingMedia = env_p->participatingMedia;

		const int mutation = (env_p->mutation < 1) ? 1 : (int)env_p->mutation;
		int output_count = 0;

		int thred_max = 1;
#ifdef _OPENMP
		thred_max = omp_get_max_threads();
#endif
		std::vector<MTRandom*> rnd;
		for (int i = 0; i < thred_max; i++)
		{
			MTRandom* r = new MTRandom;
			r->seed(prime_numbers[i + 32]);
			rnd.push_back(r);
		}
		PrimarySampleRondomVector = rnd;

		Color sumI = ZERO();

		fprintf(stderr, "Deposition Energy Estimation..");
		// edを求める
		const int ed_sample = 1;
#pragma omp parallel for schedule(dynamic, 1)/* reduction(+:num)*/
		for (int y = 0; y < cameraScreen.height; y++)
		{
#ifdef _OPENMP
			const int thread_id = omp_get_thread_num();
#else
			const int thread_id = 0;
#endif
			for (int x = 0; x < cameraScreen.width; x++)
			{
				ERPTSampler X;
				for (int ii = 0; ii < ed_sample; ii++)
				{
					sumI = sumI + generate_new_path(env_p, cameraScreen, spectrum_fnc, &X, x, y).F;
				}
			}
		}
		const double depositionEnergy = luminance(sumI / (cameraScreen.size*ed_sample)) / (double)mutation;
		fprintf(stderr, "depositionEnergy:%f\n", depositionEnergy); fflush(stderr);


		int height_last = cameraScreen.height;
		int height_init = 0;
		int myrank = 0;
		char processor_name[128] = "0";
		char processor_name2[128] = "0";
#ifdef USE_MPI
		int procnum;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		MPI_Comm_size(MPI_COMM_WORLD, &procnum);
		int namelen;
		MPI_Get_processor_name(processor_name, &namelen);
		sprintf(processor_name2, "%s(%d)", processor_name, myrank);

		int dh = height_last / procnum;

		height_init = myrank*dh;
		height_last = height_init + dh;
		if (height_last + dh > cameraScreen.height)
		{
			height_last = cameraScreen.height;
		}
		output_count = height_init;
#endif


		bool timeOver = false;
#pragma omp parallel for schedule(dynamic, 1)/* reduction(+:num)*/
		for (int y = height_init; y < height_last; y++)
		{

#ifdef _OPENMP
			if (height_init == 0 && height_last == cameraScreen.height && myrank != 0)
			{
				continue;
			}
#endif
			if (env_p->timeLimit > 0)
			{
				clock_t endTime = clock();
				if (!timeOver && (double)(endTime - startTimeD) / (double)CLOCKS_PER_SEC > env_p->timeLimit*60.0*60.0)
				{
#pragma omp critical
					{
						timeOver = true;
						std::cerr << "Time Over!!" << std::endl;
						printf("Time Over!! [%d%%]\n", ((int)(10000.0 * y / (cameraScreen.height))) / 100.0);
					}
					continue;
				}
			}

#ifdef _OPENMP
			const int thread_id = omp_get_thread_num();
#else
			const int thread_id = 0;
#endif
			std::vector<Color> tmp_image;
			tmp_image.resize(cameraScreen.size);

			std::cerr << "[" << processor_name2 << "] " << "Rendering ( " << ((int)(10000.0 * (y - height_init) / ((height_last - height_init) - 1))) / 100.0 << " )% " << ((double)(clock() - startTimeD) / CLOCKS_PER_SEC) << "sec" << std::endl;

			for (int x = 0; x < cameraScreen.width; x++)
			{
				for (int i = 0; i < supersamples*supersamples; i++)
				{
					ERPTSampler X;

					// 現在のスクリーン上のある点からのパスによる放射輝度を求める
					PathSample new_sample = generate_new_path(env_p, cameraScreen, spectrum_fnc, &X, x, y);
					const Color e = new_sample.F;

					const int new_sample_image_index = (cameraScreen.height - new_sample.y - 1) * cameraScreen.width + new_sample.x;

					// パスが光源に直接ヒットしてた場合、エネルギー分配しないで、そのまま画像に送る
					if (new_sample.direct_hit)
					{
						tmp_image[new_sample_image_index] = tmp_image[new_sample_image_index] + new_sample.F / samples;
						continue;
					}
					Color new_sample_image_color = tmp_image[new_sample_image_index] + new_sample.F / samples;

					// この辺は論文と同じ(Algorithm 2 Equal Deposition Flow)
					//どこかに遮蔽されたり背景にヒットするなどされずに非ゼロな輝度を持つのであれば
					if (luminance(e) > 0.0)
					{
						int numChains = std::floor(rnd[thread_id]->next01() + luminance(e) / (mutation * depositionEnergy));

						//if (numChains > 5)
						//{
						//	tmp_image[new_sample_image_index] = new_sample_image_color;
						//	continue;
						//}


						// Implementing Energy Redistribution Path Tracing参照
						// 周囲に分配するエネルギーがこれ
						const Color dep_value = e / luminance(e) * depositionEnergy / samples;

						//自分にも分配しておく
						//tmp_image[new_sample_image_index] = tmp_image[new_sample_image_index] + dep_value;

						//Algorithm 2 Equal Deposition Flow
						for (int nc = 0; nc < numChains; nc++)
						{
							ERPTSampler Y = X;
							PathSample Ypath = new_sample;

							// Consecutive sample filtering
							// ある点に極端にエネルギーが分配されると、スポットノイズになってしまう。
							// Unbiasedにするにはそれも仕方ないが、現実的には見苦しいのである点に対する分配回数を制限することでそのようなノイズを抑える
							// Biasedになるが、見た目は良くなる

							//Implementing Energy Redistribution Path Tracing
							//「連続サンプル」はスポットノイズになってしまう
							//反復場所の数をカウントすることによって分配回数を制限することでそのようなノイズを抑える
							//10または20のしきい値を超えた場合を制限し、それが離れて変異するまで、そのサンプルのエネルギーを堆積しない。
							const int MaxStack = 20;
							int stack_num = 0;
							int now_x = x, now_y = y;

							for (int im = 0; im < mutation; im++)
							{
								ERPTSampler Z = Y; Z.mutate();
								PathSample Zpath = generate_new_path(env_p, cameraScreen, spectrum_fnc, &Z, x, y, true);

								const double lfz = luminance(Zpath.F);
								const double lfy = luminance(Ypath.F);

								const double q = std::min(1.0, (lfy > PS_EPS16) ? lfz / lfy : 0.0);
								if (q > random()->next01())
								{
									Y = Z;
									Ypath = Zpath;
								}

								// Consecutive sample filtering
								int image_index = (cameraScreen.height - Ypath.y - 1) * cameraScreen.width + Ypath.x;
								if (now_x == Ypath.x && now_y == Ypath.y)
								{
									//同じピクセルに連続して蓄積。
									stack_num++;
								}
								else {
									//蓄積から離れて変異した。
									now_x = Ypath.x;
									now_y = Ypath.y;
									stack_num = 0;
								}

								// エネルギーをRedistributionする
								if (stack_num < MaxStack)
								{
									tmp_image[image_index] = tmp_image[image_index] + dep_value;
								}
							}
						}
					}
				}
				//現時点の中間結果を出力する
#pragma omp critical
				{
					//現時点の中間結果を出力する
					clock_t endTime = clock();
					if ((double)(endTime - startTime) / (double)CLOCKS_PER_SEC > env_p->imageDumpTime)
					{
						const int sz = cameraScreen.size;

						//オリジナルはバックアップしておく
						Color* imgsv = new Color[sz];
						memcpy(imgsv, image, sizeof(Color)*sz);

						//現在の中間結果生成
						for (int i = 0; i < sz; i++)
						{
							image[i] = image[i] + tmp_image[i] * env_p->sensor_response;;
						}

						OutputBmp(output_count);
						output_count++;
						startTime = endTime;
						std::cerr << "Image put!! [ " << output_count << " ]" << std::endl;

						//オリジナルに戻して置く
						memcpy(image, imgsv, sizeof(Color)*sz);
						delete[] imgsv;
					}
				}
			}
#pragma omp critical
			{
				const int sz = cameraScreen.size;
				for (int i = 0; i < sz; i++)
				{
					image[i] = image[i] + tmp_image[i];
				}

				clock_t endTime = clock();
				if ((double)(endTime - startTime) / (double)CLOCKS_PER_SEC > env_p->imageDumpTime)
				{
					const int sz = cameraScreen.size;

					//オリジナルはバックアップしておく
					Color* imgsv = new Color[sz];
					memcpy(imgsv, image, sizeof(Color)*sz);

					for (int i = 0; i < sz; i++)
					{
						image[i] = image[i] * env_p->sensor_response;
					}

					OutputBmp(output_count);
					output_count++;
					startTime = endTime;
					std::cerr << "Image put[" << processor_name2 << "]!! [ " << output_count << " ]" << std::endl;
					memcpy(image, imgsv, sizeof(Color)*sz);
					delete[] imgsv;

				}
			}
		}

#pragma omp parallel for schedule(dynamic, 1)
		for (int i = 0; i < 1; i++)
		{
			image[i] = image[i] * env_p->sensor_response;
		}

#ifdef USE_MPI
		if (env_p->timeLimit > 0 && timeOver)
		{
			/* empty */
		}
		else
		{
			fprintf(stderr, "MPI[%d]%s\n", myrank, processor_name);
			if (myrank == 0)
			{
				Color *image2 = new Color[cameraScreen.size];
				for (int i = 1; i < procnum; i++)
				{
					fprintf(stderr, "MPI_Recv[%d]%s\n", i, processor_name);
					MPI_Status status;
					memset(image2, '\0', sizeof(Color)*cameraScreen.size);
					MPI_Recv((void*)image2, sizeof(Color)*cameraScreen.size, MPI_BYTE, i, 99, MPI_COMM_WORLD, &status);

					const int sz = cameraScreen.size;
					for (int j = 0; j < sz; j++)
					{
						image[j] = image[j] + image2[j];
					}
				}
				delete[] image2;
			}
			else
			{
				fprintf(stderr, "MPI_Send[%d]%s\n", myrank, processor_name);
				MPI_Send((void*)image, sizeof(Color)*cameraScreen.size, MPI_BYTE, 0, 99, MPI_COMM_WORLD);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif
		for (int i = 0; i < PrimarySampleRondomVector.size(); i++)
		{
			delete PrimarySampleRondomVector[i];
		}
		PrimarySampleRondomVector.clear();

		OutputBmp(output_count);
		OutputFinish(output_count);
		return 0;
	}

	int Run()
	{
		const CameraScreenPram cameraScreen(env_p);

		int samples = env_p->samples;
		int supersamples = env_p->supersamples;
		double rate = (1.0 / supersamples);

		if (env_p->useTentFilter)
		{
			samples = samples*(supersamples * supersamples);
			supersamples = 1;
			rate = 0.5;
		}
		const double smp = 1.0 / samples / (supersamples * supersamples);


		image = new Color[cameraScreen.size];

		std::cout << cameraScreen.width << "x" << cameraScreen.height << " " << samples * (supersamples * supersamples) << " spp" << std::endl;



		spectrum_df spectrum_fnc;

		int num = 0;
		// OpenMP
		omp_set_num_threads(env_p->threads);

		clock_t startTime = clock();

		const clock_t startTimeD = startTime;

		const int nextEventEstimation = env_p->nextEventEstimation;
		const int participatingMedia = env_p->participatingMedia;

		int spp_count = 0;
		int output_count = 0;

		int thred_max = 1;
#ifdef _OPENMP
		thred_max = omp_get_max_threads();
#endif
		MTRandom *rnd = new MTRandom[thred_max];
		for (int i = 0; i < thred_max; i++)
		{
			rnd[i].seed(i+1);
		}


		int height_last = cameraScreen.height;
		int height_init = 0;
		int myrank = 0;
		char processor_name[128] = "0";
		char processor_name2[128] = "0";
#ifdef USE_MPI
		int procnum;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		MPI_Comm_size(MPI_COMM_WORLD, &procnum);
		int namelen;
		MPI_Get_processor_name(processor_name, &namelen);
		sprintf(processor_name2, "%s(%d)", processor_name, myrank);

		int dh = height_last / procnum;

		height_init = myrank*dh;
		height_last = height_init + dh;
		if (height_last + dh > cameraScreen.height)
		{
			height_last = cameraScreen.height;
		}
		output_count = height_init;
#endif

		// 一つのサブピクセルあたりsamples回サンプリングする
		for (int s = 0; s < samples; s++)
		{
#ifdef _OPENMP
			if (height_init == 0 && height_last == cameraScreen.height && myrank != 0)
			{
				continue;
			}
#endif
			std::cerr << "[" << processor_name2 << "] " << "Rendering ( " << ((int)(10000.0 * (s) / (samples))) / 100.0 << " )% " << ((double)(clock() - startTimeD) / CLOCKS_PER_SEC) << "sec" << std::endl;


#pragma omp parallel for schedule(dynamic, 1)/* reduction(+:num)*/
			for (int y = height_init; y < height_last; y++)
			{
#ifdef _OPENMP
				const int thread_id = omp_get_thread_num();
#else
				const int thread_id = 0;
#endif

				if ((height_last - height_init) > 10 && ((y - height_init) + 1) % ((height_last - height_init) / 10) == 0)
				{
					std::cerr << "\t[" << processor_name2 << "] " << "Rendering ( " << ((int)(10000.0 * (s) / (samples))) / 100.0 << " )% " << ((int)(10000.0 * (y - height_init) / ((height_last - height_init) - 1))) / 100.0 << " )%" << std::endl;
					//std::cerr << "[" << processor_name2 << "] " << "Rendering ( " << ((int)(10000.0 * (y - height_init) / ((height_last - height_init) - 1))) / 100.0 << " )%" << std::endl;
				}
				for (int x = 0; x < cameraScreen.width; x++)
				{
					//if ( cameraScreen.width > 10 && x % (cameraScreen.width/10) == 0 )
					//{
					//	std::cerr << "\t\t[" << processor_name2 << "] " << " ( " << ((int)(10000.0 * (s) / (samples))) / 100.0 << " )% " << ((int)(10000.0 * (x) / ((cameraScreen.width - 0) - 1))) / 100.0 << " )%" << std::endl;
					//}
					const int image_index = (cameraScreen.height - y - 1) * cameraScreen.width + x;

					// supersamples x supersamples のスーパーサンプリング
					for (int sy = 0; sy < supersamples; sy++)
					{
						for (int sx = 0; sx < supersamples; sx++)
						{
							Color accumulated_radiance = Color();;

							// 一つのサブピクセルあたりsamples回サンプリングする
							//for (int s = 0; s < samples; s++)
							{
								//const double r1 = sx * rate + rate * 0.5;
								//const double r2 = sy * rate + rate * 0.5;
								double dx = 0;
								double dy = 0;
								if (env_p->useTentFilter)
								{
									const double rr1 = 2.0 * rnd[thread_id].next01();
									const double rr2 = 2.0 * rnd[thread_id].next01();
									dx = rr1 < 1.0 ? sqrt(rr1) - 1.0 : 1.0 - sqrt(2.0 - rr1);
									dy = rr2 < 1.0 ? sqrt(rr2) - 1.0 : 1.0 - sqrt(2.0 - rr2);
								}
								const double r1 = (sx + 0.5 + dx)*rate;
								const double r2 = (sy + 0.5 + dy)*rate;

								// スクリーン上の位置
								const Vector3d screen_position =
									cameraScreen.screen_center +
									cameraScreen.screen_x * ((r1 + x) * cameraScreen.iw - 0.5) +
									cameraScreen.screen_y * ((r2 + y) * cameraScreen.ih - 0.5);

								// レイを飛ばす方向
								const Vector3d dir = normalize(screen_position - cameraScreen.camera_position);


								Ray ray(cameraScreen.camera_position, dir);

#if 10
								// レイがレンズにあたって屈折することでDOFが発生する。
								if (cameraScreen.lensRadius > 0.0)
								{
									double lensU, lensV;
									concentric_sample_disk(rnd[thread_id].next01(), rnd[thread_id].next01(), &lensU, &lensV);
									lensU *= cameraScreen.lensRadius;
									lensV *= cameraScreen.lensRadius;
									double ft = fabs(cameraScreen.focalDistance / ray.dir.z);
									Vector3d Pfocus = ray.org + ray.dir * ft;

									ray.org = ray.org + Vector3d(lensU, lensV, 0.0);
									ray.dir = normalize(Pfocus - ray.org);
								}
#endif
								double wavelength_pdf = 1.0;
#ifdef SPECTRUM_USE
								double wavelength/* = rnd2.next01()*0.720+0.380*/;

								wavelength = (double)samplineg_wavelength[(int)(rnd[thread_id].next01()*(double)(samplineg_wavelength_sz-1))];
								wavelength *= 0.001;

								//int wavelength_index;
								//wavelength = spectrum_fnc.wavelength_sampling(rnd[thread_id], wavelength_index, wavelength_pdf);
#else
								//屈折率を表す代表スペクトル
								const double wavelength = D_LINE_WAVELENGTH_SODIUM;
#endif
								//放射輝度の計算
								Spectrum v = radiance(0, env_p, ray, &(rnd[thread_id]), 0, wavelength, nextEventEstimation, participatingMedia)*smp;

								accumulated_radiance = accumulated_radiance + Spectrum2RGB(wavelength)*v / wavelength_pdf;

							}
							if (env_p->luminance_cutoff > 0.0)
							{
								float luminance_v = luminance(accumulated_radiance);
								if (luminance_v > env_p->luminance_cutoff)
								{
									float x = env_p->luminance_cutoff / luminance_v;
									accumulated_radiance = accumulated_radiance*x;
									//fprintf(stderr, "%f\n", accumulated_radiance.length());
								}
							}
							image[image_index] = image[image_index] + accumulated_radiance;
						}
					}
				}

#pragma omp critical
				{
					clock_t endTime = clock();
					if ((double)(endTime - startTime) / (double)CLOCKS_PER_SEC > env_p->imageDumpTime)
					{
						printf("[[%.2f%%]]\n", ((int)(10000.0 * (y - height_init) / ((height_last - height_init) - 1))) / 100.0);
						OutputBmp(output_count);
						output_count++;
						startTime = endTime;
						std::cerr << "Image put[" << processor_name2 << "]!! [ " << output_count << " ]" << std::endl;
					}
				}
			}
		}

#pragma omp parallel for schedule(dynamic, 1)
		for (int i = 0; i < cameraScreen.size; i++)
		{
			image[i] = image[i] * env_p->sensor_response;
		}

#ifdef USE_MPI
		fprintf(stderr, "MPI[%d]%s\n", myrank, processor_name);
		if (myrank == 0)
		{
			Color *image2 = new Color[cameraScreen.size];
			for (int i = 1; i < procnum; i++)
			{
				fprintf(stderr, "MPI_Recv[%d]%s\n", i, processor_name);
				MPI_Status status;
				memset(image2, '\0', sizeof(Color)*cameraScreen.size);
				MPI_Recv((void*)image2, sizeof(Color)*cameraScreen.size, MPI_BYTE, i, 99, MPI_COMM_WORLD, &status);

				int dh = cameraScreen.height / procnum;

				height_init = i*dh;
				height_last = height_init + dh;
				if (height_last + dh > cameraScreen.height)
				{
					height_last = cameraScreen.height;
				}

				for (int y = height_init; y < height_last; y++)
				{
					for (int x = 0; x < cameraScreen.width; x++)
					{
						const int image_index = (cameraScreen.height - y - 1) * cameraScreen.width + x;
						image[image_index] = image2[image_index];
					}
				}
			}
			delete[] image2;
		}
		else
		{
			fprintf(stderr, "MPI_Send[%d]%s\n", myrank, processor_name);
			MPI_Send((void*)image, sizeof(Color)*cameraScreen.size, MPI_BYTE, 0, 99, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		delete[] rnd;
		OutputBmp(output_count);
		OutputFinish(output_count);
		return 0;
	}


	void Output()
	{
		char drive[_MAX_DRIVE];	// ドライブ名
		char dir[_MAX_DIR];		// ディレクトリ名
		char fname[_MAX_FNAME];	// ファイル名
		char ext[_MAX_EXT];		// 拡張子

		_splitpath( env_p->output, drive, dir, fname, ext );

		std::cout << "Drive=" << drive << std::endl;
		std::cout << "Dir  =" << dir   << std::endl;
		std::cout << "Fname=" << fname << std::endl;
		std::cout << "Ext  =" << ext   << std::endl;

		char* suffix = getenv("SUFFIX_SYMBOL");
		if ( suffix == NULL ) suffix = "";

#ifdef	SPECTRUM_USE
		char* spectrum = "_spectrum";
#else
		char* spectrum = "";
#endif
		char filename[512];
		int spp = 1;
		if (env_p->energyRedistributionPathTracing)
		{
			sprintf(filename, "%s%s%s_%dx%d_erpt%s_%.1f%s%s%s", drive, dir, fname, env_p->samples, (int)env_p->mutation, env_p->nextEventEstimation ? "_nes" : "", renderingTime, spectrum, suffix, ext);
		}else
		if ( env_p->metropolisTransport )
		{
			spp = env_p->samples*(int)env_p->mutation;
			if ( env_p->mutation < 1 ) spp = env_p->samples;
			sprintf(filename, "%s%s%s_%dx%d_mlt%s_%.1f(%dspp)%s%s%s", drive, dir, fname, env_p->samples, (int)env_p->mutation, env_p->nextEventEstimation ? "_nes" : "", renderingTime, spp, spectrum, suffix, ext);
		}else
		{
			sprintf(filename, "%s%s%s_%dx%d%s_%.1f(%dspp)%s%s%s", drive, dir, fname, env_p->samples, env_p->supersamples, env_p->nextEventEstimation ? "_nes" : "", renderingTime, env_p->samples*env_p->supersamples*env_p->supersamples, spectrum, suffix, ext);
		}
		printf("output->[%s]\n", filename);

		if (strcmpi(ext, ".hdr") == 0)
		{
			OutputHdr(filename);
			return;
		}
		if (strcmpi(ext, ".bmp") == 0)
		{
			OutputBmp(filename);
			return;
		}

		// 出力
		const int width = env_p->image_width;
		const int height = env_p->image_height;
		save_ppm_file(std::string(filename), image, width, height);
	}

	void OutputBmp(int index = -1)
	{
		char drive[_MAX_DRIVE];	// ドライブ名
		char dir[_MAX_DIR];		// ディレクトリ名
		char fname[_MAX_FNAME];	// ファイル名
		char ext[_MAX_EXT];		// 拡張子

		_splitpath( env_p->output, drive, dir, fname, ext );

		std::cout << "Drive=" << drive << std::endl;
		std::cout << "Dir  =" << dir   << std::endl;
		std::cout << "Fname=" << fname << std::endl;
		std::cout << "Ext  =" << ext   << std::endl;

		char* suffix = getenv("SUFFIX_SYMBOL");
		if ( suffix == NULL ) suffix = "";

#ifdef	SPECTRUM_USE
		char* spectrum = "_spectrum";
#else
		char* spectrum = "";
#endif
		char filename[512];
		int spp = 1;
		if ( index < 0 )
		{
			if (env_p->energyRedistributionPathTracing)
			{
				sprintf(filename, "%s%s%s_%dx%d_erpt%s_%.1f%s%s%s", drive, dir, fname, env_p->samples, (int)env_p->mutation, env_p->nextEventEstimation ? "_nes" : "", renderingTime, spectrum, suffix, ".bmp");
			}
			else
			if (env_p->metropolisTransport)
			{
				spp = env_p->samples*(int)env_p->mutation;
				if ( env_p->mutation < 1 ) spp = env_p->samples;
				sprintf(filename, "%s%s%s_%dx%d_mlt%s_%.1f(%dspp)%s%s%s", drive, dir, fname, env_p->samples, (int)env_p->mutation, env_p->nextEventEstimation ? "_nes" : "", renderingTime, spp, spectrum, suffix, ".bmp");
			}
			else
			{
				sprintf(filename, "%s%s%s_%dx%d%s_%.1f(%dspp)%s%s%s", drive, dir, fname, env_p->samples, env_p->supersamples, env_p->nextEventEstimation ? "_nes" : "", renderingTime, env_p->samples*env_p->supersamples*env_p->supersamples, spectrum, suffix, ".bmp");
			}
			printf("output->[%s]\n", filename);
		}else
		{
			if ( drive[0] == '\0' ) strcpy(drive, ".");
			sprintf(filename, "%s%s\\image\\output_W%06d%s", drive, dir, index, ".bmp");
			printf("output->[%s]\n", filename);
		}
		// 出力
		const int width = env_p->image_width;
		const int height = env_p->image_height;

		BitMap bmp;
		bmp.Create(width, height);
		for ( int i = 0; i < height; i++ )
		{
			for ( int j = 0; j < width; j++ )
			{
				const int id = (height - i - 1) * width + j;
				bmp.cell(i,j) = Rgb(to_int(image[id].x), to_int(image[id].y), to_int(image[id].z));
			}
		}
		bmp.Write(filename);

		if(0)
		{
			float* src = ColotToFloat(image, width, height);
			float param_h = 0.01;
			float sigma = 0.01;
			
			filter_non_local_means(src, width, height, param_h, sigma);
			Color* new_image = FloatToColor(src, width, height);

			for (int i = 0; i < height; i++)
			{
				for (int j = 0; j < width; j++)
				{
					const int id = (height - i - 1) * width + j;
					bmp.cell(i, j) = Rgb(to_int(new_image[id].x), to_int(new_image[id].y), to_int(new_image[id].z));
				}
			}
			bmp.Write(filename);
			delete[] new_image;
		}
	}

	void OutputBmp(char* filename)
	{
		printf("output->[%s]\n", filename);

		// 出力
		const int width = env_p->image_width;
		const int height = env_p->image_height;

		BitMap bmp;
		bmp.Create(width, height);
		for ( int i = 0; i < height; i++ )
		{
			for ( int j = 0; j < width; j++ )
			{
				const int id = (height - i - 1) * width + j;
				bmp.cell(i,j) = Rgb(to_int(image[id].x), to_int(image[id].y), to_int(image[id].z));
			}
		}
		bmp.Write(filename);
	}

	void OutputHdr(char* filename)
	{
		printf("output->[%s]\n", filename);

		// 出力
		const int width = env_p->image_width;
		const int height = env_p->image_height;

		HDRImage hdr(width, height);
		for ( int i = 0; i < height; i++ )
		{
			for ( int j = 0; j < width; j++ )
			{
				const int id = (height - i - 1) * width + j;
				hdr.set(j, i, image[id]);
			}
		}
		hdr.save( std::string(filename));
	}

	void OutputFinish(int index)
	{
		char drive[_MAX_DRIVE];	// ドライブ名
		char dir[_MAX_DIR];		// ディレクトリ名
		char fname[_MAX_FNAME];	// ファイル名
		char ext[_MAX_EXT];		// 拡張子

		_splitpath(env_p->output, drive, dir, fname, ext);

		std::cout << "Drive=" << drive << std::endl;
		std::cout << "Dir  =" << dir << std::endl;
		std::cout << "Fname=" << fname << std::endl;
		std::cout << "Ext  =" << ext << std::endl;


		char filename[512];
		if (drive[0] == '\0') strcpy(drive, ".");
		sprintf(filename, "%s%s\\image\\Closed_%06d", drive, dir, index);
		printf("output->[%s]\n", filename);
		FILE* fp = fopen(filename, "w");
		if (fp)
		{
			fprintf(fp, "END\n");
			fclose(fp);
		}
	}
};
};

#endif
