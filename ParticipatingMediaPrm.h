#ifndef _PARTICIPATING_MEDIA_PRM_H

#define _PARTICIPATING_MEDIA_PRM_H


#include "vector3d.h"
#include "def.h"
#include "Spectrum.h"

namespace prender {

#ifdef SPECTRUM_USE
	double RGB2Spectrum(const Color rgb, const double lambda);
#endif
	
static struct material_symbol_t
{
	char* name;
	float scattering[3];
	float absorbing[3];
	float g[3];
	float eta;
}material_symbol_[]=
{
	/* Fitted data from "A Practical Model for Subsurface scattering" (Jensen et al.). No anisotropy data available. */
	{ "Apple", { 2.29f, 2.39f, 1.97f }, { 0.0030f, 0.0034f, 0.046f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Chicken1", { 0.15f, 0.21f, 0.38f }, { 0.0015f, 0.077f, 0.19f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Chicken2", { 0.19f, 0.25f, 0.32f }, { 0.0018f, 0.088f, 0.20f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Cream", { 7.38f, 5.47f, 3.15f }, { 0.0002f, 0.0028f, 0.0163f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Ketchup", { 0.18f, 0.07f, 0.03f }, { 0.061f, 0.97f, 1.45f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Marble", { 2.19f, 2.62f, 3.00f }, { 0.0021f, 0.0041f, 0.0071f }, { 0.0f, 0.0f, 0.0f }, 1.5f },
	{ "Potato", { 0.68f, 0.70f, 0.55f }, { 0.0024f, 0.0090f, 0.12f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Skimmilk", { 0.70f, 1.22f, 1.90f }, { 0.0014f, 0.0025f, 0.0142f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Skin1", { 0.74f, 0.88f, 1.01f }, { 0.032f, 0.17f, 0.48f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Skin2", { 1.09f, 1.59f, 1.79f }, { 0.013f, 0.070f, 0.145f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Spectralon", { 11.6f, 20.4f, 14.9f }, { 0.00f, 0.00f, 0.00f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Wholemilk", { 2.55f, 3.21f, 3.77f }, { 0.0011f, 0.0024f, 0.014f }, { 0.0f, 0.0f, 0.0f }, 1.3f },

	/* From "Acquiring Scattering Properties of Participating Media by Dilution"
	by Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen (SIGGRAPH 2006) */
	{ "Lowfat Milk", { 13.1157f, 15.4445f, 17.9572f }, { 0.00287f, 0.00575f, 0.01150f }, { 0.93200f, 0.90200f, 0.85900f }, 1.33f },
	{ "Reduced Milk", { 13.7335f, 15.6003f, 17.8007f }, { 0.00256f, 0.00511f, 0.01278f }, { 0.81900f, 0.79700f, 0.74600f }, 1.33f },
	{ "Regular Milk", { 18.2052f, 20.3826f, 22.3698f }, { 0.00153f, 0.00460f, 0.01993f }, { 0.75000f, 0.71400f, 0.68100f }, 1.33f },
	{ "Espresso", { 7.78262f, 8.13050f, 8.53875f }, { 4.79838f, 6.57512f, 8.84925f }, { 0.90700f, 0.89600f, 0.88000f }, 1.33f },
	{ "Mint Mocha Coffee", { 3.51133f, 4.14383f, 5.59667f }, { 3.77200f, 5.82283f, 7.82000f }, { 0.91000f, 0.90700f, 0.91400f }, 1.33f },
	{ "Lowfat Soy Milk", { 2.03838f, 2.32875f, 3.90281f }, { 0.00144f, 0.00719f, 0.03594f }, { 0.85000f, 0.85300f, 0.84200f }, 1.33f },
	{ "Regular Soy Milk", { 4.66325f, 5.20183f, 8.74575f }, { 0.00192f, 0.00958f, 0.06517f }, { 0.87300f, 0.85800f, 0.83200f }, 1.33f },
	{ "Lowfat Chocolate Milk", { 9.83710f, 11.4954f, 13.1629f }, { 0.01150f, 0.03680f, 0.15640f }, { 0.93400f, 0.92700f, 0.91600f }, 1.33f },
	{ "Regular Chocolate Milk", { 10.5685f, 13.1416f, 15.2202f }, { 0.01006f, 0.04313f, 0.14375f }, { 0.86200f, 0.83800f, 0.80600f }, 1.33f },
	{ "Coke", { 0.00254f, 0.00299f, 0.00000f }, { 0.10014f, 0.16503f, 0.24680f }, { 0.96500f, 0.97200f, 0.00000f }, 1.33f },
	{ "Pepsi", { 0.00083f, 0.00203f, 0.00000f }, { 0.09164f, 0.14158f, 0.20729f }, { 0.92600f, 0.97900f, 0.00000f }, 1.33f },
	{ "Sprite", { 0.00011f, 0.00014f, 0.00014f }, { 0.00189f, 0.00183f, 0.00200f }, { 0.94300f, 0.95300f, 0.95200f }, 1.33f },
	{ "Gatorade", { 0.03668f, 0.04488f, 0.05742f }, { 0.02479f, 0.01929f, 0.00888f }, { 0.93300f, 0.93300f, 0.93500f }, 1.33f },
	{ "Chardonnay", { 0.00021f, 0.00033f, 0.00048f }, { 0.01078f, 0.01186f, 0.02400f }, { 0.91400f, 0.95800f, 0.97500f }, 1.33f },
	{ "White Zinfandel", { 0.00022f, 0.00033f, 0.00046f }, { 0.01207f, 0.01618f, 0.01984f }, { 0.91900f, 0.94300f, 0.97200f }, 1.33f },
	{ "Merlot", { 0.00081f, 0.00000f, 0.00000f }, { 0.11632f, 0.25191f, 0.29434f }, { 0.97400f, 0.00000f, 0.00000f }, 1.33f },
	{ "Budweiser Beder", { 0.00029f, 0.00055f, 0.00059f }, { 0.01149f, 0.02491f, 0.05779f }, { 0.91700f, 0.95600f, 0.98200f }, 1.33f },
	{ "Coors Light Beer", { 0.00062f, 0.00127f, 0.00000f }, { 0.00616f, 0.01398f, 0.03498f }, { 0.91800f, 0.96600f, 0.00000f }, 1.33f },
	{ "Clorox", { 0.02731f, 0.03302f, 0.03695f }, { 0.00335f, 0.01489f, 0.02630f }, { 0.91200f, 0.90500f, 0.89200f }, 1.33f },
	{ "Apple Juice", { 0.00257f, 0.00311f, 0.00413f }, { 0.01296f, 0.02374f, 0.05218f }, { 0.94700f, 0.94900f, 0.94500f }, 1.33f },
	{ "Cranberry Juice", { 0.00196f, 0.00238f, 0.00301f }, { 0.03944f, 0.09422f, 0.12426f }, { 0.94700f, 0.95100f, 0.97400f }, 1.33f },
	{ "Grape Juice", { 0.00138f, 0.00000f, 0.00000f }, { 0.10404f, 0.23958f, 0.29325f }, { 0.96100f, 0.00000f, 0.00000f }, 1.33f },
	{ "Ruby Grapefruit Juice", { 0.15496f, 0.15391f, 0.15995f }, { 0.08587f, 0.18314f, 0.25262f }, { 0.92900f, 0.92900f, 0.93100f }, 1.33f },
	{ "White Grapefruit Juice", { 0.50499f, 0.52742f, 0.75282f }, { 0.01380f, 0.01883f, 0.05678f }, { 0.54800f, 0.54500f, 0.56500f }, 1.33f },
	{ "Shampoo", { 0.00797f, 0.00874f, 0.01127f }, { 0.01411f, 0.04569f, 0.06172f }, { 0.91000f, 0.90500f, 0.92000f }, 1.33f },
	{ "Strawberry Shampoo", { 0.00215f, 0.00245f, 0.00253f }, { 0.01449f, 0.05796f, 0.07582f }, { 0.92700f, 0.93500f, 0.99400f }, 1.33f },
	{ "Head & Shoulders Shampoo", { 0.26747f, 0.27696f, 0.29574f }, { 0.08462f, 0.15688f, 0.20365f }, { 0.91100f, 0.89600f, 0.88400f }, 1.33f },
	{ "Lemon Tea Powder", { 0.74489f, 0.83823f, 1.00158f }, { 2.42881f, 4.57573f, 7.21270f }, { 0.94600f, 0.94600f, 0.94900f }, 1.33f },
	{ "Orange Juice Powder", { 0.00193f, 0.00213f, 0.00226f }, { 0.00145f, 0.00344f, 0.00786f }, { 0.91900f, 0.91800f, 0.92200f }, 1.33f },
	{ "Pink Lemonade Powder", { 0.00123f, 0.00133f, 0.00131f }, { 0.00116f, 0.00237f, 0.00320f }, { 0.90200f, 0.90200f, 0.90400f }, 1.33f },
	{ "Cappuccino Powder", { 12.2094f, 16.4659f, 29.2727f }, { 35.8441f, 49.5470f, 61.0844f }, { 0.84900f, 0.84300f, 0.92600f }, 1.33f },
	{ "Salt Powder", { 0.13805f, 0.15677f, 0.17865f }, { 0.28415f, 0.32570f, 0.34148f }, { 0.80200f, 0.79300f, 0.82100f }, 1.33f },
	{ "Sugar Powder", { 0.00282f, 0.00315f, 0.00393f }, { 0.01264f, 0.03105f, 0.05012f }, { 0.92100f, 0.91900f, 0.93100f }, 1.33f },
	{ "Suisse Mocha Powder", { 30.0848f, 33.4452f, 38.7191f }, { 17.5020f, 27.0044f, 35.4334f }, { 0.90700f, 0.89400f, 0.88800f }, 1.33f },
	{ "Pacific Ocean Surface Water", { 0.00180f, 0.00183f, 0.00228f }, { 0.03184f, 0.03132f, 0.03015f }, { 0.90200f, 0.82500f, 0.91400f }, 1.33f },

	{ NULL, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, 0.0f }
};

#define SAMPLE_SCATTERING		0
#define SAMPLE_ABSORBING		1
#define SAMPLE_TRANSMITTANCE	2

#define SAMPLE_COEF	SAMPLE_ABSORBING

class ParticipatingMedia
{
		double cdf[3][3];
		double pdf[3][3];
		double total[3];

		double É–tSum[3];
		Color É–t_[3];
		int max_i[3], min_i[3], mid_i[3];

public:
	double scattering;			// É–sÅ@éUóê(scattering)
	double absorbing;			// É–aÅ@ãzé˚(absorption)
	double transmittance;		// É–t =É–s + É–a å∏êäåWêî(attenuation coefficient)

	Color scatteringColor;		// É–sÅ@éUóê(scattering)
	Color absorbingColor;		// É–aÅ@ãzé˚(absorption)
	Color transmittanceColor;	// É–t =É–s + É–a å∏êäåWêî(attenuation coefficient)
	Vector3d sigma_sp;
	Vector3d sigma_tp;

	double eta;					//ã¸ê‹ó¶
	Color albedoColor;
	double phase_prm[3];			// à ëää÷êîÉpÉâÉÅÅ[É^[Henyey-Greensteinä÷êî]

	double scattering_albedo;	// Scattering albedoÅBëäå›çÏópÇµÇΩÇ∆Ç´ÇÃÅAéUóêÇÃämó¶

	int use_brssdf;
	int level;					//-2:multi(bssrdf) -1:sigle+multi(bssrdf) 0Å`n:simulation(scattering=n)  
	double r_max;

	ParticipatingMedia()
	{
		r_max = 1.0;
		level = -1;
		use_brssdf = 1;
		eta = 0.0;
		phase_prm[0] = 0.0;
		phase_prm[1] = 0.0;
		phase_prm[2] = 0.0;
		scattering = 0.0;
		absorbing = 0.0;
		transmittance = 0.0;
		scatteringColor = Vector3d(0.0,0.0,0.0);
		absorbingColor  = Vector3d(0.0,0.0,0.0);;
		albedoColor = Vector3d(0.0, 0.0, 0.0);;
	}
	inline void setup_Scattering(const double material_scale=1.0)
	{
		scatteringColor = scatteringColor*material_scale;
		absorbingColor  = absorbingColor*material_scale;

		//É–t =É–s + É–a å∏êäåWêî(attenuation coefficient)
		scattering = (scatteringColor.x + scatteringColor.y + scatteringColor.z)/3,0;
		absorbing = (absorbingColor.x + absorbingColor.y + absorbingColor.z)/3,0;

		transmittance = scattering + absorbing;
		transmittanceColor = scatteringColor + absorbingColor;

		if (transmittance < 1.0e-12)
		{
			scattering = 0;
			absorbing = 0;
			transmittance = 0.0;
			scattering_albedo = 0.0;
			scatteringColor = Vector3d(1.0, 1.0, 1.0);
			absorbingColor = Vector3d(1.0, 1.0, 1.0);
			albedoColor = Vector3d(1.0, 1.0, 1.0);
		}
		else
		{
			scattering_albedo = scattering / transmittance; // Scattering albedoÅBëäå›çÏópÇµÇΩÇ∆Ç´ÇÃÅAéUóêÇÃämó¶
			albedoColor = scatteringColor / transmittanceColor;
		}
		sigma_sp = scatteringColor * (Vector3d(1.0, 1.0, 1.0) - Vector3d(phase_prm));
		sigma_tp = sigma_sp + absorbingColor;

#if USE_PARTICIPATING_MEDIA
		printf("éUóê(scattering %f) %f %f %f\n", scattering, scatteringColor.x, scatteringColor.y, scatteringColor.z);
		printf("ãzé˚(absorption %f) %f %f %f\n", absorbing, absorbingColor.x, absorbingColor.y, absorbingColor.z);
		printf("å∏êäåWêî(attenuation coefficient %f) %f %f %f\n", transmittance, transmittanceColor.x, transmittanceColor.y, transmittanceColor.z);
		printf("ÉAÉãÉxÉh(albedo coefficient %f) %f %f %f\n", scattering_albedo, albedoColor.x, albedoColor.y, albedoColor.z);
		printf("à ëää÷êîÉpÉâÉÅÅ[É^[Henyey-Greensteinä÷êî] %f %f %f\n", phase_prm[0], phase_prm[1], phase_prm[2]);
#endif
	
		Color scattering_absorbing[3] = {scatteringColor, absorbingColor, transmittanceColor};
		for ( int i = 0; i < 3; i++ )
		{
			max_i[i] = -1;
			min_i[i] = -1;
		}

		for ( int k = 0; k < 3; k++ )
		{
			double max_s = -9999999.0;
			double min_s = 9999999.0;

			total[k] = 0.0;
			for (int i = 0; i < 3; i ++) 
			{
				if (max_s < scattering_absorbing[k][i])
				{
					max_s = scattering_absorbing[k][i];
					max_i[k] = i;
				}
				if (min_s > absorbingColor[i])
				{
					min_s = absorbingColor[i];
					min_i[k] = i;
				}
			}
			fprintf(stderr, "max_s %f min_s %f\n", scattering_absorbing[k][max_i[k]], scattering_absorbing[k][min_i[k]]);
			mid_i[k] = -1;
			for (int i = 0; i < 3; i ++) 
			{
				if (max_s > scattering_absorbing[k][i] && min_s < scattering_absorbing[k][i])
				{
					mid_i[k] = i; break;
				}
			}
			if (mid_i[k] == -1)
			{
				double av = absorbingColor.Av();

				if (fabs(scattering_absorbing[k].Max() - av) < fabs(scattering_absorbing[k].Min() - av))
				{
					mid_i[k] = max_i[k];
				}
				else
				{
					mid_i[k] = min_i[k];
				}
			}
			fprintf(stderr, "mid %f %d\n", scattering_absorbing[k][mid_i[k]], mid_i[k]);

			É–tSum[k] = scattering_absorbing[k][0] + scattering_absorbing[k][1] + scattering_absorbing[k][2];

			total[k] = 0.0;
			if ( É–tSum[k] > 1.0e-14 )
			{
				É–t_[k] = Color(
					scattering_absorbing[k][min_i[k]] / É–tSum[k],
					scattering_absorbing[k][mid_i[k]] / É–tSum[k],
					scattering_absorbing[k][max_i[k]] / É–tSum[k]);

				for (int i = 0; i < 3; i ++) 
				{
					total[k] += É–t_[k][i];
					cdf[k][i] = total[k];
					pdf[k][i] = É–t_[k][i];
				}

				for (int i = 0; i < 3; i ++) 
				{
					cdf[k][i] /= total[k];
					pdf[k][i] /= total[k];
				}

			}else
			{
				É–t_[k] = Color(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
			}
		}
	}

	inline int sample(const double u, const int k, const int inv) const
	{
		if ( É–tSum[k] < 1.0e-14 )
		{
			return std::min(2, (int)(3.0*u));
		}

		int ii = 0;
		double s = 1.0;
		if (inv)
		{
			ii = min_i[k];
			s = -1.0;
		}

		if (u < É–t_[k][0]) ii = min_i[k];
		else if (u < É–t_[k][0] + É–t_[k][1]) ii = mid_i[k];
		else if (u < É–t_[k][0] + É–t_[k][1] + É–t_[k][2]) ii = max_i[k];

		if (inv)
		{
			if (-u < -É–t_[k][0]) ii = min_i[k];
			else if (-u < -É–t_[k][0] - É–t_[k][1]) ii = mid_i[k];
			else if (-u < -É–t_[k][0] - É–t_[k][1] - É–t_[k][2]) ii = max_i[k];
		}

		//fprintf(stderr, "%d %d %d -> %d\n", min_i, mid_i, max_i, ii);
		return ii;
	}

	inline int sample(const double u, const int inv) const
	{
		return sample(u, SAMPLE_COEF, inv);
	}
	inline int sample_scattering(const double u, const int inv) const
	{
		return sample(u, 0, inv);
	}
	inline int sample_absorbing(const double u, const int inv) const
	{
		return sample(u, 1, inv);
	}
	inline int sample_transmittance(const double u, const int inv) const
	{
		return sample(u, 2, inv);
	}

	inline int sample( double& pdf_, const double u, const int k ) const
	{
		for (int i = 0; i < 3; i ++)
		{
			if ( u <= cdf[k][i] )
			{
				pdf_ = pdf[k][i];
				return i;
			}
		}
		pdf_ = pdf[k][0];
		return 0;
	}

#if 0
	inline double scatteringRatio(const double d) const
	{
		return exp(-scattering * d );
	}
	inline double absorbingRatio(const double d) const
	{
		return exp(-absorbing * d );
	}
	inline double transmittanceRatio(const double d) const
	{
		return exp(-transmittance * d );
	}
#else
	inline Spectrum scatteringRatio(const double d, const double wavelength) const
	{
		//return RGB2Spectrum(Exp(-1.0*Vector3d(1,1,1)*scattering * d ), wavelength);
		return RGB2Spectrum(Exp(-1.0*scatteringColor * d ), wavelength);
	}
	inline Spectrum absorbingRatio(const double d, const double wavelength) const
	{
		//return RGB2Spectrum(Exp(-1.0*Vector3d(1,1,1)*absorbing * d ), wavelength);
		return RGB2Spectrum(Exp(-1.0*absorbingColor * d ), wavelength);
	}
	inline Spectrum transmittanceRatio(const double d, const double wavelength) const
	{
		//return RGB2Spectrum(Exp(-1.0*Vector3d(1,1,1)*transmittance * d ), wavelength);
		return RGB2Spectrum(Exp(-1.0*transmittanceColor * d ), wavelength);
	}
#endif

};

};

#endif
