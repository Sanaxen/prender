#ifndef __BITMAP_H_INCLUDED__
#define __BITMAP_H_INCLUDED__

#define FILEHEADERSIZE 14
#define INFOHEADERSIZE 40
#define HEADERSIZE (FILEHEADERSIZE+INFOHEADERSIZE)

#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include<math.h>

#include <omp.h>

//#define BMP_USE_GPU

//#include <amp.h>
//#include <amp_math.h>
//using namespace concurrency;
//#ifdef BMP_USE_GPU
//#define BMP_ACC_RESTRICTION restrict(amp)
//#define ACC_TYPE "GPU"
//#define FAST_MATH	fast_math::
//#else if
//#define BMP_ACC_RESTRICTION restrict(cpu)
//#define ACC_TYPE "CPU"
//#endif

#pragma warning(disable:4996)
#pragma warning(disable:4018)
#pragma warning(disable:4101)
#pragma warning(disable:4244)

//#define OMP_SCHEDULE_BMP schedule(guided)
#define OMP_SCHEDULE_BMP schedule(static)
#define LINELENGMAX	(4096*10)

typedef struct Rgb_{
	unsigned char b;
	unsigned char g;
	unsigned char r;
	unsigned char alp;
	~ Rgb_(){}
	Rgb_(){}
	inline Rgb_(int x, int y, int z)
	{
		r = x;
		g = y;
		b = z;
		alp = 255;
	}
	inline Rgb_(const int* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
		alp = 255;
	}
	inline Rgb_(const unsigned char* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
		alp = 255;
	}
	inline Rgb_(const unsigned char x)
	{
		r = x;
		g = x;
		b = x;
		alp = 255;
	}
} Rgb;

typedef struct iRgb_{
	unsigned int b;
	unsigned int g;
	unsigned int r;
	unsigned char alp;
	~ iRgb_(){}
	iRgb_(){}
	inline iRgb_(int x, int y, int z)
	{
		r = x;
		g = y;
		b = z;
		alp = 255;
	}
	inline iRgb_(const int* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
		alp = 255;
	}
	inline iRgb_(const unsigned int* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
		alp = 255;
	}
} iRgb;


typedef struct{
	unsigned int height;
	unsigned int width;
	Rgb *data;
}Image;


class BitMap
{
	Image* data;
	Image *Read_Bmp(char *filename);
	void Free_Image(Image *img);
	int Write_Bmp(char *filename, Image *img);
	Image *Create_Image(int width, int height);
	int WriteText(char* filename, Image* img);
	Image* Read_Text(char *filename);
	Image* Read_Csv(char* filename, double min, double max);
	int Write_Csv(char* filename, Image* img, int rgb);
	int Write_Csv(char* filename, Image* img, int rgb, double min, double max);

public:

	inline BitMap()
	{
		data = NULL;
	}

	inline ~BitMap()
	{
		Clear();
	}

	inline void Clear()
	{
		if ( data ) Free_Image(data);
		data = NULL;
	}

	void Create(int width, int height)
	{
		data = Create_Image(width, height);
	}

	inline Image* GetImage() const
	{ 
		return data;
	}
	inline int W() const
	{
		return data->width;
	}
	inline int H() const
	{
		return data->height;
	}
	inline void Copy(BitMap& bmp)
	{
		Clear();
		Create(bmp.data->width, bmp.data->height);
		memcpy(data->data, bmp.data->data, sizeof(Rgb)*bmp.data->width*bmp.data->height);
	}

	void Write(char *filename)
	{
		Write_Bmp(filename, data);
	}
	void Read(char *filename)
	{
		data = Read_Bmp( filename);
	}

	void WriteText(char *filename)
	{
		WriteText( filename, data);
	}

	void ReadText(char *filename)
	{
		Read_Text(filename);
	}

	inline Rgb& cell(const int i, const int j) const
	{
		return *(data->data+((data->height-i-1)*data->width + j));
	}

	void ToGrayScale();

	void ReadCsv(char* filename, double min, double max)
	{
		if ( data ) Free_Image(data);
		data = Read_Csv(filename, min, max);
	}
	void ReadCsv(double* value, int w, int h, double min, double max)
	{
		if ( data ) Free_Image(data);
		data = Read_Csv(value, w, h, min, max);
	}

	void WriteCsv(char* filename, int rgb)
	{
		Write_Csv( filename, data, rgb);
	}
	void WriteCsv(char* filename, int rgb, double min, double max)
	{
		Write_Csv( filename, data, rgb, min, max);
	}

	void Reverse();
	void ToGrayScale_and_Reverse();

	void Offset(int size);
	void convolve_smooth(int* mask, double mat[3][3]);
	
	void convolve_smooth(int* mask)
	{
		double conv[3][3]={{1.0, 1.0, 1.0},
		{1.0, 1.0, 1.0},
		{1.0, 1.0, 1.0}};
		convolve_smooth(mask, conv);
	}
	
	Image* Read_Csv(double* value, int w, int h, double min, double max);


	static int colortableNum;
	static unsigned char colorTbl[1024][3];
	void ColorTable();
	void ColorTable(BitMap& colormap);
	void ColorTable(int startIndex, int endIndex,unsigned char start[3], unsigned char end[3]);
	void ColorLevel( double min, double max, double* z, double zmask, unsigned char* maskcolor=NULL, int* top=NULL, double* elv=NULL);

};

void filter_non_local_means(float* image, const int width, const int height, float param_h, float sigma);

#endif /*__BITMAP_H_INCLUDED__*/
