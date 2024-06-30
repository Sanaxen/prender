#ifndef	_Vector4d_H_
#define	_Vector4d_H_

#include <cmath>

namespace prender {

struct Vector4d {
	double x, y, z, w;
	inline Vector4d(const double x = 0, const double y = 0, const double z = 0, double w = 0) : x(x), y(y), z(z), w(w) {}
	
	inline Vector4d(const double* xx)
	{
		x = xx[0];
		y = xx[1];
		z = xx[2];
		w = xx[3];
	}
	inline Vector4d(const float* xx)
	{
		x = xx[0];
		y = xx[1];
		z = xx[2];
		w = xx[3];
	}
	inline Vector4d operator+(const Vector4d &b) const 
	{
		return Vector4d(x + b.x, y + b.y, z + b.z, w + b.w);
	}
	
	inline Vector4d operator-(const Vector4d &b) const 
	{
		return Vector4d(x - b.x, y - b.y, z - b.z, w - b.w);
	}

	inline Vector4d operator*(const double b) const
	{
		return Vector4d(x * b, y * b, z * b, w * b);
	}
	
	inline Vector4d operator/(const double b) const
	{
		return Vector4d(x / b, y / b, z / b, w / b);
	}
	
	inline const double sqr() const
	{ 
		return x*x + y*y + z*z + w*w; 
	}
	
	inline const double length() const 
	{ 
		return sqrt(sqr()); 
	}
};

inline Vector4d operator*(double f, const Vector4d &v)
{ 
	return v * f; 
}

inline Vector4d normalize(const Vector4d &v)
{
	return v * (1.0 / v.length()); 
}

inline const Vector4d multiply(const Vector4d &v1, const Vector4d &v2)
{
	return Vector4d(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z, v1.w * v2.w);
}

inline const double dot(const Vector4d &v1, const Vector4d &v2) 
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}


};

#endif
