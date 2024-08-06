#ifndef	_Vector3d_H_
#define	_Vector3d_H_

#include <cmath>
#include <algorithm>
#include "constant.h"

namespace prender {

struct Vector3d {
	double x, y, z;
	inline Vector3d(const double x = 0, const double y = 0, const double z = 0) : x(x), y(y), z(z) {}
	inline Vector3d(const double* xx)
	{
		x = xx[0];
		y = xx[1];
		z = xx[2];
	}
	inline Vector3d(const float* xx)
	{
		x = xx[0];
		y = xx[1];
		z = xx[2];
	}

	inline double& operator[](  int i)
	{
		if (i == 0) return x;
		if (i == 1) return y;
		if (i == 2) return z;
	}
	inline double operator[]( int i) const
	{
		if (i == 0) return x;
		if (i == 1) return y;
		if (i == 2) return z;
	}
	inline Vector3d operator+(const Vector3d &b) const 
	{
		return Vector3d(x + b.x, y + b.y, z + b.z);
	}
	
	inline Vector3d operator-(const Vector3d &b) const 
	{
		return Vector3d(x - b.x, y - b.y, z - b.z);
	}

	inline Vector3d operator*(const double b) const
	{
		return Vector3d(x * b, y * b, z * b);
	}
	
	inline Vector3d operator/(const double b) const
	{
		return Vector3d(x / b, y / b, z / b);
	}
	
	inline Vector3d operator/(const Vector3d& v) const
	{
		return Vector3d(x / v.x, y / v.y, z / v.z);
	}

	inline const double sqr() const
	{ 
		return x*x + y*y + z*z; 
	}
	
	inline const double length() const 
	{ 
		return sqrt(sqr()); 
	}
	inline double Max() const { return max(max(x, y), z); }
	inline double Min() const { return min(min(x, y), z); }

	inline friend Vector3d operator+(double b, const Vector3d& v) { return Vector3d(b + v.x, b + v.y, b + v.z); }
	inline friend Vector3d operator-(double b, const Vector3d& v) { return Vector3d(b - v.x, b - v.y, b - v.z); }
	inline friend Vector3d operator*(double b, const Vector3d& v) { return Vector3d(b * v.x, b * v.y, b * v.z); }
	inline friend Vector3d operator/(double b, const Vector3d& v) { return Vector3d(b / v.x, b / v.y, b / v.z); }

};

//inline Vector3d operator*(double f, const Vector3d &v)
//{ 
//	return v * f; 
//}

inline Vector3d normalize(const Vector3d &v)
{
	return v * (1.0 / v.length()); 
}

inline const Vector3d multiply(const Vector3d &v1, const Vector3d &v2)
{
	return Vector3d(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline const Vector3d operator*(const Vector3d &v1, const Vector3d &v2)
{
	return multiply(v1, v2);
}

inline const double dot(const Vector3d &v1, const Vector3d &v2) 
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline const Vector3d cross(const Vector3d &v1, const Vector3d &v2)
{
	return Vector3d(
		(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector3d operator%(const Vector3d &v1, const Vector3d &v2)
{ 
	return cross(v1, v2); 
}
inline double operator^(const Vector3d &v1, const Vector3d &v2)
{ 
	return dot(v1, v2); 
}

inline double Clamp(double val, double low, double high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


//正規直交基底
inline void OrthonormalBasis(const Vector3d& orienting_normal, Vector3d& w, Vector3d& u, Vector3d& v)
{
		// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
		w = orienting_normal;
		if (fabs(w.x) > PS_EPS) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
			u = normalize(cross(Vector3d(0.0, 1.0, 0.0), w));
		else
			u = normalize(cross(Vector3d(1.0, 0.0, 0.0), w));
		v = normalize(cross(w, u));
}

inline Vector3d SphericalDirection(const double sintheta, const double costheta, const double phi, const Vector3d &x, const Vector3d &y, const Vector3d &z) 
{
	return sintheta * cos(phi) * x + sintheta * sin(phi) * y + costheta * z;
}
inline double SphericalTheta(const Vector3d &v) {
    return acos(Clamp(v.z, -1.0, 1.0));
}
inline double SphericalPhi(const Vector3d &v) {
    double p = atan2(v.y, v.x);
	return (p < 0.0) ? p + 2.0*PS_PI : p;
}


inline Vector3d Exp(const Vector3d& v )
{
	return Vector3d( exp(v.x), exp(v.y), exp(v.z) );
}
inline Vector3d Sqrt(const Vector3d& v )
{
	return Vector3d( sqrt(v.x), sqrt(v.y), sqrt(v.z) );
}
inline Vector3d Sin(const Vector3d& v )
{
	return Vector3d( sin(v.x), sin(v.y), sin(v.z) );
}
inline Vector3d Cos(const Vector3d& v )
{
	return Vector3d( cos(v.x), cos(v.y), cos(v.z) );
}

};

#endif
