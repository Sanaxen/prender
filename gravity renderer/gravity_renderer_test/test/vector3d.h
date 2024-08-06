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

#define FOR_LOOP4	for ( int i = 0; i < 4; i++ )
#define FOR_LOOP16	for ( int i = 0; i < 16; i++ )

#define FOR_LOOP3	for ( int i = 0; i < 3; i++ )

class Matrix4D
{
public:
	double val[16];

	inline Matrix4D()
	{
		FOR_LOOP16 val[i] = 0;
	}
	inline Matrix4D(
		double x00, double x01, double x02, double x03,
		double x10, double x11, double x12, double x13,
		double x20, double x21, double x22, double x23,
		double x30, double x31, double x32, double x33
		)
	{
		val[0] = x00; val[1] = x01; val[2] = x02; val[2] = x03;
		val[3] = x10; val[4] = x11; val[5] = x12; val[6] = x13;
		val[7] = x20; val[8] = x21; val[9] = x22; val[10] = x23;
		val[11] = x30; val[12] = x31; val[13] = x32; val[14] = x33;
	}
	inline Matrix4D(double m[16])
	{
		memcpy(val, m, 16 * sizeof(double));
	}
	inline void getVal(double m[16])
	{
		memcpy(m, val, 16 * sizeof(double));
	}

	inline Matrix4D& operator=(const Matrix4D& A)
	{
		FOR_LOOP16 this->val[i] = A.val[i];
		return *this;
	}
	inline Matrix4D operator+()
	{
		return *this;
	}
	inline Matrix4D operator-()
	{
		Matrix4D A;
		FOR_LOOP16 A.val[i] = -val[i];
		return A;
	}
	inline Matrix4D& operator+=(const Matrix4D& A)
	{
		FOR_LOOP4 val[i] += A.val[i];
		return *this;
	}
	inline Matrix4D& operator-=(const Matrix4D& A)
	{
		FOR_LOOP4 val[i] -= A.val[i];
		return *this;
	}
	inline Matrix4D& operator*=(double k)
	{
		FOR_LOOP4 val[i] *= k;
		return *this;
	}
	inline Matrix4D& operator/=(double k)
	{
		FOR_LOOP4 val[i] /= k;
		return *this;
	}

	inline bool operator==(const Matrix4D& A) const
	{
		bool result = true;
		FOR_LOOP4 result &= (val[i] == A.val[i]);
		return result;
	}
	inline bool operator!=(const Matrix4D& A) const
	{
		return !((*this) == A);
	}
	inline double& operator()(int i, int j)
	{
		return val[4 * i + j];
	}
	inline void LoadIdentity()
	{
		size_t	s = sizeof(val);
		memset(val, 0, s);
		val[0] = 1.0;
		val[5] = 1.0;
		val[10] = 1.0;
		val[15] = 1.0;
	}
	inline double Trace(const Matrix4D& A) const
	{
		double tr = 0;
		FOR_LOOP4 tr += A.val[4 * i + i];
		return tr;
	}

	inline Matrix4D Transpose(const Matrix4D& A) const
	{
		Matrix4D AT;
		for (int i = 0; i < 4; i++)for (int j = 0; j < 4; j++) AT.val[4 * i + j] = A.val[4 * j + i];
		return AT;
	}

	//拡大・縮小マトリクスの生成
	Matrix4D& Scale(const double scale_x, const double scale_y, const double scale_z)
	{
		Matrix4D&	m = *this;
		LoadIdentity();
		m(0, 0) = scale_x;
		m(1, 1) = scale_y;
		m(2, 2) = scale_z;
		return *this;
	}

	//並行移動マトリクスの生成
	Matrix4D& Translation(const double pos_x, const double pos_y, const double pos_z)
	{
		Matrix4D&	m = *this;
		LoadIdentity();
		m(3, 0) = pos_x;
		m(3, 1) = pos_y;
		m(3, 2) = pos_z;
		return *this;
	}

	//各軸の回転マトリクスを生成
	Matrix4D& RotationX(const double angle)
	{
		Matrix4D&	m = *this;
		LoadIdentity();
		m(1, 1) =
			m(2, 2) = (cos(angle));
		m(1, 2) = (sin(angle));
		m(2, 1) = -m(1, 2);
		return *this;
	}
	Matrix4D& RotationY(const double angle)
	{
		Matrix4D&	m = *this;
		LoadIdentity();
		m(0, 0) =
			m(2, 2) = (cos(angle));
		m(2, 0) = (sin(angle));
		m(0, 2) = -m(2, 0);
		return *this;
	}
	Matrix4D& RotationZ(const double angle)
	{
		Matrix4D&	m = *this;
		LoadIdentity();
		m(0, 0) =
			m(1, 1) = (cos(angle));
		m(0, 1) = (sin(angle));
		m(1, 0) = -m(0, 1);
		return *this;
	}

	//任意軸の回転マトリクス生成
	Matrix4D& AxisRotation(Vector3d& axis, const double angle)
	{
		Matrix4D&	m = *this;
		Vector3d	v = normalize(axis);
		double	s = (sin(angle));
		double	c = (cos(angle));
		m(0, 0) = v.x * v.x * (1.0 - c) + c;
		m(0, 1) = v.x * v.y * (1.0 - c) - v.z * s;
		m(0, 2) = v.z * v.x * (1.0 - c) + v.y * s;
		m(0, 3) = 0.0;

		m(1, 0) = v.x * v.y * (1.0 - c) + v.z * s;
		m(1, 1) = v.y * v.y * (1.0 - c) + c;
		m(1, 2) = v.y * v.z * (1.0 - c) - v.x * s;
		m(1, 3) = 0.0;

		m(2, 0) = v.z * v.x * (1.0 - c) - v.y * s;
		m(2, 1) = v.y * v.z * (1.0 - c) + v.x * s;
		m(2, 2) = v.z * v.z * (1.0 - c) + c;
		m(2, 3) = 0.0;

		m(3, 0) = 0.0;
		m(3, 1) = 0.0;
		m(3, 2) = 0.0;
		m(3, 3) = 1.0;

		return *this;
	}


	inline Matrix4D InvertMatrix();
};

inline  Matrix4D operator+(const Matrix4D& A, const Matrix4D& B)
{
	 Matrix4D C;
	FOR_LOOP16 C.val[i] = A.val[i] + B.val[i];
	return C;
}
 inline  Matrix4D operator-(const Matrix4D& A, const Matrix4D& B)
{
	 Matrix4D C;
	FOR_LOOP16 C.val[i] = A.val[i] - B.val[i];
	return C;
}

 inline  Matrix4D operator*(double k, const  Matrix4D& A)
{
	 Matrix4D B;
	FOR_LOOP16 B.val[i] = A.val[i] * k;
	return B;
}
 inline  Matrix4D operator*(const Matrix4D& A, double k)
{
	 Matrix4D B;
	FOR_LOOP16 B.val[i] = A.val[i] * k;
	return B;
}
 inline Matrix4D operator/(const Matrix4D& A, double k)
{
	 Matrix4D B;
	FOR_LOOP16 B.val[i] = A.val[i] / k;
	return B;
}
inline Matrix4D operator*(const Matrix4D& A, const Matrix4D& B)
{
	 Matrix4D C;
	C.val[0] = A.val[0] * B.val[0] + A.val[4] * B.val[1] + A.val[8] * B.val[2] + A.val[12] * B.val[3];
	C.val[1] = A.val[1] * B.val[0] + A.val[5] * B.val[1] + A.val[9] * B.val[2] + A.val[13] * B.val[3];
	C.val[2] = A.val[2] * B.val[0] + A.val[6] * B.val[1] + A.val[10] * B.val[2] + A.val[14] * B.val[3];
	C.val[3] = A.val[3] * B.val[0] + A.val[7] * B.val[1] + A.val[11] * B.val[2] + A.val[15] * B.val[3];
	C.val[4] = A.val[0] * B.val[4] + A.val[4] * B.val[5] + A.val[8] * B.val[6] + A.val[12] * B.val[7];
	C.val[5] = A.val[1] * B.val[4] + A.val[5] * B.val[5] + A.val[9] * B.val[6] + A.val[13] * B.val[7];
	C.val[6] = A.val[2] * B.val[4] + A.val[6] * B.val[5] + A.val[10] * B.val[6] + A.val[14] * B.val[7];
	C.val[7] = A.val[3] * B.val[4] + A.val[7] * B.val[5] + A.val[11] * B.val[6] + A.val[15] * B.val[7];
	C.val[8] = A.val[0] * B.val[8] + A.val[4] * B.val[9] + A.val[8] * B.val[10] + A.val[12] * B.val[11];
	C.val[9] = A.val[1] * B.val[8] + A.val[5] * B.val[9] + A.val[9] * B.val[10] + A.val[13] * B.val[11];
	C.val[10] = A.val[2] * B.val[8] + A.val[6] * B.val[9] + A.val[10] * B.val[10] + A.val[14] * B.val[11];
	C.val[11] = A.val[3] * B.val[8] + A.val[7] * B.val[9] + A.val[11] * B.val[10] + A.val[15] * B.val[11];
	C.val[12] = A.val[0] * B.val[12] + A.val[4] * B.val[13] + A.val[8] * B.val[14] + A.val[12] * B.val[15];
	C.val[13] = A.val[1] * B.val[12] + A.val[5] * B.val[13] + A.val[9] * B.val[14] + A.val[13] * B.val[15];
	C.val[14] = A.val[2] * B.val[12] + A.val[6] * B.val[13] + A.val[10] * B.val[14] + A.val[14] * B.val[15];
	C.val[15] = A.val[3] * B.val[12] + A.val[7] * B.val[13] + A.val[11] * B.val[14] + A.val[15] * B.val[15];
	return C;
}

 inline bool InvertMatrix_(const Matrix4D& A, Matrix4D& B, double& det)
{
	double inv[16];

	det = 0.0;

	inv[0] = A.val[5] * A.val[10] * A.val[15] -
		A.val[5] * A.val[11] * A.val[14] -
		A.val[9] * A.val[6] * A.val[15] +
		A.val[9] * A.val[7] * A.val[14] +
		A.val[13] * A.val[6] * A.val[11] -
		A.val[13] * A.val[7] * A.val[10];

	inv[4] = -A.val[4] * A.val[10] * A.val[15] +
		A.val[4] * A.val[11] * A.val[14] +
		A.val[8] * A.val[6] * A.val[15] -
		A.val[8] * A.val[7] * A.val[14] -
		A.val[12] * A.val[6] * A.val[11] +
		A.val[12] * A.val[7] * A.val[10];

	inv[8] = A.val[4] * A.val[9] * A.val[15] -
		A.val[4] * A.val[11] * A.val[13] -
		A.val[8] * A.val[5] * A.val[15] +
		A.val[8] * A.val[7] * A.val[13] +
		A.val[12] * A.val[5] * A.val[11] -
		A.val[12] * A.val[7] * A.val[9];

	inv[12] = -A.val[4] * A.val[9] * A.val[14] +
		A.val[4] * A.val[10] * A.val[13] +
		A.val[8] * A.val[5] * A.val[14] -
		A.val[8] * A.val[6] * A.val[13] -
		A.val[12] * A.val[5] * A.val[10] +
		A.val[12] * A.val[6] * A.val[9];

	inv[1] = -A.val[1] * A.val[10] * A.val[15] +
		A.val[1] * A.val[11] * A.val[14] +
		A.val[9] * A.val[2] * A.val[15] -
		A.val[9] * A.val[3] * A.val[14] -
		A.val[13] * A.val[2] * A.val[11] +
		A.val[13] * A.val[3] * A.val[10];

	inv[5] = A.val[0] * A.val[10] * A.val[15] -
		A.val[0] * A.val[11] * A.val[14] -
		A.val[8] * A.val[2] * A.val[15] +
		A.val[8] * A.val[3] * A.val[14] +
		A.val[12] * A.val[2] * A.val[11] -
		A.val[12] * A.val[3] * A.val[10];

	inv[9] = -A.val[0] * A.val[9] * A.val[15] +
		A.val[0] * A.val[11] * A.val[13] +
		A.val[8] * A.val[1] * A.val[15] -
		A.val[8] * A.val[3] * A.val[13] -
		A.val[12] * A.val[1] * A.val[11] +
		A.val[12] * A.val[3] * A.val[9];

	inv[13] = A.val[0] * A.val[9] * A.val[14] -
		A.val[0] * A.val[10] * A.val[13] -
		A.val[8] * A.val[1] * A.val[14] +
		A.val[8] * A.val[2] * A.val[13] +
		A.val[12] * A.val[1] * A.val[10] -
		A.val[12] * A.val[2] * A.val[9];

	inv[2] = A.val[1] * A.val[6] * A.val[15] -
		A.val[1] * A.val[7] * A.val[14] -
		A.val[5] * A.val[2] * A.val[15] +
		A.val[5] * A.val[3] * A.val[14] +
		A.val[13] * A.val[2] * A.val[7] -
		A.val[13] * A.val[3] * A.val[6];

	inv[6] = -A.val[0] * A.val[6] * A.val[15] +
		A.val[0] * A.val[7] * A.val[14] +
		A.val[4] * A.val[2] * A.val[15] -
		A.val[4] * A.val[3] * A.val[14] -
		A.val[12] * A.val[2] * A.val[7] +
		A.val[12] * A.val[3] * A.val[6];

	inv[10] = A.val[0] * A.val[5] * A.val[15] -
		A.val[0] * A.val[7] * A.val[13] -
		A.val[4] * A.val[1] * A.val[15] +
		A.val[4] * A.val[3] * A.val[13] +
		A.val[12] * A.val[1] * A.val[7] -
		A.val[12] * A.val[3] * A.val[5];

	inv[14] = -A.val[0] * A.val[5] * A.val[14] +
		A.val[0] * A.val[6] * A.val[13] +
		A.val[4] * A.val[1] * A.val[14] -
		A.val[4] * A.val[2] * A.val[13] -
		A.val[12] * A.val[1] * A.val[6] +
		A.val[12] * A.val[2] * A.val[5];

	inv[3] = -A.val[1] * A.val[6] * A.val[11] +
		A.val[1] * A.val[7] * A.val[10] +
		A.val[5] * A.val[2] * A.val[11] -
		A.val[5] * A.val[3] * A.val[10] -
		A.val[9] * A.val[2] * A.val[7] +
		A.val[9] * A.val[3] * A.val[6];

	inv[7] = A.val[0] * A.val[6] * A.val[11] -
		A.val[0] * A.val[7] * A.val[10] -
		A.val[4] * A.val[2] * A.val[11] +
		A.val[4] * A.val[3] * A.val[10] +
		A.val[8] * A.val[2] * A.val[7] -
		A.val[8] * A.val[3] * A.val[6];

	inv[11] = -A.val[0] * A.val[5] * A.val[11] +
		A.val[0] * A.val[7] * A.val[9] +
		A.val[4] * A.val[1] * A.val[11] -
		A.val[4] * A.val[3] * A.val[9] -
		A.val[8] * A.val[1] * A.val[7] +
		A.val[8] * A.val[3] * A.val[5];

	inv[15] = A.val[0] * A.val[5] * A.val[10] -
		A.val[0] * A.val[6] * A.val[9] -
		A.val[4] * A.val[1] * A.val[10] +
		A.val[4] * A.val[2] * A.val[9] +
		A.val[8] * A.val[1] * A.val[6] -
		A.val[8] * A.val[2] * A.val[5];

	det = A.val[0] * inv[0] + A.val[1] * inv[4] + A.val[2] * inv[8] + A.val[3] * inv[12];
	if (det == 0.0)  return false;

	det = 1.0 / det;

	FOR_LOOP16 B.val[i] = inv[i] * det;

	return true;
}

Matrix4D Matrix4D::InvertMatrix()
{
	Matrix4D B;
	double det;

	bool res = InvertMatrix_(*this, B, det);
	if (!res)
	{
		printf("InvertMatrix_ Error.\n");
	}

	return B;
}

inline Vector3d operator * (Matrix4D mtx, Vector3d vec)
{
	Vector3d	ret;
	ret.x = vec.x * mtx(0, 0) + vec.y * mtx(1, 0) + vec.z * mtx(2, 0) + mtx(3, 0);
	ret.y = vec.x * mtx(0, 1) + vec.y * mtx(1, 1) + vec.z * mtx(2, 1) + mtx(3, 1);
	ret.z = vec.x * mtx(0, 2) + vec.y * mtx(1, 2) + vec.z * mtx(2, 2) + mtx(3, 2);
	return ret;

}
};

#endif
