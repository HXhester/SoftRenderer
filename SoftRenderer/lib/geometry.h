#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
class Matrix;

template<class t> struct Vec2 {
	union 
	{
		struct { t x, y; };
		struct { t u, v; };

		t raw[2];
	};

	Vec2() : u(0), v(0) {}
	Vec2(t _u, t _v) :u(_u), v(_v) {}

	inline Vec2<t> operator +(const Vec2<t>& V) const { return Vec2<t>(u + V.u, v + V.v); }
	inline Vec2<t> operator -(const Vec2<t>& V) const { return Vec2<t>(u - V.u, v - V.v); }
	inline Vec2<t> operator *(float f) const { return Vec2<t>(f * u, f * v); }
	inline t operator *(const Vec2<t>& V) const { return u * V.u + v * V.v; }
	inline t operator [](const int i) const { assert(i >= 0 && i < 2); return raw[i]; }	//read
	inline t &operator [](const int i)  { assert(i >= 0 && i < 2); return raw[i]; }		//write
	template <class > friend std::ostream& operator <<(std::ostream& s, Vec2<t>& v);
};

template<class t> struct Vec3 {
	union
	{
		struct { t x, y, z; };

		t raw[3];
	};

	Vec3() : x(0), y(0), z(0) {}
	Vec3(t _x, t _y, t _z) :x(_x), y(_y), z(_z) {}
	Vec3<t>(Matrix m);

	inline Vec3<t> operator +(const Vec3<t>& V) const { return Vec3<t>(x + V.x, y + V.y, z + V.z); }
	inline Vec3<t> operator -(const Vec3<t>& V) const { return Vec3<t>(x - V.x, y - V.y, z - V.z); }
	inline Vec3<t> operator*(float f) const { return Vec3<t>(f * x, f * y, f * z); }	
	inline Vec3<t> operator^(const Vec3<t>& V) const { return Vec3<t>(y * V.z - z * V.y, z * V.x - x * V.z, x * V.y - y * V.x); }	// Cross product
	inline t operator*(const Vec3<t>& V) const { return x * V.x + y * V.y + z * V.z; }	// Dot product
	inline t operator [](const int i) const { assert(i >= 0 && i < 3); return raw[i]; }
	template <class > friend std::ostream& operator <<(std::ostream& s, Vec3<t>& v);

	t magnitude2() { return (*this) * (*this); }
	t magnitude() { return std::sqrt(magnitude2()); }
	Vec3 & normalize() { *this = *this *(1./ magnitude()); return *this; }
};

typedef Vec2<int> Vec2i;
typedef Vec2<float> Vec2f;
typedef Vec3<int> Vec3i;
typedef Vec3<float> Vec3f;

const int DEFAULT_SIZE = 4;

class Matrix {
private:
	std::vector<std::vector<float>> m;
	int rows, cols;

public:
	Matrix(int r = DEFAULT_SIZE, int c = DEFAULT_SIZE);
	Matrix(Vec3f v);
	inline int nrows() { return rows; }
	inline int ncols() { return cols; }

	std::vector<float>& operator[](const int i);
	Matrix operator*(const Matrix& a);

	static Matrix identity(int dimension);
	Matrix inverse();
	Matrix transpose();

	friend std::ostream& operator<<(std::ostream& s, Matrix& m);
};

#endif