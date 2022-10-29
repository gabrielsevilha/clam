/*
	Clam (C Linear Algebra Math) By Gabriel Sevilha.
	
	Library created for my own use, but if you like it, you can use it too, feel free to study the code.
	
	-------------------------------------------------------------------------------
	
	MIT License

	Copyright (c) 2020 Gabriel Sevilha

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
	
	-------------------------------------------------------------------------------
*/

#ifndef CLAM_LIBRARY_H
#define CLAM_LIBRARY_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//Types

#define INIT_MATRIX_3X3 {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0}
#define INIT_ZERO_MATRIX_3X3 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}

#define INIT_MATRIX_4X4 {1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0}
#define INIT_ZERO_MATRIX_4X4 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}

#define INIT_QUATERNION {0.0,0.0,0.0,1.0}
#define INIT_ZERO_QUATERNION {0.0,0.0,0.0,0.0}

typedef struct Vector2{

	float x, y;

}Vector2;

typedef struct Vector3{

	float x, y, z;

}Vector3;

typedef struct Vector4{

	float x, y, z, w;

}Vector4;

typedef float Matrix3x3[9];

typedef float Matrix4x4[16];

typedef float* Matrix3x3P;

typedef float* Matrix4x4P;

typedef struct Quaternion{

	float x, y, z, w;

}Quaternion;

//Vector2

Vector2 createVector2(float x, float y){

	return (Vector2){x, y};

}

Vector2 createCopyVector2(Vector2 v){

	return (Vector2){v.x, v.y};

}

Vector2 negateVector2(Vector2 v){

	return (Vector2){-v.x, -v.y};

}

Vector2 sumVector2(Vector2 v1, Vector2 v2){

	return (Vector2){v1.x+v2.x, v1.y+v2.y};

}

Vector2 subVector2(Vector2 v1, Vector2 v2){

	return (Vector2){v1.x-v2.x, v1.y-v2.y};

}

Vector2 mulVector2(Vector2 v1, Vector2 v2){

	return (Vector2){v1.x*v2.x, v1.y*v2.y};

}

Vector2 divVector2(Vector2 v1, Vector2 v2){

	return (Vector2){v1.x/v2.x, v1.y/v2.y};

}

Vector2 scaleVector2(Vector2 v, float s){
	
	return (Vector2){v.x*s, v.y*s};
	
}

Vector2 lerpVector2(Vector2 from, Vector2 dest, float t){

	Vector2 v = from;
	
	v.x += (dest.x - from.x) * t;
	v.y += (dest.y - from.y) * t;
	
	return v;

}

Vector2 normalizeVector2(Vector2 v){
	
	float length = sqrt((v.x*v.x)+(v.y*v.y));
	
	return (Vector2) {v.x / length, v.y / length};
	
}

float crossVector2(Vector2 v1, Vector2 v2){

	return (v1.x*v2.y - v1.y*v2.x);

}

Vector2 reflectVector2(Vector2 d, Vector2 n){
	
	float dot = d.x*n.x + d.y*n.y;
	
	return (Vector2){
		d.x - 2.0 * (dot * n.x),
		d.y - 2.0 * (dot * n.y)
	};
	
}

float dotVector2(Vector2 v1, Vector2 v2){

	return (float)(v1.x*v2.x + v1.y*v2.y);

}

float distanceVector2(Vector2 v1, Vector2 v2){

	return sqrt( pow(v1.x-v2.x,2.0f) + pow(v1.y-v2.y,2.0f) );

}

float lengthVector2(Vector2 v){
	
	return sqrt((v.x*v.x)+(v.y*v.y));
	
}

float angleVector2(Vector2 v1, Vector2 v2){
	
	float dot = v1.x*v2.x + v1.y*v2.y;
	float v1_lenght = sqrt(v1.x*v1.x + v1.y*v1.y);
	float v2_lenght = sqrt(v2.x*v2.x + v2.y*v2.y);
	
	return acosf(dot * (1.0f / (v1_lenght * v2_lenght))) * 180.0f/3.14159265358979323846f;
	
}

//Vector3

Vector3 createVector3(float x, float y, float z){

	return (Vector3){x, y ,z};

}

Vector3 createCopyVector3(Vector3 v){

	return (Vector3){v.x, v.y ,v.z};

}

Vector3 negateVector3(Vector3 v){

	return (Vector3){-v.x, -v.y, -v.z};

}

Vector3 sumVector3(Vector3 v1, Vector3 v2){

	return (Vector3){v1.x+v2.x, v1.y+v2.y, v1.z+v2.z};

}

Vector3 subVector3(Vector3 v1, Vector3 v2){

	return (Vector3){v1.x-v2.x, v1.y-v2.y, v1.z-v2.z};

}

Vector3 mulVector3(Vector3 v1, Vector3 v2){

	return (Vector3){v1.x*v2.x, v1.y*v2.y, v1.z*v2.z};

}

Vector3 divVector3(Vector3 v1, Vector3 v2){

	return (Vector3){v1.x/v2.x, v1.y/v2.y, v1.z/v2.z};

}

Vector3 scaleVector3(Vector3 v, float s){

	return (Vector3){v.x*s, v.y*s, v.z*s};

}

Vector3 lerpVector3(Vector3 from, Vector3 dest, float t){

	Vector3 v = from;
	
	v.x += (dest.x - from.x) * t;
	v.y += (dest.y - from.y) * t;
	v.z += (dest.z - from.z) * t;
	
	return v;

}

Vector3 normalizeVector3(Vector3 v){
	
	float length = sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
	
	return (Vector3) {v.x / length, v.y / length, v.z / length};
	
}

Vector3 crossVector3(Vector3 v1, Vector3 v2){

	return (Vector3){v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x};

}

Vector3 reflectVector3(Vector3 d, Vector3 n){
	
	float dot = d.x*n.x + d.y*n.y + d.z*n.z;
	
	return (Vector3){
		d.x - 2.0 * (dot * n.x),
		d.y - 2.0 * (dot * n.y),
		d.z - 2.0 * (dot * n.z)
	};
	
}

float dotVector3(Vector3 v1, Vector3 v2){

	return (float)(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);

}

float distanceVector3(Vector3 v1, Vector3 v2){

	return sqrt( pow(v1.x-v2.x,2.0f) + pow(v1.y-v2.y,2.0f) + pow(v1.z-v2.z,2.0f) );

}

float lengthVector3(Vector3 v){
	
	return sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
	
}

float angleVector3(Vector3 v1, Vector3 v2){
	
	float dot = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	float v1_lenght = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
	float v2_lenght = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
	
	return acosf(dot * (1.0f / (v1_lenght * v2_lenght))) * 180.0f/3.14159265358979323846f;
	
}

//Vector4

Vector4 createVector4(float x, float y, float z, float w){

	return (Vector4){x, y, z, w};

}

Vector4 createCopyVector4(Vector4 v){

	return (Vector4){v.x, v.y, v.z, v.w};

}

Vector4 negateVector4(Vector4 v){

	return (Vector4){-v.x, -v.y, -v.z, -v.w};

}

Vector4 sumVector4(Vector4 v1, Vector4 v2){

	return (Vector4){v1.x+v2.x, v1.y+v2.y, v1.z+v2.z, v1.w+v2.w};

}

Vector4 subVector4(Vector4 v1, Vector4 v2){

	return (Vector4){v1.x-v2.x, v1.y-v2.y, v1.z-v2.z, v1.w-v2.w};

}

Vector4 mulVector4(Vector4 v1, Vector4 v2){

	return (Vector4){v1.x*v2.x, v1.y*v2.y, v1.z*v2.z, v1.w*v2.w};

}

Vector4 divVector4(Vector4 v1, Vector4 v2){

	return (Vector4){v1.x/v2.x, v1.y/v2.y, v1.z/v2.z, v1.w/v2.w};

}

Vector4 scaleVector4(Vector4 v, float s){

	return (Vector4){v.x*s, v.y*s, v.z*s, v.w*s};

}

Vector4 lerpVector4(Vector4 from, Vector4 dest, float t){

	Vector4 v = from;
	
	v.x += (dest.x - from.x) * t;
	v.y += (dest.y - from.y) * t;
	v.z += (dest.z - from.z) * t;
	v.w += (dest.w - from.w) * t;
	
	return v;

}

Vector4 normalizeVector4(Vector4 v){
	
	float length = sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z)+(v.w*v.w));
	
	return (Vector4) {v.x / length, v.y / length, v.z / length, v.w / length};
	
}

Vector4 reflectVector4(Vector4 d, Vector4 n){
	
	float dot = d.x*n.x + d.y*n.y + d.z*n.z + d.w*d.w;
	
	return (Vector4){
		d.x - 2.0 * (dot * n.x),
		d.y - 2.0 * (dot * n.y),
		d.z - 2.0 * (dot * n.z),
		d.w - 2.0 * (dot * n.w)
	};
	
}

float dotVector4(Vector4 v1, Vector4 v2){

	return (float)(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w);

}

float distanceVector4(Vector4 v1, Vector4 v2){

	return sqrt( pow(v1.x-v2.x,2.0f) + pow(v1.y-v2.y,2.0f) + pow(v1.z-v2.z,2.0f) + pow(v1.w-v2.w,2.0f) );

}

float lengthVector4(Vector4 v){
	
	return sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z)+(v.w*v.w));
	
}

float angleVector4(Vector4 v1, Vector4 v2){
	
	float dot = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
	float v1_lenght = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z + v1.w*v1.w);
	float v2_lenght = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z + v2.w*v2.w);
	
	return acosf(dot * (1.0f / (v1_lenght * v2_lenght))) * 180.0f/3.14159265358979323846f;
	
}

//Matrix3x3

Matrix3x3P createMatrix3x3(float major_value){

	Matrix3x3P m = (Matrix3x3P) calloc(9, sizeof(float)*9);
	
	m[0] = major_value;
	m[4] = major_value;
	m[8] = major_value;
	
	return m;

}

Matrix3x3P createCopyMatrix3x3(Matrix3x3P m){
	
	Matrix3x3P r = (Matrix3x3P) malloc(sizeof(float)*9);
	
	r[0] = m[0], r[1] = m[1], r[2] = m[2],
	r[3] = m[3], r[4] = m[4], r[5] = m[5],
	r[6] = m[6], r[7] = m[7], r[8] = m[8];
	
	return r;

}

void copyMatrix3x3(Matrix3x3P src, Matrix3x3P dest){
	
	dest[0] = src[0], dest[1] = src[1], dest[2] = src[2],
	dest[3] = src[3], dest[4] = src[4], dest[5] = src[5],
	dest[6] = src[6], dest[7] = src[7],	dest[8] = src[8];
	
}

void identityMatrix3x3(Matrix3x3P m){

	m[0] = 1.0f, m[1] = 0.0f, m[2] = 0.0f,
	m[3] = 0.0f, m[4] = 1.0f, m[5] = 0.0f,
	m[6] = 0.0f, m[7] = 0.0f, m[8] = 1.0f;
	
}

void multiplyMatrix3x3(Matrix3x3P m1, Matrix3x3P m2, Matrix3x3P dest){
	
	float a0 = m1[0], a1 = m1[1], a2 = m1[2],
	a3 = m1[3],	a4 = m1[4], a5 = m1[5],
	a6 = m1[6], a7 = m1[7],	a8 = m1[8],

	b0 = m2[0], b1 = m2[1], b2 = m2[2],
	b3 = m2[3],	b4 = m2[4], b5 = m2[5],
	b6 = m2[6], b7 = m2[7],	b8 = m2[8];

	dest[0] = a0 * b0 + a3 * b1 + a6 * b2;
	dest[1] = a1 * b0 + a4 * b1 + a7 * b2;
	dest[2] = a2 * b0 + a5 * b1 + a8 * b2;
	
	dest[3] = a0 * b3 + a3 * b4 + a6 * b5;
	dest[4] = a1 * b3 + a4 * b4 + a7 * b5;
	dest[5] = a2 * b3 + a5 * b4 + a8 * b5;
	
	dest[6] = a0 * b6 + a3 * b7 + a6 * b8;
	dest[7] = a1 * b6 + a4 * b7 + a7 * b8;
	dest[8] = a2 * b6 + a5 * b7 + a8 * b8;
		
}

Vector2 multiplyMatrix3x3Vector2(Matrix3x3P m, Vector2 v, float z){

	return (Vector2){
		m[0] * v.x + m[3] * v.y + m[6] * z,
		m[1] * v.x + m[4] * v.y + m[7] * z
	};
	
}

Vector3 multiplyMatrix3x3Vector3(Matrix3x3P m, Vector3 v){

	return (Vector3){
		m[0] * v.x + m[3] * v.y + m[6] * v.z,
		m[1] * v.x + m[4] * v.y + m[7] * v.z,
		m[2] * v.x + m[5] * v.y + m[8] * v.z
	};
	
}

void transposeMatrix3x3(Matrix3x3P m){
	
	float m1 = m[1], m2 = m[2],
	m3 = m[3], m5 = m[5],
	m6 = m[6], m7 = m[7];
		
	m[1] = m3, m[2] = m6,
	m[3] = m1, m[5] = m7,
	m[6] = m2, m[7] = m5;
	
}

float determinantMatrix3x3(Matrix3x3P m){
	
	float m0 = m[0], m1 = m[1], m2 = m[2],
		m3 = m[3], m4 = m[4], m5 = m[5],
		m6 = m[6], m7 = m[7], m8 = m[8];
	
	return (m0*m4*m8 + m1*m5*m6 + m2*m3*m7) - (m2*m4*m6 + m0*m5*m7 + m1*m3*m8);
	
}

void inverseMatrix3x3(Matrix3x3P m){
	
	float m0 = m[0], m1 = m[1], m2 = m[2],
		m3 = m[3], m4 = m[4], m5 = m[5],
		m6 = m[6], m7 = m[7], m8 = m[8];
	
	float r[9];
	
	r[0] = m4*m8 - m5*m7;
	r[3] = -(m3*m8 - m5*m6);
	r[6] = m3*m7 - m4*m6;
	
	r[1] = -(m1*m8 - m2*m7);
	r[4] = m0*m8 - m2*m6;
	r[7] = -(m0*m7 - m1*m6);
	
	r[2] = m1*m5 - m2*m4;
	r[5] = -(m0*m5 - m2*m3);
	r[8] = m0*m4 - m1*m3;
	
	float d = 1.0f / ( (m0*m4*m8 + m1*m5*m6 + m2*m3*m7) - (m2*m4*m6 + m0*m5*m7 + m1*m3*m8) );
	
	m[0] = r[0] * d, m[1] = r[1] * d, m[2] = r[2] * d,
	m[3] = r[3] * d, m[4] = r[4] * d, m[5] = r[5] * d,
	m[6] = r[6] * d, m[7] = r[7] * d, m[8] = r[8] * d;
	
}

void translateMatrix3x3(Matrix3x3P m, Vector2 v){
	
	float r[] = {
		1.0,0.0,0.0,
		0.0,1.0,0.0,
		v.x,v.y,1.0
	};
	
	multiplyMatrix3x3(m,r,m);
	
}

void scaleMatrix3x3(Matrix3x3P m, Vector2 v){

	float r[] = {
		v.x,0.0,0.0,
		0.0,v.y,0.0,
		0.0,0.0,1.0
	};
	
	multiplyMatrix3x3(m,r,m);
	
}

void rotateMatrix3x3(Matrix3x3P m, float angle){
	
	float c = cosf(angle);
	float s = sinf(angle);
	
	float r[] = {
		c  , s  ,0.0,
		-s , c  ,0.0,
		0.0,0.0,1.0
	};
	
	multiplyMatrix3x3(m,r,m);

}

//Matrix4x4

Matrix4x4P createMatrix4x4(float major_value){

	Matrix4x4P m = (Matrix4x4P) calloc(16, sizeof(float)*16);
	
	m[0] = major_value;
	m[5] = major_value;
	m[10] = major_value;
	m[15] = major_value;
	
	return m;

}

Matrix4x4P createCopyMatrix4x4(Matrix4x4P m){
	
	Matrix4x4P r = (Matrix4x4P) malloc(sizeof(float)*16);
	
	r[0] = m[0], r[1] = m[1], r[2] = m[2], r[3] = m[3],
	r[4] = m[4], r[5] = m[5], r[6] = m[6], r[7] = m[7],
	r[8] = m[8], r[9] = m[9], r[10] = m[10], r[11] = m[11],
	r[12] = m[12], r[13] = m[13], r[14] = m[14], r[15] = m[15];
	
	return r;

}

void copyMatrix4x4(Matrix4x4P src, Matrix4x4P dest){
	
	dest[0] = src[0], dest[1] = src[1], dest[2] = src[2], dest[3] = src[3],
	dest[4] = src[4], dest[5] = src[5], dest[6] = src[6], dest[7] = src[7],
	dest[8] = src[8], dest[9] = src[9], dest[10] = src[10], dest[11] = src[11],
	dest[12] = src[12], dest[13] = src[13], dest[14] = src[14], dest[15] = src[15];
	
}

void identityMatrix4x4(Matrix4x4P m){

	m[0] = 1.0f, m[1] = 0.0f, m[2] = 0.0f, m[3] = 0.0f,
	m[4] = 0.0f, m[5] = 1.0f, m[6] = 0.0f, m[7] = 0.0f,
	m[8] = 0.0f, m[9] = 0.0f, m[10] = 1.0f, m[11] = 0.0f,
	m[12] = 0.0f, m[13] = 0.0f, m[14] = 0.0f, m[15] = 1.0f;
	
}

void multiplyMatrix4x4(Matrix4x4P m1, Matrix4x4P m2, Matrix4x4P dest){
	
	float a0 = m1[0], a1 = m1[1], a2 = m1[2], a3 = m1[3],
	a4 = m1[4], a5 = m1[5], a6 = m1[6], a7 = m1[7],
	a8 = m1[8], a9 = m1[9], a10 = m1[10], a11 = m1[11],
	a12 = m1[12], a13 = m1[13], a14 = m1[14], a15 = m1[15],

	b0 = m2[0], b1 = m2[1], b2 = m2[2], b3 = m2[3],
	b4 = m2[4], b5 = m2[5], b6 = m2[6], b7 = m2[7],
	b8 = m2[8], b9 = m2[9], b10 = m2[10], b11 = m2[11],
	b12 = m2[12], b13 = m2[13], b14 = m2[14], b15 = m2[15];

	dest[0] = a0 * b0 + a4 * b1 + a8 * b2 + a12 * b3;
	dest[1] = a1 * b0 + a5 * b1 + a9 * b2 + a13 * b3;
	dest[2] = a2 * b0 + a6 * b1 + a10 * b2 + a14 * b3;
	dest[3] = a3 * b0 + a7 * b1 + a11 * b2 + a15 * b3;

	dest[4] = a0 * b4 + a4 * b5 + a8 * b6 + a12 * b7;
	dest[5] = a1 * b4 + a5 * b5 + a9 * b6 + a13 * b7;
	dest[6] = a2 * b4 + a6 * b5 + a10 * b6 + a14 * b7;
	dest[7] = a3 * b4 + a7 * b5 + a11 * b6 + a15 * b7;

	dest[8] = a0 * b8 + a4 * b9 + a8 * b10 + a12 * b11;
	dest[9] = a1 * b8 + a5 * b9 + a9 * b10 + a13 * b11;
	dest[10] = a2 * b8 + a6 * b9 + a10 * b10 + a14 * b11;
	dest[11] = a3 * b8 + a7 * b9 + a11 * b10 + a15 * b11;

	dest[12] = a0 * b12 + a4 * b13 + a8 * b14 + a12 * b15;
	dest[13] = a1 * b12 + a5 * b13 + a9 * b14 + a13 * b15;
	dest[14] = a2 * b12 + a6 * b13 + a10 * b14 + a14 * b15;
	dest[15] = a3 * b12 + a7 * b13 + a11 * b14 + a15 * b15;
		
}

Vector3 multiplyMatrix4x4Vector3(Matrix4x4P m, Vector3 v, float w){

	return (Vector3){
		m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * w,
		m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * w,
		m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * w
	};
	
}

Vector4 multiplyMatrix4x4Vector4(Matrix4x4P m, Vector4 v){

	return (Vector4){
		m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * v.w,
		m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * v.w,
		m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w,
		m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w
	};
	
}

void transposeMatrix4x4(Matrix4x4P m){
	
	float  m1 = m[1], m2 = m[2], m3 = m[3],
		m4 = m[4], m6 = m[6], m7 = m[7],
		m8 = m[8], m9 = m[9], m11 = m[11],
		m12 = m[12], m13 = m[13], m14 = m[14];
		
	m[1] = m4, m[2] = m8, m[3] = m12,
	m[4] = m1, m[6] = m9, m[7] = m13,
	m[8] = m2, m[9] = m6, m[11] = m14,
	m[12] = m3, m[13] = m7, m[14] = m11;
	
}

float determinantMatrix4x4(Matrix4x4P m){
	
	float m0 = m[0], m1 = m[1], m2 = m[2], m3 = m[3],
		m4 = m[4], m5 = m[5], m6 = m[6], m7 = m[7],
		m8 = m[8], m9 = m[9], m10 = m[10], m11 = m[11],
		m12 = m[12], m13 = m[13], m14 = m[14], m15 = m[15];
	
	float a0 = m10*m15 - m11*m14,
		a1 = m9*m15 - m11*m13,
		a2 = m9*m14 - m10*m13,
		a3 = m8*m15 - m11*m12,
		a4 = m8*m14 - m10*m12,
		a5 = m8*m13 - m9*m12;
		
	
	float r1 = (m5 * a0) - (m6 * a1) + (m7 * a2),
		r2 = (m4 * a0) - (m6 * a3) + (m7 * a4),
		r3 = (m4 * a1) - (m5 * a3) + (m7 * a5),
		r4 = (m4 * a2) - (m5 * a4) + (m6 * a5);
	
	return m0*r1 - m1*r2 + m2*r3 - m3*r4;
	
}

void inverseMatrix4x4(Matrix4x4P m){

	float r[16];
	
	float   m0 = m[0], m1 = m[1], m2 = m[2], m3 = m[3],
			m4 = m[4], m5 = m[5], m6 = m[6], m7 = m[7],
			m8 = m[8], m9 = m[9], m10 = m[10], m11 = m[11],
			m12 = m[12], m13 = m[13], m14 = m[14], m15 = m[15];
	
	r[0] = (m5*m10*m15 + m6*m11*m13 + m7*m9*m14) - (m7*m10*m13 + m5*m11*m14 + m6*m9*m15);
	r[1] = - ( (m1*m10*m15 + m2*m11*m13 + m3*m9*m14) - (m3*m10*m13 + m1*m11*m14 + m2*m9*m15) );
	r[2] = (m1*m6*m15 + m2*m7*m13 + m3*m5*m14) - (m3*m6*m13 + m1*m7*m14 + m2*m5*m15);
	r[3] = - ( (m1*m6*m11 + m2*m7*m9 + m3*m5*m10) - (m3*m6*m9 + m1*m7*m10 + m2*m5*m11) );
	
	r[4] = - ( (m4*m10*m15 + m6*m11*m12 + m7*m8*m14) - (m7*m10*m12 + m4*m11*m14 + m6*m8*m15) );
	r[5] = (m0*m10*m15 + m2*m11*m12 + m3*m8*m14) - (m3*m10*m12 + m0*m11*m14 + m2*m8*m15);
	r[6] = - ( (m0*m6*m15 + m2*m7*m12 + m3*m4*m14) - (m3*m6*m12 + m0*m7*m14 + m2*m4*m15) );
	r[7] = (m0*m6*m11 + m2*m7*m8 + m3*m4*m10) - (m3*m6*m8 + m0*m7*m10 + m2*m4*m11);
	
	r[8] = (m4*m9*m15 + m5*m11*m12 + m7*m8*m13) - (m7*m9*m12 + m4*m11*m13 + m5*m8*m15);
	r[9] = - ( (m0*m9*m15 + m1*m11*m12 + m3*m8*m13) - (m3*m9*m12 + m0*m11*m13 + m1*m8*m15) );
	r[10] = (m0*m5*m15 + m1*m7*m12 + m3*m4*m13) - (m3*m5*m12 + m0*m7*m13 + m1*m4*m15);
	r[11] = - ((m0*m5*m11 + m1*m7*m8 + m3*m4*m9) - (m3*m5*m8 + m0*m7*m9 + m1*m4*m11) );
	
	r[12] = - ( (m4*m9*m14 + m5*m10*m12 + m6*m8*m13) - (m6*m9*m12 + m4*m10*m13 + m5*m8*m14) );
	r[13] = (m0*m9*m14 + m1*m10*m12 + m2*m8*m13) - (m2*m9*m12 + m0*m10*m13 + m1*m8*m14);
	r[14] = - ( (m0*m5*m14 + m1*m6*m12 + m2*m4*m13) - (m2*m5*m12 + m0*m6*m13 + m1*m4*m14) );
	r[15] = (m0*m5*m10 + m1*m6*m8 + m2*m4*m9) - (m2*m5*m8 + m0*m6*m9 + m1*m4*m10);
	
	float d = 1.0 / (m0*r[0] + m1*r[4] + m2*r[8] + m3*r[12]);
	
	m[0] = r[0] * d, m[1] = r[1] * d, m[2] = r[2] * d, m[3] = r[3] * d,
	m[4] = r[4] * d, m[5] = r[5] * d, m[6] = r[6] * d, m[7] = r[7] * d,
	m[8] = r[8] * d, m[9] = r[9] * d, m[10] = r[10] * d, m[11] = r[11] * d,
	m[12] = r[12] * d, m[13] = r[13] * d, m[14] = r[14] * d, m[15] = r[15] * d;
	
}

void translateMatrix4x4(Matrix4x4P m, Vector3 v){
	
	float r[] = {
		1.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,
		0.0,0.0,1.0,0.0,
		v.x,v.y,v.z,1.0
	};
	
	multiplyMatrix4x4(m,r,m);
	
}

void scaleMatrix4x4(Matrix4x4P m, Vector3 v){

	float r[] = {
		v.x,0.0,0.0,0.0,
		0.0,v.y,0.0,0.0,
		0.0,0.0,v.z,0.0,
		0.0,0.0,0.0,1.0
	};
	
	multiplyMatrix4x4(m,r,m);
	
}

void rotateXMatrix4x4(Matrix4x4P m, float angle){

	float c = cosf(angle);
	float s = sinf(angle);
	
	float r[] = {
		1.0,0.0,0.0,0.0,
		0.0,c,s,0.0,
		0.0,-s,c,0.0,
		0.0,0.0,0.0,1.0
	};
	
	multiplyMatrix4x4(m,r,m);

}

void rotateYMatrix4x4(Matrix4x4P m, float angle){

	float c = cosf(angle);
	float s = sinf(angle);
	
	float r[] = {
		c,0.0,-s,0.0,
		0.0,1.0,0.0,0.0,
		s,0.0,c,0.0,
		0.0,0.0,0.0,1.0
	};
	
	multiplyMatrix4x4(m,r,m);

}

void rotateZMatrix4x4(Matrix4x4P m, float angle){
	
	float c = cosf(angle);
	float s = sinf(angle);
	
	float r[] = {
		c,s,0.0,0.0,
		-s,c,0.0,0.0,
		0.0,0.0,1.0,0.0,
		0.0,0.0,0.0,1.0
	};
	
	multiplyMatrix4x4(m,r,m);

}

void rotateMatrix4x4(Matrix4x4P m, float angle, Vector3 v){

	float c = cosf(angle);
	float s = sinf(angle);
	float t = (1.0f - c);
	
	float r[] = {
		1.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,
		0.0,0.0,1.0,0.0,
		0.0,0.0,0.0,1.0
	};

	v = normalizeVector3(v);
	
	r[0] = c + (v.x*v.x) * t;
	r[1] = t * v.x * v.y + s * v.z;
	r[2] = t * v.x * v.z - s * v.y;
	
	r[4] = t * v.x * v.y - s * v.z;
	r[5] = t * (v.y*v.y) + c;
	r[6] = t * v.y * v.z + s * v.x;
	
	r[8] = t * v.x * v.z + s * v.y;
	r[9] = t * v.y * v.z - s * v.x;
	r[10] = t * (v.z*v.z) + c;
	
	multiplyMatrix4x4(m,r,m);

}

//Quaternions

Quaternion createQuaternion(float x, float y, float z, float w){
	
	return (Quaternion) {x, y, z, w};

}

Quaternion createCopyQuaternion(Quaternion q){

	return (Quaternion) {q.x, q.y, q.z, q.w};
	
}

Quaternion initQuaternionAxis(float angle, Vector3 axis){
	
	Vector3 n = normalizeVector3(axis);
	
	float a = angle * 0.5f;
	float s = sinf(a);
	float c = cosf(a);
	
	return (Quaternion){ s * n.x, s * n.y, s * n.z, c };
	
}

float dotQuaternion(Quaternion q1, Quaternion q2){
	
	Vector4 v0 = (Vector4){q1.x, q1.y, q1.z, q1.w};
	Vector4 v1 = (Vector4){q2.x, q2.y, q2.z, q2.w};
	
	return dotVector4(v0,v1);
	
}

Quaternion lerpQuaternion(Quaternion from, Quaternion dest, float t){
	
	Vector4 v0 = (Vector4){from.x, from.y, from.z, from.w};
	Vector4 v1 = (Vector4){dest.x, dest.y, dest.z, dest.w};
	
	Vector4 r = lerpVector4(v0, v1, t);
	
	return (Quaternion){r.x, r.y, r.z, r.w};
	
}

Quaternion slerpQuaternion(Quaternion from, Quaternion dest, float t){

	Vector4 v0 = (Vector4){from.x, from.y, from.z, from.w};
	Vector4 v1 = (Vector4){dest.x, dest.y, dest.z, dest.w};
	
	float cos_theta = dotQuaternion(from, dest);
	
	if(fabsf(cos_theta) >= 1.0f){
		return from;
	}
	
	if(cos_theta < 0.0){
		negateVector4(v0);
		cos_theta = -cos_theta;
	}
	
	float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
	
	if(fabsf(sin_theta) < 0.0001f){
		return lerpQuaternion(from,dest,t);
	}
	
	float angle = acosf(cos_theta);
	
	v0 = scaleVector4(v0, sinf((1.0f - t) * angle));
	v1 = scaleVector4(v1, sinf(t * angle));
	v0 = sumVector4(v0, v1);
	v0 = scaleVector4(v0, 1.0f/sin_theta);
	
	return (Quaternion){v0.x, v0.y, v0.z, v0.w};
	
}

void quaternionToMatrix4x4(Quaternion q, Matrix4x4 dest){
	
	float x = q.x, y = q.y, z = q.z, w = q.w;
	
	dest[3] = 0.0f,
	dest[7] = 0.0f,
	dest[11] = 0.0f,
	dest[12] = 0.0f,
	dest[13] = 0.0f,
	dest[14] = 0.0f,
	dest[15] = 1.0f;
	
}

//Camera

void orthographicMatrix(float left, float right, float bottom, float top, float near, float far, Matrix4x4P m){

	float dif_right_left = right - left;
	float dif_top_bottom = top - bottom;
	float dif_far_near = far - near;

	m[0] = 2.0f / (dif_right_left);
	m[5] = 2.0f / (dif_top_bottom);
	m[10] = -2.0f / (dif_far_near);
	m[12] = -( (right+left)/(dif_right_left) );
	m[13] = -( (top+bottom)/(dif_top_bottom) );
	m[14] = -( (far+near)/(dif_far_near) );
	m[15] = 1.0f;

}

void perspectiveMatrix(float fov, float aspect, float near, float far, Matrix4x4P m){

	fov *= 3.14159265358979323846f/180.0f;

	float tan_half_fov = tanf( fov / 2.0f );
	float dif_far_near = far - near;
	float scale_x = 1.0f / (aspect * tan_half_fov );
	float scale_y = 1.0f / (tan_half_fov);
	
	m[0] = scale_x;
	m[5] = scale_y;
	m[10] = -( (far + near) / dif_far_near );
	m[11] = -1.0f;
	m[14] = -( (2.0f * near * far) / dif_far_near);
	m[15] = 0.0f;

}

void lookatMatrix(Vector3 eye, Vector3 center, Vector3 up, Matrix4x4P m){

	Vector3 z_axis = subVector3(center,eye);
	z_axis = normalizeVector3(z_axis);
	
	Vector3 x_axis = crossVector3(z_axis,up);
	x_axis = normalizeVector3(x_axis);
	
	Vector3 y_axis = crossVector3(x_axis,z_axis);
	
	m[0] = x_axis.x, m[1] = y_axis.x, m[2] = -z_axis.x, m[3] = 0.0f,
	m[4] = x_axis.y, m[5] = y_axis.y, m[6] = -z_axis.y, m[7] = 0.0f;
	m[8] = x_axis.z, m[9] = y_axis.z, m[10] = -z_axis.z, m[11] = 0.0f;	
	m[12] = -dotVector3(x_axis, eye),
	m[13] = -dotVector3(y_axis, eye),
	m[14] = -dotVector3((Vector3){-z_axis.x,-z_axis.y,-z_axis.z}, eye),
	m[15] = 1.0f;
	
}

//Utils

float degToRad(float angle){

	return angle * 3.14159265358979323846f / 180.0f;

}

float radToDeg(float angle){

	return angle * 180.0f / 3.14159265358979323846f;

}

Vector2 arrayToVector2(float v[2]){
	
	return (Vector2){v[0], v[1]};	
	
}

Vector3 arrayToVector3(float v[3]){
	
	return (Vector3){v[0], v[1], v[2]};	
	
}

Vector4 arrayToVector4(float v[4]){
	
	return (Vector4){v[0], v[1], v[2], v[3]};	
	
}

int compVector2(Vector2 v1, Vector2 v2){
	
	return (v1.x == v2.x && v1.y == v2.y);
	
}

int compVector3(Vector3 v1, Vector3 v2){
	
	return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
	
}

int compVector4(Vector4 v1, Vector4 v2){
	
	return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z && v1.w == v2.w);
	
}

int compMatrix3x3(Matrix3x3P m1, Matrix3x3P m2){

	return (m1[0] == m2[0] && m1[1] == m2[1] && m1[2] == m2[2]
	&& m1[3] == m2[3] && m1[4] == m2[4] && m1[5] == m2[5]
	&& m1[6] == m2[6] && m1[7] == m2[7]	&& m1[8] == m2[8]);

}

int compMatrix4x4(Matrix4x4P m1, Matrix4x4P m2){

	return (m1[0] == m2[0] && m1[1] == m2[1] && m1[2] == m2[2] && m1[3] == m2[3]
	&& m1[4] == m2[4] && m1[5] == m2[5] && m1[6] == m2[6] && m1[7] == m2[7]
	&& m1[8] == m2[8] && m1[9] == m2[9] && m1[10] == m2[10] && m1[11] == m2[11]
	&& m1[12] == m2[12] && m1[13] == m2[13] && m1[14] == m2[14] && m1[15] == m2[15]);

}

void printVector2(Vector2 v){

	printf("%f %f\n",v.x,v.y);
	
}

void printVector3(Vector3 v){

	printf("%f %f %f\n",v.x,v.y,v.z);
	
}

void printVector4(Vector4 v){

	printf("%f %f %f %f\n",v.x,v.y,v.z,v.w);
	
}

void printMatrix3x3(Matrix3x3P m){

	printf("%f %f %f\n%f %f %f\n%f %f %f\n",
		m[0],m[1],m[2],
		m[3],m[4],m[5],
		m[6],m[7],m[8]
	);

}

void printMatrix4x4(Matrix4x4P m){

	printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n",
		m[0],m[1],m[2],m[3],
		m[4],m[5],m[6],m[7],
		m[8],m[9],m[10],m[11],
		m[12],m[13],m[14],m[15]
	);

}

void printQuaternion(Quaternion q){
	
	printf("%f %f %f %f\n",q.x,q.y,q.z,q.w);
	
}

//Extra

int rayPlane(Vector3 origin, Vector3 direction, Vector3 plane_position, Vector3 plane_normal, float* distance){
	
	float d = dotVector3(direction, plane_normal);
	float p1 = dotVector3(plane_position, negateVector3(plane_normal));
	float p2 = dotVector3(origin, plane_normal);
	
	*distance = -(p1 + p2) / d;
	
	if(*distance >= 0.0f) return 1;
	
	return 0;
	
}

int rayTriangle(Vector3 origin, Vector3 direction, Vector3 v0, Vector3 v1, Vector3 v2, float* d, float* u, float* v){

	float epsilon = 0.0000001;
	
	Vector3 edge1 = subVector3(v1,v0);
	Vector3 edge2 = subVector3(v2,v0);
	
	Vector3 h = crossVector3(direction,edge2);
	
	float a = dotVector3(edge1,h);
	
	if(a > -epsilon && a < epsilon)
		return 0;
	
	float f = 1.0 / a;
	
	Vector3 s = subVector3(origin,v0);
	*u = f * dotVector3(s,h);

	if(*u < 0.0f || *u > 1.0f)
		return 0;
	
	Vector3 q = crossVector3(s,edge1);
	*v = f * dotVector3(direction,q);
	
	if(*v < 0.0f || *u + *v > 1.0f)
		return 0;
	
	*d = f * dotVector3(edge2,q);
	
	if(*d > epsilon)
		return 1;
	else
		return 0;

}

Vector3 rgbToHsv(float r, float g, float b){

    r /= 255.0f, g /= 255.0f, b /= 255.0f;
    
    float cmax = fmax(r, fmax(g, b));
    float cmin = fmin(r, fmin(g, b));
    float diff = cmax - cmin;
    float h = -1, s = -1;
    
    if (cmax == cmin)
        h = 0;
    else if (cmax == r) 
        h = (int)(60 * ((g - b) / diff) + 360) % 360;
    else if (cmax == g) 
        h = (int)(60 * ((b - r) / diff) + 120) % 360;
    else if (cmax == b) 
        h = (int)(60 * ((r - g) / diff) + 240) % 360;

    if (cmax == 0) 
        s = 0; 
    else
        s = (diff / cmax) * 100; 
        
    float v = cmax * 100; 
    
    return (Vector3) {h, s, v};
    
}

Vector3 hsvToRgb(float h, float s,float v){

	if(h > 360 || h < 0 || s > 100 || s < 0 || v > 100 || v < 0)
		return (Vector3){0,0,0};

	float _s = s / 100.0f, _v = v / 100.0f;
	float c = _s * _v;
	float x = c * ( 1 - fabs( fmod( h/60.0, 2) - 1 ) );

	float r, g, b;

	if(h >= 0 && h < 60)
		r = c, g = x, b = 0;
	else if(h >= 60 && h < 120)
		r = x, g = c, b = 0;
	else if(h >= 120 && h < 180)
		r = 0, g = c, b = x;
	else if(h >= 180 && h < 240)
		r = 0, g = x, b = c;
	else if(h >= 240 && h < 300)
		r = x, g = 0, b = c;
	else
		r = c, g = 0, b = x;
		
	float m = _v - c;
		
	return (Vector3) { (r+m)*255, (g+m)*255, (b+m)*255 };
    
}

#endif
