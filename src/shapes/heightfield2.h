
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "shape.h"
#include "shapes/trianglemesh.h"



class Heightfield2;


class Plane3D : public Shape{
public:
	Plane3D(Heightfield2* hf, const Point& p1, const Point& p2, const Point& p3, const Normal& n1, const Normal& n2, const Normal& n3);
	inline bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg, float minT, float maxT);
	inline bool IntersectP(const Ray &ray,float minT, float maxT);
	
	void GetShadingGeometry(const Transform &obj2world,
		const DifferentialGeometry &dg,
		DifferentialGeometry *dgShading) const;
	
	BBox ObjectBound() const { return BBox(); };
private:
	Point p1, p2, p3;
	Normal n1, n2, n3;
	Vector normal, e1, e2;
	float uvs[3][2];
	Heightfield2* hf;
	Vector dpdu, dpdv;
	Normal dndu, dndv;
	float lowestZ;
	float highestZ;
};

class HVoxel {
public:
	HVoxel(){
		planes[0] = NULL;
		planes[1] = NULL;
	}
	HVoxel(Heightfield2* hf, const Point& p1, const Point& p2, const Point& p3, const Point& p4, const Normal& n1, const Normal& n2, const Normal& n3, const Normal& n4);

	inline void freePlane(){
		for (int i = 0;i < 2;i++) {
			if (planes[i]) {
				delete planes[i];
				planes[i] = NULL;
			}
		}
	}
	inline bool IntersectP(const Ray &ray, float minT, float maxT);
	inline bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg, float minT, float maxT);
	
	Plane3D* planes[2];
};


// Heightfield2 Declarations
class Heightfield2 : public Shape {
public:
    // Heightfield2 Public Methods
    Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield2();
    bool CanIntersect() const;
    BBox ObjectBound() const;
	bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &ray) const;
	void GetShadingGeometry(const Transform &obj2world,
		const DifferentialGeometry &dg,
		DifferentialGeometry *dgShading) const {
		dg.shape->GetShadingGeometry(obj2world, dg, dgShading);
	}
private:
	inline int offset(int x, int y) const {
		return y*nVoxels[0] + x;
	}
    // Heightfield2 Private Data
    float *z;
    int nx, ny;
	int nv, np;
	int nVoxels[2];
	float width[2];
	float invWidth[2];
	BBox bound;
	Point *points;
	Normal* normals;
	HVoxel *voxels;
	float* uvs;
};


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD2_H
