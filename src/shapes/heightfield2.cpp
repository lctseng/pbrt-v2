
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


 // shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

Plane3D::Plane3D(Heightfield2* vhf, const Point& pp1, const Point& pp2, const Point& pp3) {
	hf = vhf;
	for (int i = 0;i < 2;i++) {
		uvs[0][i] = pp1[i];
	}
	for (int i = 0;i < 2;i++) {
		uvs[1][i] = pp2[i];
	}
	for (int i = 0;i < 2;i++) {
		uvs[2][i] = pp3[i];
	}
	p1 = (*hf->ObjectToWorld)(pp1);
	p2 = (*hf->ObjectToWorld)(pp2);
	p3 = (*hf->ObjectToWorld)(pp3);
	// compute normal
	e1 = p2 - p1;
	e2 = p3 - p1;
	normal = Cross(e1, e2);
	// pre-process
	// Compute triangle partial derivatives
	float du1 = uvs[0][0] - uvs[2][0];
	float du2 = uvs[1][0] - uvs[2][0];
	float dv1 = uvs[0][1] - uvs[2][1];
	float dv2 = uvs[1][1] - uvs[2][1];
	// Compute deltas for triangle partial derivatives
	Vector dp1, dp2;
	dp1 = p1 - p3;
	dp2 = p2 - p3;
	float determinant = du1 * dv2 - dv1 * du2;
	if (determinant == 0.f) {
		// Handle zero determinant for triangle partial derivative matrix
		CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
	}
	else {
		float invdet = 1.f / determinant;
		dpdu = (dv2 * dp1 - dv1 * dp2) * invdet;
		dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
	}
}
bool Plane3D::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg, float minT, float maxT) {
	float t;

	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector s = ray.o - p1;
	float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	Vector s2 = Cross(s, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	t = Dot(e2, s2) * invDivisor;
	if (t < minT || t > maxT)
		return false;
	

	// Interpolate $(u,v)$ triangle parametric coordinates
	float b0 = 1 - b1 - b2;
	float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
	float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];

	// Fill in _DifferentialGeometry_ from triangle hit
	*dg = DifferentialGeometry(ray(t), dpdu, dpdv,
		Normal(0, 0, 0), Normal(0, 0, 0),
		tu, tv, hf);
	*tHit = t;
	*rayEpsilon = 1e-3f * *tHit;
	
	
	return true;
}
bool Plane3D::IntersectP(const Ray &ray,float minT, float maxT) {
	float t;

	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector s = ray.o - p1;
	float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	Vector s2 = Cross(s, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	t = Dot(e2, s2) * invDivisor;
	if (t < minT || t > maxT)
		return false;

	return true;
}


HVoxel::HVoxel(Heightfield2 *hf, const Point& p1, const Point& p2, const Point& p3, const Point& p4)
{
	planes[0] = Plane3D(hf, p1, p2, p4);
	planes[1] = Plane3D(hf, p1, p4, p3);
}
bool HVoxel::IntersectP(const Ray &ray, float minT, float maxT) {
	for (auto& plane : planes) {
		if (plane.IntersectP(ray,minT, maxT)) {
			return true;
		}
	}
	return false;
}
bool HVoxel::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg, float minT, float maxT) {
	for (auto& plane : planes) {
		if (plane.Intersect(ray, tHit, rayEpsilon, dg, minT, maxT)) {
			return true;
		}
	}
	return false;
}

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
	bool ro, int x, int y, const float *zs)
	: Shape(o2w, w2o, ro) {
	nx = x;
	ny = y;
	z = new float[nx*ny];
	memcpy(z, zs, nx*ny * sizeof(float));
	// additional: accerlation structures
	bound = WorldBound();
	// points
	invWidth[0] = (float)(nx - 1);
	invWidth[1] = (float)(ny - 1);
	width[0] = 1.0f / invWidth[0];
	width[1] = 1.0f / invWidth[1];
	np = x * y;
	points = new Point[np];//AllocAligned<Point>(np);
	uvs = new float[np * 2];
	int pos = 0;
	for (int vy = 0; vy < ny; ++vy) {
		for (int vx = 0; vx < nx; ++vx) {
			points[pos].x = uvs[2 * pos] = (float)vx * width[0];
			points[pos].y = uvs[2 * pos + 1] = (float)vy * width[1];
			points[pos].z = z[pos];
			pos++;
		}
	}
	// voxels
	nVoxels[0] = nx - 1;
	nVoxels[1] = ny - 1;
	nv = nVoxels[0] * nVoxels[1];
	voxels = new HVoxel[nv];//AllocAligned<HVoxel>(nv);
	for (int j = 0;j < nVoxels[1];j++) {
		for (int i = 0;i < nVoxels[0];i++) {
#define VERT(vx,vy) ((vx)+(vy)*nx)
			voxels[offset(i, j)] = HVoxel(
				this, 
				points[VERT(i, j)],
				points[VERT(i+1, j)],
				points[VERT(i, j+1)],
				points[VERT(i+1, j+1)]
			);
#undef VERT
		}
	}
}


Heightfield2::~Heightfield2() {
	delete[] uvs;
	delete[] voxels;//FreeAligned(voxels);
	delete[] points;//FreeAligned(points);
	delete[] z;
}


BBox Heightfield2::ObjectBound() const {
	float minz = z[0], maxz = z[0];
	for (int i = 1; i < nx*ny; ++i) {
		if (z[i] < minz) minz = z[i];
		if (z[i] > maxz) maxz = z[i];
	}
	return BBox(Point(0, 0, minz), Point(1, 1, maxz));
}


bool Heightfield2::CanIntersect() const {
	return true;
}


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const ParamSet &params) {
	int nu = params.FindOneInt("nu", -1);
	int nv = params.FindOneInt("nv", -1);
	int nitems;
	const float *Pz = params.FindFloat("Pz", &nitems);
	Assert(nitems == nu*nv);
	Assert(nu != -1 && nv != -1 && Pz != NULL);
	return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


bool Heightfield2::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg) const {
	// Check ray against overall grid bounds
	float rayT;
	if (bound.Inside(ray(ray.mint))) {
		rayT = ray.mint;
	}
	else if (!bound.IntersectP(ray, &rayT))
	{
		return false;
	}
	Point gridIntersect = ray(rayT);
	// Set up 2D DDA for ray
	float NextCrossingT[2], DeltaT[2];
	int Step[2], Out[2], Pos[2];
	Ray rayObj = (*WorldToObject)(ray);
	Point gridIntersectObj = (*WorldToObject)(gridIntersect);
	for (int axis = 0; axis < 2; ++axis) {
		int v = gridIntersectObj[axis] * invWidth[axis];
		Pos[axis] = Clamp(v, 0, nVoxels[axis] - 1);
		// Compute current voxel for axis
		if (rayObj.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			float l = (Pos[axis] + 1) * width[axis] - gridIntersectObj[axis];
			NextCrossingT[axis] = rayT + l / rayObj.d[axis];
			DeltaT[axis] = width[axis] / rayObj.d[axis];
			Step[axis] = 1;
			Out[axis] = nVoxels[axis];
		}
		else {
			// Handle ray with negative direction for voxel stepping
			float l = (Pos[axis]) * width[axis] - gridIntersectObj[axis];
			NextCrossingT[axis] = rayT + l / rayObj.d[axis];
			DeltaT[axis] = -width[axis] / rayObj.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}
	// Walk ray through voxel grid
	for (;;) {

		// Check for intersection in current voxel and advance to next
		HVoxel& voxel = voxels[offset(Pos[0], Pos[1])];
		int stepAxis;
		float LimitT;
		// find LimitT and _stepAxis_
		if (NextCrossingT[0] < NextCrossingT[1]) {
			stepAxis = 0;
		}
		else {
			stepAxis = 1;
		}
		LimitT = NextCrossingT[stepAxis];
		if (voxel.Intersect(ray, tHit, rayEpsilon, dg, rayT, min(LimitT, ray.maxt))) {
			return true;
		}
		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		if (ray.maxt < LimitT) {
			//printf("Break due to limit\n");
			break;
		}

		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;
		rayT = NextCrossingT[stepAxis];
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return false;
}


bool Heightfield2::IntersectP(const Ray &ray) const {
	float rayT;
	if (bound.Inside(ray(ray.mint))) {
		rayT = ray.mint;
	}
	else if (!bound.IntersectP(ray, &rayT))
	{
		PBRT_GRID_RAY_MISSED_BOUNDS();
		return false;
	}
	Point gridIntersect = ray(rayT);
	// Set up 2D DDA for ray
	float NextCrossingT[2], DeltaT[2];
	int Step[2], Out[2], Pos[2];
	Ray rayObj = (*WorldToObject)(ray);
	Point gridIntersectObj = (*WorldToObject)(gridIntersect);

	for (int axis = 0; axis < 2; ++axis) {
		int v = gridIntersectObj[axis] * invWidth[axis];
		Pos[axis] = Clamp(v, 0, nVoxels[axis] - 1);
		// Compute current voxel for axis
		if (rayObj.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			float l = (Pos[axis] + 1) * width[axis] - gridIntersectObj[axis];
			NextCrossingT[axis] = rayT + l / rayObj.d[axis];
			DeltaT[axis] = width[axis] / rayObj.d[axis];
			Step[axis] = 1;
			Out[axis] = nVoxels[axis];
		}
		else {
			// Handle ray with negative direction for voxel stepping
			float l = (Pos[axis]) * width[axis] - gridIntersectObj[axis];
			NextCrossingT[axis] = rayT + l / rayObj.d[axis];
			DeltaT[axis] = -width[axis] / rayObj.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}
	// Walk ray through voxel grid
	for (;;) {
		// Check for intersection in current voxel and advance to next
		HVoxel& voxel = voxels[offset(Pos[0], Pos[1])];
		int stepAxis;
		float LimitT;
		// find LimitT and _stepAxis_
		if (NextCrossingT[0] < NextCrossingT[1]) {
			stepAxis = 0;
		}
		else {
			stepAxis = 1;
		}
		LimitT = NextCrossingT[stepAxis];
		if (voxel.IntersectP(ray, rayT, min(LimitT, ray.maxt))) {
			return true;
		}
		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		if (ray.maxt < LimitT)
			break;
		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;
		rayT = NextCrossingT[stepAxis];
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return false;
}



