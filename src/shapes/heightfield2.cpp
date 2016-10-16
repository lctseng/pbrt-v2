
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

Plane3D::Plane3D(const Point& p1, const Point& p2, const Point& p3) {
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;
	//printf("Plane: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", p1->x, p1->y, p1->z, p2->x, p2->y, p2->z, p3->x, p3->y, p3->z);
	// compute normal
	Vector v1 = points[1] - points[0];
	Vector v2 = points[2] - points[0];
	normal = Cross(v1, v2);
}
bool Plane3D::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg, float minT, float maxT) {
	float t, lower, upper;
	// compute intersection
	// compute t, formula: t = (n1*(x0-xr) + n2*(y0-yr) + n3*(z0-zr))/(n1*a+n2*b+n3*c)
	// n1*a+n2*b+n3*c, this should be non-zero
	lower = normal.x*ray.d.x + normal.y*ray.d.y + normal.z*ray.d.z;
	if (abs(lower) < 0.0005f) {
		return false;
	}
	upper = normal.x*(points[0].x - ray.o.x) +
		normal.y*(points[0].y - ray.o.y) +
		normal.z*(points[0].z - ray.o.z);
	t = upper / lower;
	// check t range
	if (t < minT || t > maxT) {
		return false;
	}
	// already intersected, test it with the rendered triangle

	//*tHit = t;
	//*rayEpsilon = 1e-3f * *tHit;

	return true;
}


HVoxel::HVoxel(const Point& p1, const Point& p2, const Point& p3, const Point& p4)
{
	fromTriangle = false;
	planes[0] = Plane3D(p1, p2, p4);
	planes[1] = Plane3D(p1, p4, p3);
}
void HVoxel::SetTriangle(TriangleMesh* m) {
	fromTriangle = true;
	triangle = m;
	m->Refine(refined);
}
bool HVoxel::IntersectP(const Ray &ray) {
	if (fromTriangle) {
		for (auto& shape : refined) {
			if (shape->IntersectP(ray)) {
				return true;
			}
		}
		return false;
	}
	return false;
}
bool HVoxel::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg, float minT, float maxT) {
	for (auto& plane : planes) {
		if (plane.Intersect(ray, tHit, rayEpsilon, dg, minT, maxT)) {
			for (auto& shape : refined) {
				if (shape) {
					if (shape->Intersect(ray, tHit, rayEpsilon, dg)) {
						return true;
					}
				}
			}
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
			int ntris = 2;
			int verts[] = {
				VERT(i, j),
				VERT(i + 1, j),
				VERT(i + 1, j + 1),

				VERT(i, j),
				VERT(i + 1, j + 1),
				VERT(i, j + 1)
			};
			ParamSet paramSet;
			paramSet.AddInt("indices", verts, 3 * ntris);
			paramSet.AddFloat("uv", uvs, 2 * np);
			paramSet.AddPoint("P", points, np);
			auto t = CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet);
			int vindex = offset(i, j);
			voxels[vindex] = HVoxel(
				(*ObjectToWorld)(points[VERT(i, j)]),
				(*ObjectToWorld)(points[VERT(i + 1, j)]),
				(*ObjectToWorld)(points[VERT(i, j + 1)]),
				(*ObjectToWorld)(points[VERT(i + 1, j + 1)])
			);
			voxels[vindex].SetTriangle(t);
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


void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
	/*
	int ntris = 2*(nx-1)*(ny-1);
	refined.reserve(ntris);
	int *verts = new int[3*ntris];
	Point *P = new Point[nx*ny];
	float *uvs = new float[2*nx*ny];
	int nverts = nx*ny;
	int x, y;
	// Compute heightfield vertex positions
	int pos = 0;
	for (y = 0; y < ny; ++y) {
		for (x = 0; x < nx; ++x) {
			P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
			P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
			P[pos].z = z[pos];
			++pos;
		}
	}

	// Fill in heightfield vertex offset array
	int *vp = verts;
	for (y = 0; y < ny-1; ++y) {
		for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
			*vp++ = VERT(x, y);
			*vp++ = VERT(x+1, y);
			*vp++ = VERT(x+1, y+1);

			*vp++ = VERT(x, y);
			*vp++ = VERT(x+1, y+1);
			*vp++ = VERT(x, y+1);
		}
#undef VERT
	}
	ParamSet paramSet;
	paramSet.AddInt("indices", verts, 3*ntris);
	paramSet.AddFloat("uv", uvs, 2 * nverts);
	paramSet.AddPoint("P", P, nverts);
	refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
	delete[] P;
	delete[] uvs;
	delete[] verts;
	*/
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
	//printf("Ray input d: (%f, %f, %f)\n", ray.d.x, ray.d.y, ray.d.z);
	float rayT;
	if (bound.Inside(ray(ray.mint))) {
		//printf("Ray inside!\n");
		rayT = ray.mint;
	}
	else if (!bound.IntersectP(ray, &rayT))
	{
		//printf("Ray out of bound!\n");
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
		Assert(Pos[axis] >= 0);
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
		//rayT = NextCrossingT[stepAxis];
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return false;
	/*
	PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
	// Compute $\VEC{s}_1$

	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
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
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
		return false;

	// Compute triangle partial derivatives
	Vector dpdu, dpdv;
	float uvs[3][2];
	GetUVs(uvs);

	// Compute deltas for triangle partial derivatives
	float du1 = uvs[0][0] - uvs[2][0];
	float du2 = uvs[1][0] - uvs[2][0];
	float dv1 = uvs[0][1] - uvs[2][1];
	float dv2 = uvs[1][1] - uvs[2][1];
	Vector dp1 = p1 - p3, dp2 = p2 - p3;
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

	// Interpolate $(u,v)$ triangle parametric coordinates
	float b0 = 1 - b1 - b2;
	float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
	float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];

	// Test intersection against alpha texture, if present
	if (ray.depth != -1) {
		if (mesh->alphaTexture) {
			DifferentialGeometry dgLocal(ray(t), dpdu, dpdv,
				Normal(0, 0, 0), Normal(0, 0, 0),
				tu, tv, this);
			if (mesh->alphaTexture->Evaluate(dgLocal) == 0.f)
				return false;
		}
	}

	// Fill in _DifferentialGeometry_ from triangle hit
	*dg = DifferentialGeometry(ray(t), dpdu, dpdv,
		Normal(0, 0, 0), Normal(0, 0, 0),
		tu, tv, this);
	*tHit = t;
	*rayEpsilon = 1e-3f * *tHit;
	PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
	*/
	return true;

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
		if (voxel.IntersectP(ray)) {
			return true;
		}
		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		if (ray.maxt < LimitT)
			break;
		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;
		//rayT = NextCrossingT[stepAxis];
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return false;
}



