
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"
#include <deque>



// len info
class RealisticLen {
public:
	RealisticLen() {};
	RealisticLen(float radius, float n, float axisPos, float aperture);
	float radius, n, axisPos, aperture;
private:
};


// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	typedef std::deque<RealisticLen> LenSet;
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
 

	Point RasterToCamera(const Point& p) const;

private:
	// RealisticCamera Public Methods
	LenSet lens;
	float xRes, yRes;
	float filmdistance, aperture_diameter;

	void ReadLens(const string& filename);
	Ray ProcessSnellsLaw(const RealisticLen& len, float n1, float n2, const Ray& inRay) const;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H