
#include "stdafx.h"
#include "cameras/realistic.h"
#include <fstream>
#include <sstream>

RealisticLen::RealisticLen(float radius, float n, float axisPos, float aperture)
:radius(radius), n(n), axisPos(axisPos), aperture(aperture){
	
}

void RealisticCamera::ReadLens(const string& filename) {
	string line;
	std::ifstream fin(filename);
	while (getline(fin, line)) {
		// filter comments
		int firstChar = line.find_first_not_of(' ');
		if (line[firstChar] != '#') {
			// form sstream
			std::stringstream ss;
			ss << line;
			// read 4 data
			float radius, axisPos, n, aperture;
			ss >> radius >> axisPos >> n >> aperture;
			// fix for camera aperture diameter
			if (radius == 0 && aperture < aperture_diameter) {
				aperture = aperture_diameter;
			}
			// fix for zero n
			if (n == 0) {
				n = 1.f;
			}
			// create new len
			lens.emplace_front(radius, axisPos, n, aperture);
		}
	}
}


RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
	float hither, float yon,
	float sopen, float sclose,
	float filmdistance, float aperture_diameter, string specfile,
	float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f), // pbrt-v2 doesnot specify hither and yon
	filmdistance(filmdistance), aperture_diameter(aperture_diameter)
{
	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.
	ReadLens(specfile);
	// film distance is the axisPos of nearest len
	lens.front().axisPos = filmdistance;
	// compute effective resolution of x and y
	float aspectRatio = f->xResolution / f->yResolution;
	yRes = sqrt((filmdiag*filmdiag) / (1.0 + (aspectRatio*aspectRatio)));
	xRes = aspectRatio * yRes;
	
}

Point RealisticCamera::RasterToCamera(const Point& p) const {
	float x = xRes - p.x*xRes / (float)film->xResolution - xRes / 2.f;
	float y = p.y * yRes / (float)film->yResolution - yRes / 2.f;
	float z = filmdistance * -1; // we treat the first len as plane Z
	return Point(x, y, z);

}

Ray RealisticCamera::ProcessSnellsLaw(const RealisticLen& len, float n1, float n2, const Ray& inRay) const {
	return inRay;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// YOUR CODE HERE -- make that ray!

	// use sample->imageX and sample->imageY to get raster-space coordinates
	// of the sample point on the film.
	// use sample->lensU and sample->lensV to get a sample position on the lens
	// Generate raster and camera samples
	
	// we need to map the sample xy into xRes and yRes
	// this is the ray start point (in screen space)
	Point pointOnFlim = RasterToCamera(Point(sample.imageX, sample.imageY, 0));
	float xDisk, yDisk;


	// sample a point from disk
	ConcentricSampleDisk(sample.lensU, sample.lensV, &xDisk, &yDisk);

	// scale to aperture/2 of first len
	float halfAperture = lens.front().aperture / 2.f;
	xDisk *= halfAperture;
	yDisk *= halfAperture;

	// create a initial ray
	// from film to disk
	Ray lenRay(pointOnFlim, Vector(Point(xDisk, yDisk, 0) - pointOnFlim), 0.f, INFINITY);
	
	float currentN, nextN;
	currentN = 1.0f; // Air
	// traveral around the len set
	for (int i = 0;i < lens.size();i++) {
		nextN = lens[i].n;
		lenRay = ProcessSnellsLaw(lens[i], currentN, nextN, lenRay );
	}

	CameraToWorld(lenRay, ray);
	return 1.f;
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
	const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance != -1);
	if (specfile == "") {
		Severe("No lens spec file supplied!\n");
	}
	return new RealisticCamera(cam2world, hither, yon,
		shutteropen, shutterclose, filmdistance, fstop,
		specfile, filmdiag, film);
}
