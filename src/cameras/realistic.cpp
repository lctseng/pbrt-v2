
#include "stdafx.h"
#include "cameras/realistic.h"
#include <fstream>
#include <sstream>

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0) x0 = x1 = -0.5 * b / a;
	else {
		float q = (b > 0) ?
			-0.5 * (b + sqrt(discr)) :
			-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}
	return true;
}


RealisticLen::RealisticLen(float radius, float axisPos, float n, float aperture)
:radius(radius), n(n), axisPos(axisPos), aperture(aperture){
	float apertureRadius = aperture / 2.f;
	apertureRadius2 = apertureRadius * apertureRadius;
}

bool RealisticLen::Intersect(const Ray& ray, float* t) const {
	// perform ray-sphere intersection test
	// ref: http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
	float tFinal;
	if (radius != 0.f) {
		// a sphere surface
		// analytic solution
		float t0, t1;
		Vector L = ray.o - centerPoint;
		float a = Dot(ray.d, ray.d);
		float b = 2 * Dot(ray.d, L);
		float c = Dot(L, L) - radius * radius;
		if (!solveQuadratic(a, b, c, t0, t1)) {
			return false;
		}
		// t0 <= t1
		if (t0 > t1) {
			swap(t0, t1);
		}
		// if radius is neg, use t0
		if (radius < 0.f) {
			tFinal = t0;
		}
		else {
			tFinal = t1;
		}
	}
	else {
		// is a pure plane
		tFinal = (axisPoint.z - ray.o.z) / ray.d.z;
	}
	// check aperture
	Point pIntersect = ray(tFinal);
	float disToCenterSquare = pIntersect.x * pIntersect.x + pIntersect.y * pIntersect.y;
	if (disToCenterSquare > (apertureRadius2)) {
		// this ray is blocked by aperture
		return false;
	}
	*t = tFinal;
	return true;
}

void RealisticCamera::ReadLens(const string& filename) {
	string line;
	std::ifstream fin(filename);
	while (getline(fin, line))
	{
		std::istringstream ss(line);
		if (line.find('#') !=  std::string::npos) {
			// find #
			continue;
		}
		float radius, axisPos, n, aperture;
		ss >> radius;
		ss >> axisPos;
		ss >> n;
		if (n == 0.f) {
			n = 1.f;
		}
		ss >> aperture;
		if (radius == 0 && aperture_diameter < aperture)
		{
			aperture = aperture_diameter;
		}
		lens.push_front(RealisticLen(radius, axisPos, n, aperture));
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
	// setup center point for sphere lens
	float currentDistance = 0;
	for (auto& len : lens) {
		currentDistance += len.axisPos;
		len.axisToFilm = currentDistance;
	}
	filmPlaneZ = -currentDistance;
	for (auto& len : lens) {
		len.axisPoint = Point(0, 0, filmPlaneZ + len.axisToFilm);
		if (len.radius != 0.f) {
			len.centerPoint = Point(0, 0, len.axisPoint.z - len.radius);
		}
	}
	// compute effective resolution of x and y
	float aspectRatio = f->xResolution / f->yResolution;
	yRes = sqrt((filmdiag*filmdiag) / (1.0 + (aspectRatio*aspectRatio)));
	xRes = aspectRatio * yRes;
	// compute first len area for weight
	firstLenArea =  pow(lens.front().aperture / 2.f, 2.f) * M_PI;

	
}

Point RealisticCamera::RasterToCamera(const Point& p) const {
	float x = xRes - p.x*xRes / (float)film->xResolution - xRes / 2.f;
	float y = p.y * yRes / (float)film->yResolution - yRes / 2.f;
	float z = filmPlaneZ;
	return Point(x, y, z);

}

bool RealisticCamera::ProcessSnellsLaw(const RealisticLen& len, const Point& pIntersect, float n1, float n2, const Ray& inRay, Ray* outRay) const {
	// implement using Heckbert¡¦s method
	float n = n1 / n2;
	// compute normal of surface
	Vector N;
	if(len.radius < 0.f){
		N = Normalize(pIntersect - len.centerPoint);
	}
	else if (len.radius > 0.f) {
		N = Normalize(len.centerPoint - pIntersect);
	}
	else {
		// for planar, N is neg Z
		N = Vector(0, 0, -1);
	}	

	float c1 = -Dot(inRay.d, N);
	float c2Squared = 1.f - (n*n)*(1.f - (c1*c1));
	if (c2Squared < 0.f) {
		return false;
	}
	float c2 = sqrt(c2Squared);
	*outRay = Ray(pIntersect, Normalize(n * inRay.d + (n*c1 - c2) * N), 0.f, INFINITY);
	
	return true;
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

	// compute z coordinate of disk
	// using triangle formula: a2 + b2 = c2
	float zDisk;
	if(lens.front().radius != 0.f){
		float d = sqrt(pow(lens.front().radius, 2.f) - pow(halfAperture, 2.f));
		if (lens.front().radius < 0.f) {
			zDisk = lens.front().centerPoint.z - d;
		}
		else {
			zDisk = lens.front().centerPoint.z + d;
		}
	}
	else {
		zDisk = lens.front().axisPoint.z;
	}


	// create a initial ray
	// from film to disk
	Point pointOnDisk(xDisk, yDisk, zDisk);
	Vector filmToDisk = Normalize(pointOnDisk - pointOnFlim);
	Ray lenRay(pointOnFlim, filmToDisk, 0.f, INFINITY);
	
	float currentN, nextN;
	// traveral around the len set
	// from inner to outer
	for (int i = 0;i < lens.size();i++) {
		currentN = lens[i].n;
		if (i == lens.size() - 1) {
			// this is last len, shoot into air
			nextN = 1.f;
		}
		else {
			nextN = lens[i + 1].n;
		}
		float intersectT;
		if (lens[i].Intersect(lenRay, &intersectT)) {
			Point pIntersect = lenRay(intersectT);
			Ray outRay;
			if (ProcessSnellsLaw(lens[i], pIntersect, currentN, nextN, lenRay, &outRay)) {
				lenRay = outRay;
			}
			else {
				ray = NULL;
				return 0.f;
			}
		}
		else {
			// block by this len
			ray = NULL;
			return 0.f;
		}
	}
	*ray = CameraToWorld(lenRay);
	ray->d = Normalize(ray->d);
	// compute weight
	float cosTheta = Dot(Vector(0,0,1), filmToDisk);
	return (firstLenArea / pow(abs(filmPlaneZ - zDisk), 2.f)) * pow(cosTheta, 4.f);
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
