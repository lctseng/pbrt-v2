
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


 // lights/infinite.cpp*
#include "stdafx.h"
#include "lights/MedianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"
#include "list"

// InfiniteAreaLight Utility Classes
struct MedianCutEnvironmentCube {
	// MedianCutEnvironmentCube Public Methods
	MedianCutEnvironmentCube(const MedianCutEnvironmentLight *l, const Scene *s,
		float t, bool cv, float pe)
		: light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
	Spectrum operator()(int, int, const Point &p, const Vector &w) {
		Ray ray(p, w, pEpsilon, INFINITY, time);
		if (!computeVis || !scene->IntersectP(ray))
			return light->Le(RayDifferential(ray));
		return 0.f;
	}
	const MedianCutEnvironmentLight *light;
	const Scene *scene;
	float time, pEpsilon;
	bool computeVis;
};



// InfiniteAreaLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
	delete distribution;
	delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
	const Spectrum &L, int ns, const string &texmap)
	: Light(light2world, ns) {
	int width = 0, height = 0;
	RGBSpectrum *texels = NULL;
	// Read texel data from _texmap_ into _texels_
	if (texmap != "") {
		texels = ReadImage(texmap, &width, &height);
		if (texels)
			for (int i = 0; i < width * height; ++i)
				texels[i] *= L.ToRGBSpectrum();
	}
	if (!texels) {
		width = height = 1;
		texels = new RGBSpectrum[1];
		texels[0] = L.ToRGBSpectrum();
	}
	
	radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
	
	// Initialize sampling PDFs for infinite area light

	// Compute scalar-valued image _img_ from environment map
	float filter = 1.f / max(width, height);
	float *img = new float[width*height];
	for (int v = 0; v < height; ++v) {
		float vp = (float)v / (float)height;
		float sinTheta = sinf(M_PI * float(v + .5f) / float(height));
		for (int u = 0; u < width; ++u) {
			float up = (float)u / (float)width;
			img[u + v*width] = radianceMap->Lookup(up, vp, filter).y();
			img[u + v*width] *= sinTheta;
		}
	}

	// Compute sampling distributions for rows and columns of image
	distribution = new Distribution2D(img, width, height);
	CreatePointLights(texels, img, ns, width, height);
	
	delete[] texels;
	delete[] img;
}

#ifdef _DEBUG
#define SHOW_DEBUG_INFO 1
#else
#define SHOW_DEBUG_INFO 0
#endif
#define INDEX_AT(i,j) ((i)*width+(j))

// find a area sum
float FindAreaSum(float* table, int width, int height, MedianRect& rect) {
	float bottomRight = 0.0, topLeft = 0.0, topRight = 0.0, bottomLeft = 0.0;
	// compute topLeft
	if (rect.topLeft.x > 0 && rect.topLeft.y > 0) {
		topLeft = table[INDEX_AT(rect.topLeft.y - 1, rect.topLeft.x - 1)];
	}
	// compute topRight
	if (rect.topRight.y > 0) {
		topRight = table[INDEX_AT(rect.topRight.y - 1, rect.topRight.x)];
	}
	// compute bottomLeft
	if (rect.bottomLeft.x > 0) {
		bottomLeft = table[INDEX_AT(rect.bottomLeft.y, rect.bottomLeft.x - 1)];
	}
	// compute bottomRight
	bottomRight = table[INDEX_AT(rect.bottomRight.y, rect.bottomRight.x)];

	return bottomRight - topRight - bottomLeft + topLeft;
}

// convert specturm to float
float RGBSpectrumToFloat(RGBSpectrum& s) {
	float rgb[3];
	s.ToRGB(rgb);
	return 0.2125 * rgb[0] + 0.7154 * rgb[1] + 0.0721 * rgb[2];
}
// find split axis among rect
// param: table, rect, axis
// return: position index on that axis
int FindSplitIndexInRect(float* table, int width, int height, MedianRect& rect, int axis) {
	int minIndex;
	if (axis == MEDIAN_AXIS_X) {
		// find a column to split
		float minDiff = INFINITY;
		minIndex = rect.topLeft.x;

		MedianRect leftRect = rect, rightRect = rect;
		for (int i = rect.topLeft.x + 1;i <= rect.topRight.x;i++) {
			// determine the right-x for leftRect
			leftRect.topRight.x = leftRect.bottomRight.x = i - 1;
			// determine the left-x for rightRect
			rightRect.topLeft.x = rightRect.bottomLeft.x = i;

			float leftArea = FindAreaSum(table, width, height, leftRect);
			float rightArea = FindAreaSum(table, width, height, rightRect);
			float areaDiff = abs(leftArea - rightArea);
			if (areaDiff < minDiff) {
				minDiff = areaDiff;
				minIndex = i;
			}
		}
	}
	else {
		// find a row to split
		// find a column to split
		float minDiff = INFINITY;
		minIndex = rect.topLeft.y;

		MedianRect topRect = rect, bottomRect = rect;
		for (int i = rect.topLeft.y + 1;i <= rect.bottomLeft.y;i++) {
			// determine the bottom-y for topRect
			topRect.bottomRight.y = topRect.bottomLeft.y = i - 1;
			// determine the top-y for bottomRect
			bottomRect.topLeft.y = bottomRect.topRight.y = i;

			float topArea = FindAreaSum(table, width, height, topRect);
			float bottomArea = FindAreaSum(table, width, height, bottomRect);
			float areaDiff = abs(topArea - bottomArea);
			if (areaDiff < minDiff) {
				minDiff = areaDiff;
				minIndex = i;
			}
		}
	}
	return minIndex;
}


void MedianCutEnvironmentLight::CreatePointLights(RGBSpectrum* pixels, float* img, int numberOfLights,  int width, int height) {



	// Remember to scale the light intensity with the areas (solid angles)
	float solidAngleScale = ((2.f * M_PI) / (width - 1)) * ((M_PI) / (height - 1));
	/*
	for (int v = 0; v < height; v++) {
		float sinTheta = sinf(M_PI * float(v + .5f) / float(height));
		for (int u = 0; u < width; u++)
			pixels[u + v*width] = pixels[u + v*width] * solidAngleScale * sinTheta;
	}
	*/


	// build summed area table
	float* raw_table = new float[width * height];
	float* summed_table = new float[width * height];
	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			//raw_table[INDEX_AT(i, j)] = RGBSpectrumToFloat(pixels[INDEX_AT(i, j)]) * solidAngleScale;
			raw_table[INDEX_AT(i, j)] = img[INDEX_AT(i, j)] * solidAngleScale;
		}
	}

	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			// sum from left to current
			int cur_index = INDEX_AT(i, j);
			summed_table[cur_index] = 0.f;
			for (int k = 0;k <= j;k++) {
				summed_table[cur_index] += raw_table[INDEX_AT(i, k)];
			}
			// sum the previous row
			if (i > 0) {
				int prev_index = INDEX_AT(i - 1, j);
				summed_table[cur_index] += summed_table[prev_index];
			}
		}
	}


	MedianRect initRegion;
	initRegion.topLeft = {0, 0};
	initRegion.topRight = { width - 1, 0 };
	initRegion.bottomLeft = { 0, height - 1 };
	initRegion.bottomRight = { width -1 , height - 1 };

	std::vector<MedianRect> regionsA;
	std::vector<MedianRect> regionsB;
	std::vector<MedianRect> *pWorking, *pResult;
	pWorking = &regionsA;
	pResult = &regionsB;
	pWorking->push_back(initRegion);
	while(true){
		// split all region
		pResult->clear();
		for (auto& region : *pWorking) {
			int rw = region.width(), rh = region.height();
			if (rw > rh) {
				// split at x
				if (rw > 1) {
					// do split
					int split_x = FindSplitIndexInRect(summed_table, width, height, region, MEDIAN_AXIS_X);
					MedianRect leftRegion = region, rightRegion = region;
					// determine the right-x for leftRegion
					leftRegion.topRight.x = leftRegion.bottomRight.x = split_x - 1;
					// determine the left-x for rightRegion
					rightRegion.topLeft.x = rightRegion.bottomLeft.x = split_x;

					pResult->push_back(leftRegion);
					pResult->push_back(rightRegion);
				}
				else {
					// store back again
					pResult->push_back(region);
				}
			}
			else {
				// split at y
				if (rh > 1) {
					// do split
					int split_y = FindSplitIndexInRect(summed_table, width, height, region, MEDIAN_AXIS_Y);
					MedianRect topRegion = region, bottomRegion = region;
					// determine the bottom-y for topRegion
					topRegion.bottomRight.y = topRegion.bottomLeft.y = split_y - 1;
					// determine the top-y for bottomRegion
					bottomRegion.topLeft.y = bottomRegion.topRight.y = split_y;

					pResult->push_back(topRegion);
					pResult->push_back(bottomRegion);
				}
				else {
					// store back again
					pResult->push_back(region);
				}
			}
		}
		// swap index if we have next
		if (pResult->size() < numberOfLights) {
			swap(pWorking, pResult);
		}
		else {
			break;
		}
	}
#if SHOW_DEBUG_INFO

	// copy pixel array for debug
	RGBSpectrum* clonePixels = new RGBSpectrum[width * height];

	for (int i = 0;i < width * height;i++) {
		clonePixels[i] = pixels[i];
	}

	printf("There are %d regions\n", pResult->size());

	for (auto& region : *pResult) {
		DrawRectOutlineOnSpectrum(clonePixels, width, height, region);
	}
	WriteSpectrumToFile("lightMap.exr", clonePixels, width, height);


	delete[] clonePixels;

#endif
	// region division is ready
	// create points light from each region
	
	// u, v must to scale to [0,1]
	float scale_u = 1.f / width, scale_v = 1.f / height;
	for (auto& region : *pResult) {
		// new point light
		MedianPointLight light;
		// simple implementation: put in the center
		light.uv[0] = (region.topLeft.x + region.topRight.x) * 0.5 * scale_u;
		light.uv[1] = (region.topLeft.y + region.bottomLeft.y) * 0.5 * scale_v;
		// sum the spectrum
		RGBSpectrum spectrum;
		for (int i = region.topLeft.y;i <= region.bottomLeft.y; i++) {
			for (int j = region.topLeft.x;j <= region.topRight.x;j++) {
				spectrum += pixels[INDEX_AT(i, j)];
			}
		}
		light.intensity = spectrum;
		// store new light
		pointLights.push_back(light);
	}

	lightPdf = 1.0f / numberOfLights;

	delete[] raw_table;
	delete[] summed_table;
}


void DrawRectOutlineOnSpectrum(RGBSpectrum* pixels, int width, int height, MedianRect& rect) {
	int line_width = 1;
	int halfLineWidth = line_width / 2;
	MedianRect dRect;
	// top line: horz
	dRect = rect;
	dRect.topLeft.y -= halfLineWidth;
	dRect.topRight.y -= halfLineWidth;
	dRect.bottomLeft.y = dRect.topLeft.y + line_width;
	dRect.bottomRight.y = dRect.topRight.y + line_width;
	DrawRectAreaOnSpectrum(pixels, width, height, dRect);
	// bottom line: horz
	dRect = rect;
	dRect.topLeft.y = dRect.bottomLeft.y - halfLineWidth;
	dRect.topRight.y = dRect.bottomRight.y - halfLineWidth;
	dRect.bottomLeft.y = dRect.topLeft.y + line_width;
	dRect.bottomRight.y = dRect.topRight.y + line_width;
	DrawRectAreaOnSpectrum(pixels, width, height, dRect);
	// left line: vert
	dRect = rect;
	dRect.topLeft.x -= halfLineWidth;
	dRect.bottomLeft.x -= halfLineWidth;
	dRect.topRight.x = dRect.topLeft.x + line_width;
	dRect.bottomRight.x = dRect.bottomLeft.x + line_width;
	DrawRectAreaOnSpectrum(pixels, width, height, dRect);
	// right line: vert
	dRect = rect;
	dRect.topLeft.x = dRect.topRight.x - halfLineWidth;
	dRect.bottomLeft.x = dRect.bottomRight.x - halfLineWidth;
	dRect.topRight.x = dRect.topLeft.x + line_width;
	dRect.bottomRight.x = dRect.bottomLeft.x + line_width;
	DrawRectAreaOnSpectrum(pixels, width, height, dRect);
}

void DrawRectAreaOnSpectrum(RGBSpectrum* pixels, int width, int height, MedianRect& rect) {
	float color[] = { 0, 1, 0 };
	for (int i = rect.topLeft.y;i <= rect.bottomLeft.y;i++) {
		for (int j = rect.topLeft.x; j <= rect.topRight.x;j++) {
			// range check
			if (i >= 0 && i < height && j >= 0 && j < width) {
				pixels[INDEX_AT(i, j)] = RGBSpectrum::FromRGB(color);
			}
		}
	}
}

void WriteSpectrumToFile(const string& filename, RGBSpectrum* pixels, int width, int height) {
	float* rgbArray = new float[width *height * 3];
	float* alphaArray = new float[width * height];
	// copy to rgbArray
	for (int i = 0;i < width * height;i++) {
		pixels[i].ToRGB(rgbArray + i * 3);
		alphaArray[i] = 1.0;
	}
	WriteImage(filename, rgbArray, alphaArray, width, height, width, height, 0, 0);
	printf("Image: %s written\n", filename.c_str());
	delete[] rgbArray;
	delete[] alphaArray;
}

#undef INDEX_AT



Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
	Point worldCenter;
	float worldRadius;
	scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
	return M_PI * worldRadius * worldRadius *
		Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
	Vector wh = Normalize(WorldToLight(r.d));
	float s = SphericalPhi(wh) * INV_TWOPI;
	float t = SphericalTheta(wh) * INV_PI;
	return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
	int lmax, const Scene *scene, bool computeLightVis,
	float time, RNG &rng, Spectrum *coeffs) const {
	// Project _InfiniteAreaLight_ to SH using Monte Carlo if visibility needed
	if (computeLightVis) {
		Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
			time, rng, coeffs);
		return;
	}
	for (int i = 0; i < SHTerms(lmax); ++i)
		coeffs[i] = 0.f;
	int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
	if (min(ntheta, nphi) > 50) {
		// Project _InfiniteAreaLight_ to SH from lat-long representation

		// Precompute $\theta$ and $\phi$ values for lat-long map projection
		float *buf = new float[2 * ntheta + 2 * nphi];
		float *bufp = buf;
		float *sintheta = bufp;  bufp += ntheta;
		float *costheta = bufp;  bufp += ntheta;
		float *sinphi = bufp;    bufp += nphi;
		float *cosphi = bufp;
		for (int theta = 0; theta < ntheta; ++theta) {
			sintheta[theta] = sinf((theta + .5f) / ntheta * M_PI);
			costheta[theta] = cosf((theta + .5f) / ntheta * M_PI);
		}
		for (int phi = 0; phi < nphi; ++phi) {
			sinphi[phi] = sinf((phi + .5f) / nphi * 2.f * M_PI);
			cosphi[phi] = cosf((phi + .5f) / nphi * 2.f * M_PI);
		}
		float *Ylm = ALLOCA(float, SHTerms(lmax));
		for (int theta = 0; theta < ntheta; ++theta) {
			for (int phi = 0; phi < nphi; ++phi) {
				// Add _InfiniteAreaLight_ texel's contribution to SH coefficients
				Vector w = Vector(sintheta[theta] * cosphi[phi],
					sintheta[theta] * sinphi[phi],
					costheta[theta]);
				w = Normalize(LightToWorld(w));
				Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
					SPECTRUM_ILLUMINANT);
				SHEvaluate(w, lmax, Ylm);
				for (int i = 0; i < SHTerms(lmax); ++i)
					coeffs[i] += Le * Ylm[i] * sintheta[theta] *
					(M_PI / ntheta) * (2.f * M_PI / nphi);
			}
		}

		// Free memory used for lat-long theta and phi values
		delete[] buf;
	}
	else {
		// Project _InfiniteAreaLight_ to SH from cube map sampling
		SHProjectCube(MedianCutEnvironmentCube(this, scene, time, computeLightVis,
			pEpsilon),
			p, 200, lmax, coeffs);
	}
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
	const ParamSet &paramSet) {
	Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
	Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
	string texmap = paramSet.FindOneFilename("mapname", "");
	int nSamples = paramSet.FindOneInt("nsamples", 1);
	if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
	return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
	const LightSample &ls, float time, Vector *wi, float *pdf,
	VisibilityTester *visibility) const {
	PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
	if (pointLights.empty()) {
		printf("===Warning: No light available!===");
		*pdf = 0.0f;
		return 0;
	}
	// pick a light
	auto& light = pointLights[rand() % pointLights.size()];

	// Convert infinite light sample point to direction
	float theta = light.uv[1] * M_PI, phi = light.uv[0] * 2.f * M_PI;
	float costheta = cosf(theta), sintheta = sinf(theta);
	float sinphi = sinf(phi), cosphi = cosf(phi);
	*wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
		costheta));

	// Set PDF
	*pdf = lightPdf;

	// Set visibility
	visibility->SetRay(p, pEpsilon, *wi, time);

	// Return radiance value for infinite light direction	
	Spectrum Ls = Spectrum(light.intensity,
		SPECTRUM_ILLUMINANT);
	PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();

	return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
	PBRT_INFINITE_LIGHT_STARTED_PDF();
	Vector wi = WorldToLight(w);
	float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
	float sintheta = sinf(theta);
	if (sintheta == 0.f) return 0.f;
	float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
		(2.f * M_PI * M_PI * sintheta);
	PBRT_INFINITE_LIGHT_FINISHED_PDF();
	return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
	const LightSample &ls, float u1, float u2, float time,
	Ray *ray, Normal *Ns, float *pdf) const {
	PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
	// Compute direction for infinite light sample ray

	// Find $(u,v)$ sample coordinates in infinite light texture
	float uv[2], mapPdf;
	distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
	if (mapPdf == 0.f) return Spectrum(0.f);

	float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
	float costheta = cosf(theta), sintheta = sinf(theta);
	float sinphi = sinf(phi), cosphi = cosf(phi);
	Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
		costheta));
	*Ns = (Normal)d;

	// Compute origin for infinite light sample ray
	Point worldCenter;
	float worldRadius;
	scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
	Vector v1, v2;
	CoordinateSystem(-d, &v1, &v2);
	float d1, d2;
	ConcentricSampleDisk(u1, u2, &d1, &d2);
	Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
	*ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

	// Compute _InfiniteAreaLight_ ray PDF
	float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
	float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
	*pdf = directionPdf * areaPdf;
	if (sintheta == 0.f) *pdf = 0.f;
	Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
	PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
	return Ls;
}


