
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
	CreatePointLights(texels, width, height);
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
    delete[] texels;
    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
        }
    }

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
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
	return rgb[0] + rgb[1] + rgb[2];
	// return 0.2125 * rgb[0] + 0.7154 * rgb[1] + 0.0721 * rgb[2];
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
			rightRect.topLeft.x = rightRect.bottomLeft.x = i ;
			
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


void MedianCutEnvironmentLight::CreatePointLights(RGBSpectrum* pixels, int width, int height) {
	

	// build summed area table
	float* raw_table = new float[width * height];
	float* summed_table = new float[width * height];
	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			raw_table[INDEX_AT(i, j)] = RGBSpectrumToFloat(pixels[INDEX_AT(i, j)]);
		}
	}

	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			// sum from left to current
			int cur_index = INDEX_AT(i,j);
			summed_table[cur_index] = 0.f;
			for (int k = 0;k <= j;k++) {
				summed_table[cur_index] += raw_table[INDEX_AT(i,k)];
			}
			// sum the previous row
			if (i > 0) {
				int prev_index = INDEX_AT(i-1, j);
				summed_table[cur_index] += summed_table[prev_index];
			}
		}
	}
#if SHOW_DEBUG_INFO

	printf("Partial row image:\n");
	for (int i = 0;i < 10;i++) {
		for (int j = 0;j < 10;j++) {
			printf("%.6f ", raw_table[INDEX_AT(i, j)]);
		}
		printf("\n");
	}

	printf("Partial summed table:\n");
	for (int i = 0;i < 10;i++) {
		for (int j = 0;j < 10;j++) {
			printf("%.6f ", summed_table[INDEX_AT(i,j)]);
		}
		printf("\n");
	}
	// Verifying
	printf("Verifying table\n");
	float br = summed_table[INDEX_AT(3, 4)];
	float tl = summed_table[INDEX_AT(1, 1)];
	float tr = summed_table[INDEX_AT(1, 4)];
	float bl = summed_table[INDEX_AT(3, 1)];
	float ans = raw_table[INDEX_AT(2, 2)] + raw_table[INDEX_AT(2, 3)] + raw_table[INDEX_AT(2, 4)] +
		raw_table[INDEX_AT(3, 2)] + raw_table[INDEX_AT(3, 3)] + raw_table[INDEX_AT(3, 4)];
	printf("Predicted: %.6f, Real ans: %.6f\n", br - tr - bl + tl , ans);
	

	// Verifying SummedArea Function
	MedianRect rect;
	rect.topLeft = {2, 2};
	rect.topRight = { 4, 2 };
	rect.bottomLeft = {2, 3};
	rect.bottomRight = { 4, 3 };
	printf("Predict Funcrion: %.6f, Real ans: %.6f\n", FindAreaSum(summed_table, width, height, rect), ans);

	// Verifying Split Function
	// whole
	rect.topLeft = { 0 ,0 };
	rect.topRight = { 9 ,0 };
	rect.bottomLeft = { 0 ,9 };
	rect.bottomRight = { 9 ,9 };
	printf("Split whole at x: %d\n", FindSplitIndexInRect(summed_table, width, height,rect, MEDIAN_AXIS_X));
	printf("Split whole at y: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_Y));
	
	// partial x
	rect.topLeft = { 4 ,0 };
	rect.topRight = { 7 ,0 };
	rect.bottomLeft = { 4 ,9 };
	rect.bottomRight = { 7 ,9 };
	printf("Split partial-x at x: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_X));
	printf("Split partial-x at y: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_Y));

	// partial y
	rect.topLeft = { 0 ,4 };
	rect.topRight = { 9 ,4 };
	rect.bottomLeft = { 0 ,7 };
	rect.bottomRight = { 9 ,7 };
	printf("Split partial-y at x: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_X));
	printf("Split partial-y at y: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_Y));

	// partial both
	rect.topLeft = { 4 ,1 };
	rect.topRight = { 7 ,1 };
	rect.bottomLeft = { 4 ,4 };
	rect.bottomRight = { 7 ,4 };
	printf("Split partial-both at x: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_X));
	printf("Split partial-both at y: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_Y));

	// boundary-test
	rect.topLeft = { 3 ,1 };
	rect.topRight = { 3 ,1 };
	rect.bottomLeft = { 3 ,4 };
	rect.bottomRight = { 3 ,4 };
	printf("Split boundary-test at x: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_X));
	rect.topLeft = { 4 ,4 };
	rect.topRight = { 7 ,4 };
	rect.bottomLeft = { 4 ,4 };
	rect.bottomRight = { 7 ,4 };
	printf("Split boundary-test at y: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_Y));

	// corner-test
	rect.topLeft = { 5 ,6 };
	rect.topRight = { 5 ,6 };
	rect.bottomLeft = { 5 ,6 };
	rect.bottomRight = { 5 ,6 };
	printf("Split cornor-test at x: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_X));
	printf("Split cornor-test at y: %d\n", FindSplitIndexInRect(summed_table, width, height, rect, MEDIAN_AXIS_Y));
	
	// copy pixel array for debug
	RGBSpectrum* clonePixels = new RGBSpectrum[width * height];
	
	for (int i = 0;i < width * height;i++) {
		clonePixels[i] = pixels[i];
	}

	WriteSpectrumToFile("lightMap.exr", pixels, width, height);


	delete[] clonePixels;

#endif



	delete[] raw_table;
	delete[] summed_table;
}


void WriteSpectrumToFile(const string& filename, RGBSpectrum* pixels, int width, int height) {
	float* rgbArray = new float[width *height * 3];
	// copy to rgbArray
	for (int i = 0;i < width * height;i++) {
		pixels[i].ToRGB(rgbArray+i*3);
	}
	WriteImage(filename, rgbArray, rgbArray, width, height, width, height, 0 , 0);
	delete[] rgbArray;
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
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
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
    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return 0.f;

    // Convert infinite light sample point to direction
    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    // Compute PDF for sampled infinite light direction
    *pdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    if (sintheta == 0.f) *pdf = 0.f;

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    Spectrum Ls = Spectrum(radianceMap->Lookup(uv[0], uv[1]),
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


