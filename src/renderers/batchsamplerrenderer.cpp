
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


// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/batchsamplerrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t hash = 0, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
} 

// BatchSamplerRendererTask Definitions
void BatchSamplerRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _BatchSamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();

	BatchSamplerRendererQueue queue(this, maxSamples, rng);
	Sample* samples;
	RayDifferential* rays;
	float* rayWeights;

    // Get samples from _Sampler_ and update image
    int sampleCount;
	while (samples = queue.RequireSampleSpace(),(sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
		rays = queue.RequireRaySpace();
		rayWeights = queue.RequireRayWeightSpace();
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
			rayWeights[i] = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
        }
		queue.CommitLiRequest(sampleCount);
    }
	queue.Flush();

    // Clean up after _BatchSamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// BatchSamplerRenderer Method Definitions
BatchSamplerRenderer::BatchSamplerRenderer(Sampler *s, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 bool visIds) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


BatchSamplerRenderer::~BatchSamplerRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


void BatchSamplerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

    // Create and launch _BatchSamplerRendererTask_s for rendering image

    // Compute number of _BatchSamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
	nTasks = 1;//nPixels / 10000;
	nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i)
        renderTasks.push_back(new BatchSamplerRendererTask(scene, this, camera,
                                                      reporter, sampler, sample, 
                                                      visualizeObjectIds, 
                                                      nTasks-1-i, nTasks));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();
}


Spectrum BatchSamplerRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
	Assert(ray.time == sample->time);
	Assert(!ray.HasNaNs());
	// Allocate local variables for _isect_ and _T_ if needed
	Spectrum localT;
	if (!T) T = &localT;
	Intersection localIsect;
	if (!isect) isect = &localIsect;
	Spectrum Li = 0.f;
	if (scene->Intersect(ray, isect))
		Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
			rng, arena);
	else {
		// Handle ray that doesn't intersect any geometry
		for (uint32_t i = 0; i < scene->lights.size(); ++i)
			Li += scene->lights[i]->Le(ray);
	}
	Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
		T, arena);
	return *T * Li + Lvi;
}


Spectrum BatchSamplerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}

BatchSamplerRendererQueue::BatchSamplerRendererQueue(BatchSamplerRendererTask* render_task, int maxSample, RNG& rng):
render_task(render_task),maxSample(maxSample),rng(rng)
{
	int minSpace = max(maxSample, BATCH_RENDER_SIZE);

	renderer = nullptr;
	renderer = dynamic_cast<const BatchSamplerRenderer*>(render_task->renderer);
	assert(renderer != nullptr);

	samples = render_task->origSample->Duplicate(minSpace);
	rays = new RayDifferential[minSpace];
	isects = new Intersection[minSpace];
	rayWeights = new float[minSpace];
	hits = new bool[minSpace];
	arenas = new MemoryArena[minSpace];
	Lis = new Spectrum[minSpace];
	Lvis = new Spectrum[minSpace];
	Ls = new Spectrum[minSpace];
	Ts = new Spectrum[minSpace];
	taskNum = 0;
}

BatchSamplerRendererQueue::~BatchSamplerRendererQueue() {
	delete[] samples;
	delete[] rays;
	delete[] isects;
	delete[] rayWeights;
	delete[] hits;
	delete[] arenas;
	delete[] Lis;
	delete[] Lvis;
	delete[] Ls;
	delete[] Ts;
}

RayDifferential* BatchSamplerRendererQueue::RequireRaySpace() {
	return rays + taskNum;
}
Intersection* BatchSamplerRendererQueue::RequireIntersectionSpace() {
	return isects + taskNum;
}
Sample* BatchSamplerRendererQueue::RequireSampleSpace() {
	return samples + taskNum;
}
float* BatchSamplerRendererQueue::RequireRayWeightSpace() {
	return rayWeights + taskNum;
}
bool* BatchSamplerRendererQueue::RequireHitSpace() {
	return hits + taskNum;
}

void BatchSamplerRendererQueue::CommitLiRequest(int count) {
	taskNum += count;
	if (taskNum >= BATCH_RENDER_SIZE) {
		LaunchLiProcess();
	}
}
void BatchSamplerRendererQueue::WaitForRequest() {
	
}

void BatchSamplerRendererQueue::Flush() {
	LaunchLiProcess();
	WaitForRequest();
}
void BatchSamplerRendererQueue::LaunchLiProcess() {
	if (taskNum > 0) {
		LaunchPrimaryRayIntersection();
		// do no-weight ray
		LaunchNoWeightRayProcess();
		// do Li for hit
		LaunchSurfaceIntegration();
		// do Li for NotHit
		LaunchMissedRayProcess();
		// volume integration for weighted rays
		LaunchVolumeIntegration();
		// combine
		LaunchCombineProcess();
		// add sample and clean up
		for (int i = 0;i < taskNum;i++) {
			render_task->camera->film->AddSample(samples[i], Ls[i]);
			arenas[i].FreeAll();
		}
		taskNum = 0;
	}
}



void BatchSamplerRendererQueue::LaunchPrimaryRayIntersection() {
	BEGIN_TIMING(LaunchPrimaryRayIntersection);
	for (int i = 0;i < taskNum;i++) {
		hits[i] = render_task->scene->Intersect(rays[i], &isects[i]);
	}
	END_TIMING(LaunchPrimaryRayIntersection);
}

void BatchSamplerRendererQueue::LaunchSurfaceIntegration() {
	BEGIN_TIMING(LaunchSurfaceIntegration);
	for (int i = 0;i < taskNum;i++) {
		if (rayWeights[i] > 0.f) {
			if (hits[i]) {
				Spectrum& Li = Lis[i];
				Li = renderer->surfaceIntegrator->Li(render_task->scene, renderer, rays[i], isects[i], &samples[i],
					rng, arenas[i]);
			}
		}
	}
	END_TIMING(LaunchSurfaceIntegration);
}

void BatchSamplerRendererQueue::LaunchNoWeightRayProcess() {
	BEGIN_TIMING(LaunchNoWeightRayProcess);
	for (int i = 0;i < taskNum;i++) {
		if (rayWeights[i] <= 0.f) {
			Ls[i] = 0.f;
			Ts[i] = 1.f;
		}
	}
	END_TIMING(LaunchNoWeightRayProcess);
}
void BatchSamplerRendererQueue::LaunchMissedRayProcess() {
	BEGIN_TIMING(LaunchMissedRayProcess);
	for (int i = 0;i < taskNum;i++) {
		if (rayWeights[i] > 0.f) {
			if (!hits[i]) {
				Spectrum& Li = Lis[i];
				Li = 0;
				for (uint32_t j = 0; j < render_task->scene->lights.size(); ++j) {
					Li += render_task->scene->lights[j]->Le(rays[i]);
				}
			}
		}
	}
	END_TIMING(LaunchMissedRayProcess);
}
void BatchSamplerRendererQueue::LaunchVolumeIntegration() {
	BEGIN_TIMING(LaunchVolumeIntegration);
	for (int i = 0;i < taskNum;i++) {
		if (rayWeights[i] > 0.f) {
			Lvis[i] = renderer->volumeIntegrator->Li(render_task->scene, renderer, rays[i], &samples[i], rng,
				&Ts[i], arenas[i]);
		}
	}
	END_TIMING(LaunchVolumeIntegration);
}

void BatchSamplerRendererQueue::LaunchCombineProcess() {
	BEGIN_TIMING(LaunchCombineProcess);
	for (int i = 0;i < taskNum;i++) {
		if (rayWeights[i] > 0.f) {
			Ls[i] = Ts[i] * Lis[i] + Lvis[i];
			// Issue warning if unexpected radiance value returned
			if (Ls[i].HasNaNs()) {
				Error("Not-a-number radiance value returned "
					"for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}
			else if (Ls[i].y() < -1e-5) {
				Error("Negative luminance value, %f, returned "
					"for image sample.  Setting to black.", Ls[i].y());
				Ls[i] = Spectrum(0.f);
			}
			else if (isinf(Ls[i].y())) {
				Error("Infinite luminance value returned "
					"for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}
		}
	}
	END_TIMING(LaunchCombineProcess);
}