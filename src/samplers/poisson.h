#if defined(_MSC_VER)
#pragma once
#endif


#ifndef _POISSON_H_INCLUDED
#define _POISSON_H_INCLUDED


#include "sampler.h"
#include "paramset.h"
#include "film.h"

#include <vector>

#define PRINT_POINT(point, PREFIX) { \
  printf(PREFIX);\
  for(int i=0;i<DIM;i++){\
    printf("%f, ", point[i]);\
  }\
  printf("\n");\
}

#define GENERATE_FROM_IMAGE 1
#define GENERATE_FROM_SAMPLE 2
#define GENERATE_FROM_RANDOM 3

#define CAMERA_SAMPLE_GENERATE GENERATE_FROM_IMAGE
#define TIME_SAMPLE_GENERATE GENERATE_FROM_RANDOM

template<int DIM>
class PoissonGridPoint {
public:
	PoissonGridPoint();
	PoissonGridPoint(const float(&coor)[DIM]);
	float coordinate[DIM];
	bool m_valid;

	inline float ComputeDistance(const PoissonGridPoint<DIM>& other) const;

	inline float& operator[](int index) { return coordinate[index]; }
	inline float operator[](int index) const { return coordinate[index]; }
};

template<int DIM>
class PoissonGenerator {
public:
	PoissonGenerator(int nSamples, float* rangeLimit = nullptr, float minDistance = -1, int k = 30);
	~PoissonGenerator();
	int PlaceSamples(float* samples, int offset = 0, int step = DIM);

	

	void SetPRNG(RNG* pRng);

	bool CheckCoordinateFit(const PoissonGridPoint<DIM>& p) const;

	static const int GRID_CHECK_SIZE = 3;
	float GetMinDistanceAdjustValue() const;

private:

	inline PoissonGridPoint<DIM> PopRandom(std::vector<PoissonGridPoint<DIM> >& activeList);
	inline PoissonGridPoint<DIM> GenerateRandomPoint();
	inline void GenerateRandomNumber();
	inline void ComputeGridIndex(const PoissonGridPoint<DIM>& p, int(&result)[DIM]);
	inline int ComputeAccumulatedGridIndex(const int(&indice)[DIM]) const;
	bool HasNeighbor(const PoissonGridPoint<DIM>& p);
	void GenerateCoordinateAround(const PoissonGridPoint<DIM>& p);
	inline PoissonGridPoint<DIM> GenerateRandomPointAround(const PoissonGridPoint<DIM>& p);
	inline void InsertIntoGrid(PoissonGridPoint<DIM>&& p);

	int m_k;
	int m_nSamples;

	float m_minDistance;
	float m_cellSize;
	int m_gridWidth;

	int m_dimWidth[DIM];

	float m_numberBuffer[DIM];
	std::vector<PoissonGridPoint<DIM> > m_activeList;

	float m_rangeLimit[DIM];

	RNG localRng;
	RNG* pRng;


	PoissonGridPoint<DIM>* m_grid;
};



// BestCandidateSampler Declarations
class PoissonDiskSampler : public Sampler {
public:

	enum SampleMode { mode_repeat, mode_single };

	// PoissonDiskSampler Public Methods
	PoissonDiskSampler(int xstart, int xend, int ystart, int yend,
		int nPixelSamples, float sopen, float sclose, int nMaxSample, SampleMode mode = mode_repeat);
	~PoissonDiskSampler();
	Sampler *GetSubSampler(int num, int count);
	int RoundSize(int size) const {
		return RoundUpPow2(size);
	}
	int MaximumSampleCount() { return 1; }
	int GetMoreSamples(Sample *sample, RNG &rng);
private:
	// PoissonDiskSampler Private Data

	void PrepareNewSamples(RNG& rng);

	SampleMode m_SampleMode;

	int nMaxSample;

	RNG rng;

	//int table.

	int nTotalSamplesRequired;
	int nTotalSamplesPrepared;
	int nValidSamples;
	int nCurrentSampleIndex;

	int nEmittedSamples;

	float* samples;

	PoissonGenerator<1>* pGenerator_1D;
	PoissonGenerator<2>* pGenerator_2D;

	float xTileWitdh, yTileWitdh;

};


PoissonDiskSampler *CreatePoissonDiskSampler(const ParamSet &params, const Film *film,
	const Camera *camera);



#endif
