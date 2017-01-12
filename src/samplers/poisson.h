#ifndef _POISSON_H_INCLUDED
#define _POISSON_H_INCLUDED

#include <vector>

#define POISSON_VERBOSE

#define PRINT_POINT(point, PREFIX) { \
  printf(PREFIX);\
  for(int i=0;i<DIM;i++){\
    printf("%f, ", point[i]);\
  }\
  printf("\n");\
}

template<int DIM>
class PoissonGridPoint{
public:
  PoissonGridPoint();
  PoissonGridPoint(const float (&coor)[DIM]);
  float coordinate[DIM]; 
  bool m_valid;

  bool CheckCoordinateFit() const;
  inline float ComputeDistance(const PoissonGridPoint<DIM>& other) const;

  inline float& operator[](int index){ return coordinate[index]; }
  inline float operator[](int index) const { return coordinate[index]; }
};

template<int DIM>
class PoissonGenerator{
public:
  PoissonGenerator(int nSamples, float minDistance = -1 ,int k = 30);
  ~PoissonGenerator();
  int PlaceSamples(float* samples, int offset = 0, int step = DIM);
  
  static const int GRID_CHECK_SIZE = 3;


private:

  inline PoissonGridPoint<DIM> PopRandom(std::vector<PoissonGridPoint<DIM> >& activeList);
  inline PoissonGridPoint<DIM> GenerateRandomPoint();
  inline void GenerateRandomNumber();
  inline void ComputeGridIndex(const PoissonGridPoint<DIM>& p, int (&result)[DIM]);
  inline int ComputeAccumulatedGridIndex(const int (&indice)[DIM]) const;
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

  PoissonGridPoint<DIM>* m_grid;
};


#endif
