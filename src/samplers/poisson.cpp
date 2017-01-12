#include "stdafx.h"

#include "poisson.h"
#include <cmath>

#include <cstdio>
#include <stdexcept>
#include <string>


using std::vector;

#include <random>
#include <ctime>
class DefaultPRNG
{
public:
  DefaultPRNG()
    : m_Gen( std::random_device()() )
      , m_Dis( 0.0f, 1.0f )
  {
  // prepare PRNG
    m_Gen.seed( time( nullptr ) );
  }

  explicit DefaultPRNG( uint32_t seed )
  : m_Gen( seed )
    , m_Dis( 0.0f, 1.0f )
  {
  }

  float RandomFloat()
  {
    return static_cast<float>( m_Dis( m_Gen ) );
  }

  int RandomInt( int Max )
  {
    std::uniform_int_distribution<> DisInt( 0, Max );
    return DisInt( m_Gen );
  }

private:
  std::mt19937 m_Gen;
  std::uniform_real_distribution<float> m_Dis;
};

DefaultPRNG rng;


template<int DIM>
PoissonGridPoint<DIM>::PoissonGridPoint()
:m_valid(false)
{
}

template<int DIM>
PoissonGridPoint<DIM>::PoissonGridPoint(const float (&coor)[DIM])
:m_valid(true)
{
  for(int i=0;i<DIM;i++){
    coordinate[i] = coor[i];
  }
}

template<int DIM>
float PoissonGridPoint<DIM>::ComputeDistance(const PoissonGridPoint<DIM>& other) const{
  float result = 0.0;
  for(int i=0;i<DIM;i++){
    float diff = (*this)[i] - other[i];
    result += diff * diff;
  }
  return sqrt(result);
}
template<>
float PoissonGridPoint<1>::ComputeDistance(const PoissonGridPoint<1>& other) const{
  return fabs(other[0] - (*this)[0]);
}

template<int DIM>
bool PoissonGridPoint<DIM>::CheckCoordinateFit() const{
  for(int i=0;i<DIM;i++){
    if(coordinate[i] < 0.f || coordinate[i] > 1.f){
      return false;
    }
  }
  return true;
}

template<int DIM>
PoissonGenerator<DIM>::PoissonGenerator(int nSamples, float minDistance,  int k)
:m_nSamples(nSamples), m_k(k)
{
  if(minDistance < 0.f){
    m_minDistance = 1.f / powf(nSamples,1.f/DIM) / 1.25f;
  }
  else{
    m_minDistance = minDistance;
  }
  m_cellSize = m_minDistance / sqrtf(DIM);
  m_gridWidth = ceil(1.0f / m_cellSize);
 
  m_dimWidth[DIM - 1] = m_gridWidth;
  for(int i=DIM - 2;i>=0;i--){
    m_dimWidth[i] = m_dimWidth[i + 1] * m_gridWidth;
  }

#ifdef POISSON_VERBOSE
  printf("# dim: %d, min distance: %f, cell size: %f, grid width: %d\n", DIM, m_minDistance, m_cellSize, m_gridWidth);
#endif

  m_grid = new PoissonGridPoint<DIM>[m_dimWidth[0]];
}

template<int DIM>
PoissonGenerator<DIM>::~PoissonGenerator(){
  delete[] m_grid;
}

template<int DIM>
void PoissonGenerator<DIM>::GenerateRandomNumber(){
  for(int i=0;i<DIM;i++){
    m_numberBuffer[i] = rng.RandomFloat();
  }
}

template<int DIM>
PoissonGridPoint<DIM> PoissonGenerator<DIM>::GenerateRandomPoint(){
  GenerateRandomNumber();
  return PoissonGridPoint<DIM>(m_numberBuffer);
}

template<int DIM>
PoissonGridPoint<DIM> PoissonGenerator<DIM>::PopRandom(vector<PoissonGridPoint<DIM> >& activeList){
  int index = rng.RandomInt(activeList.size() - 1);  
  auto p = activeList[index];
  activeList.erase(activeList.begin() + index);
  return std::move(p);
}

template<int DIM>
bool PoissonGenerator<DIM>::HasNeighbor(const PoissonGridPoint<DIM>& p){
  // compute the corresponding grid
  int gridIndex[DIM];
  ComputeGridIndex(p, gridIndex);
  // start the loop
  int indice[DIM] = {gridIndex[0] - GRID_CHECK_SIZE};
  int indent = 0;
  for(;;){
    if(indice[indent] < (gridIndex[indent] + GRID_CHECK_SIZE)){
      if(indice[indent] >= 0 && indice[indent] < m_gridWidth){
        // run for current indent
        if(indent < DIM -1){
          // not the deepest
          // prepare for next indent
          ++indent;
          indice[indent] = gridIndex[indent] - GRID_CHECK_SIZE;
        }
        else{
          // deepest
          int index = ComputeAccumulatedGridIndex(indice); 
          auto& check_p = m_grid[index];
          if(check_p.m_valid){
            float distance = p.ComputeDistance(check_p);
            //printf("distance: %f\n",distance);
            if(distance < m_minDistance ){
              return true;
            }
          }
          ++indice[indent];
        }
      }
      else{
        // current index is invalid
        // just go to next iteration
        ++indice[indent];
      }
    }
    else{
      // current indent is finished
      --indent;
      if(indent < 0){
	// out of loop
	break;
      }
      else{
	++indice[indent];
      }
    }
  }
  return false;
}

template<int DIM>
void PoissonGenerator<DIM>::ComputeGridIndex(const PoissonGridPoint<DIM>& p, int (&result)[DIM]){
  int temp;
  for(int i=0;i<DIM;i++){
    temp = (int)(p[i]/m_cellSize);
    if(temp < 0){
      temp = 0;
    }
    if(temp >= m_gridWidth){
      temp = m_gridWidth - 1;
    }
    result[i] = temp;
  }
}

template<int DIM>
int PoissonGenerator<DIM>::ComputeAccumulatedGridIndex(const int (&indice)[DIM]) const{
  int result = indice[DIM - 1];
  for(int i=DIM - 2;i>=0;i--){
    result += indice[i] * m_dimWidth[i+1];
  }
  return result;
}



template<int DIM>
inline void PoissonGenerator<DIM>::InsertIntoGrid(PoissonGridPoint<DIM>&& p){
  int gridIndex[DIM];
  ComputeGridIndex(p, gridIndex);
  int index = ComputeAccumulatedGridIndex(gridIndex);
  m_grid[index] = p;
}

template<int DIM>
int PoissonGenerator<DIM>::PlaceSamples(float* samples, int offset, int step){
  if(m_nSamples <= 0){
    return 0;
  }
  // clear
  m_activeList.clear();
  for(int i=0;i<m_dimWidth[0];i++){
    m_grid[i].m_valid = false;
  }

  float* sampleIndex = samples + offset;

  int sampleCount = 1;
  // generate first point
  auto point = GenerateRandomPoint();
  // put that point in output and active
  m_activeList.push_back(point);
  for(int j=0;j<DIM;j++){
    sampleIndex[j] = point[j];
  }
  sampleIndex += step;
  // until active is empty
  do{
    // randomly pop one point
    auto pop_point = PopRandom(m_activeList);
    // generate k points around this point
    for(int i=0;i<m_k;i++){
      auto new_point = GenerateRandomPointAround(pop_point);
      // check the point with grid
      if(new_point.CheckCoordinateFit()){
        if(!HasNeighbor(new_point)){
          m_activeList.push_back(new_point);
          InsertIntoGrid(std::move(new_point));
          for(int j=0;j<DIM;j++){
            sampleIndex[j] = new_point[j];
          }
          sampleIndex += step;
          sampleCount++;
          if(sampleCount >= m_nSamples){
            break;
          }
        }
      }
    }
  }while(!m_activeList.empty() && sampleCount < m_nSamples);
#ifdef POISSON_VERBOSE
  printf("# Total %d samples\n", sampleCount);
#endif
  return sampleCount;
}

template<int DIM>
void PoissonGenerator<DIM>::GenerateCoordinateAround(const PoissonGridPoint<DIM>& p){
  char buf[1024];
  sprintf(buf, "GenerateCoordinateAround is not implemented in this dimension: %d\n", DIM);
  throw(std::runtime_error(buf));
}

template<int DIM>
PoissonGridPoint<DIM> PoissonGenerator<DIM>::GenerateRandomPointAround(const PoissonGridPoint<DIM>& p){
  GenerateCoordinateAround(p);
  return PoissonGridPoint<DIM>(m_numberBuffer);
}

// User provide: coordinate generation for 2D
template<>
void PoissonGenerator<2>::GenerateCoordinateAround(const PoissonGridPoint<2>& p){
  GenerateRandomNumber();
  // random radius
  float radius = m_minDistance * (m_numberBuffer[0] + 1.f);
  // random angle
  float angle = 2 * M_PI * m_numberBuffer[1];
  // save the coordinate
  m_numberBuffer[0] = p[0] + radius * cos(angle);
  m_numberBuffer[1] = p[1] + radius * sin(angle);
}
template class PoissonGenerator<2>;

// User provide: coordinate generation for 1D
template<>
void PoissonGenerator<1>::GenerateCoordinateAround(const PoissonGridPoint<1>& p){
  GenerateRandomNumber();
  // random radius
  if(m_numberBuffer[0] < 0.5f){
    m_numberBuffer[0] = p[0] - m_minDistance * (m_numberBuffer[0]*2.f + 1.f);
  }
  else{
    m_numberBuffer[0] = p[0] + m_minDistance * ((m_numberBuffer[0] - 0.5f)*2.f + 1.f);
  }
}
template class PoissonGenerator<1>;
