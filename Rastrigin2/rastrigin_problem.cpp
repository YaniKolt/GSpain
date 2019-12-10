#include "rastrigin_problem.h"

#define _USE_MATH_DEFINES
#include <math.h>

// ------------------------------------------------------------------------------------------------
TRastriginProblem::TRastriginProblem()
{
  mIsInitialized = false;
  mDimension = 1;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::SetConfigPath(const std::string& configPath)
{
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::SetDimension(int dimension)
{
  if(dimension > 0 && dimension <= mMaxDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::Initialize()
{
  if (mDimension > 0)
  {
    mIsInitialized = true;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
void TRastriginProblem::GetBounds(double* lower, double *upper)
{
  if (mIsInitialized)
    for (int i = 0; i < mDimension; i++)
    {
      lower[i] = -1.8;
      upper[i] = 2.2;
    }
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::GetOptimumValue(double& value) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  value = 0.0;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::GetOptimumPoint(double* point) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  for (int i = 0; i < mDimension; i++)
    point[i] = 0.0;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int TRastriginProblem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double TRastriginProblem::CalculateFunctionals(const double* x, int fNumber)
{
  double sum = 0.;
  for (int j = 0; j < mDimension; j++)
    sum += x[j] * x[j] - 10. * cos(2.0 * M_PI * x[j]) + 10.0;
  return sum;
}

// ------------------------------------------------------------------------------------------------
TRastriginProblem::~TRastriginProblem()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new TRastriginProblem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------
