#include "be_function.h"

#define _USE_MATH_DEFINES
#include <math.h>

// ------------------------------------------------------------------------------------------------
TBeFun::TBeFun()
{
  mIsInitialized = false;
  mDimension = 1;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::SetConfigPath(const std::string& configPath)
{
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::SetDimension(int dimension)
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
int TBeFun::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::Initialize()
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
void TBeFun::GetBounds(double* lower, double *upper)
{
  if (mIsInitialized)
    for (int i = 0; i < mDimension; i++)
    {
      lower[i] = -3; //-1.8;
      upper[i] = 2; //2.2;
    }
}

// ------------------------------------------------------------------------------------------------
int TBeFun::GetOptimumValue(double& value) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  value = 0.0;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::GetOptimumPoint(double* point) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  for (int i = 0; i < mDimension; i++)
    point[i] = 0.0;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int TBeFun::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double TBeFun::CalculateFunctionals(const double* x, int fNumber)
{
  double sum = 0.;
  for (int j = 0; j < mDimension; j++)
    //sum += x[j] * x[j] - 10. * cos(2.0 * M_PI * x[j]) + 10.0;
    sum += x[j] * x[j];
  return sum;
}

// ------------------------------------------------------------------------------------------------
TBeFun::~TBeFun()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new TBeFun();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------
