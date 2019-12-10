#include "grishagin_problem.h"
#include "pugixml.hpp"

#include <math.h>
#include <string>

// ------------------------------------------------------------------------------------------------
TGrishaginProblem::TGrishaginProblem()
{
  mIsInitialized = false;
  mDimension = mSupportedDimension;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::SetConfigPath(const std::string& configPath)
{
  mConfigPath = std::string(configPath);
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::SetDimension(int dimension)
{
  if(dimension == mSupportedDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::Initialize()
{
  if (mDimension == mSupportedDimension && !mIsInitialized)
  {
    int funcNumber = 0;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(mConfigPath.c_str());
    if (result.status != pugi::status_ok)
      return IProblem::ERROR;

    pugi::xml_node config = doc.child("config");
    if (config.child("function_number"))
    {
      try
      {
        funcNumber = std::stoi(config.child("function_number").child_value());
        if (funcNumber < 1 || funcNumber > 100)
          throw std::invalid_argument("");
      }
      catch (std::invalid_argument except)
      {
        return IProblem::ERROR;
      }
    }
    else
      return IProblem::ERROR;

    mFunction = vagrisfunction::GrishaginFunction();
    mFunction.SetFunctionNumber(funcNumber);
    mIsInitialized = true;

    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
void TGrishaginProblem::GetBounds(double* lower, double *upper)
{
  if (mIsInitialized)
  {
    lower[0] = lower[1] = 0.0;
    upper[0] = upper[1] = 1.0;
  }
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::GetOptimumValue(double& value) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  double optPoint[2];
  mFunction.GetOptimumCoordinates(optPoint);
  value = mFunction.EvaluateFunction(optPoint);

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::GetOptimumPoint(double* point) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  mFunction.GetOptimumCoordinates(point);

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int TGrishaginProblem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double TGrishaginProblem::CalculateFunctionals(const double* y, int fNumber)
{
  return mFunction.EvaluateFunction(y);
}

// ------------------------------------------------------------------------------------------------
TGrishaginProblem::~TGrishaginProblem()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new TGrishaginProblem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------
