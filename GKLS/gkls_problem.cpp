#include "gkls_problem.h"
#include "pugixml.hpp"

#include <math.h>
#include <string>

// ------------------------------------------------------------------------------------------------
TGKLSProblem::TGKLSProblem()
{
  mIsInitialized = false;
  mDimension = 0;
  mPFunction = NULL;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::SetConfigPath(const std::string& configPath)
{
  mConfigPath = std::string(configPath);
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::SetDimension(int dimension)
{
  if(dimension > 1 && dimension <= mMaxDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::Initialize()
{
  if (!mIsInitialized)
  {
    int funcNumber = 0;
    int dimension = 0;
    int numMinuma = 0;
    double global_dist = 0.0;
    double global_radius = 0.0;

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(mConfigPath.c_str());
    if (result.status != pugi::xml_parse_status::status_ok)
      return IProblem::ERROR;

    pugi::xml_node config = doc.child("config");
    try
    {
      funcNumber = std::stoi(config.child("function_number").child_value());
      dimension = std::stoi(config.child("dimension").child_value());
      numMinuma = std::stoi(config.child("num_minima").child_value());
      global_dist = std::stod(config.child("global_dist").child_value());
      global_radius = std::stod(config.child("global_radius").child_value());
      if (funcNumber < 1 || funcNumber > 100)
        throw std::invalid_argument("");
      if (dimension < 1 || dimension > mMaxDimension)
        throw std::invalid_argument("");
      mDimension = dimension;
    }
    catch (std::invalid_argument except)
    {
      return IProblem::ERROR;
    }

    mPFunction = new gklsfunction::GKLSFunction();
    mPFunction->SetGlobalDistance(global_dist);
    mPFunction->SetGlobalRadius(global_radius);
    mPFunction->SetGlobalMinimumValue(-1.0);
    mPFunction->SetNumberOfLocalMinima(numMinuma);
    mPFunction->SetDimention(mDimension);

    if (mPFunction->CheckParameters() != GKLS_OK)
      return IProblem::ERROR;

    mPFunction->SetFunctionNumber(funcNumber);
    mIsInitialized = true;

    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
void TGKLSProblem::GetBounds(double* lower, double *upper)
{
  if (mIsInitialized)
  {
    for (unsigned i = 0; i < mPFunction->GetDimention(); i++)
    {
      lower[i] = -1.0;
      upper[i] =  1.0;
    }
  }
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::GetOptimumValue(double& value) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  value = mPFunction->GetGlobalMinimumValue();

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::GetOptimumPoint(double* point) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  mPFunction->GetGlobalMinimumPoint(point);

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateFunctionals(const double* y, int fNumber)
{
  return mPFunction->EvaluateDFunction(y);
}

// ------------------------------------------------------------------------------------------------
TGKLSProblem::~TGKLSProblem()
{
  if (mIsInitialized)
    delete mPFunction;
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new TGKLSProblem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------
