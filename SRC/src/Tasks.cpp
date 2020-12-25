#include "Tasks.h"
#include <iostream>
#include "pugixml.hpp"
#include "../sample_src/grishagin_function.hpp"
#include "../sample_src/ShekelProblem.hpp"
#include "../sample_src/HillProblem.hpp"
//#include "../GlobalSearch/AGPHeap.h"
//#include "../build/SRC/AHeap.h"


// ------------------------------------------------------------------------------------------------
MyProblem::MyProblem()
{
}

// ------------------------------------------------------------------------------------------------
int MyProblem::SetConfigPath(const std::string& configPath)
{
  mConfigPath = std::string(configPath);
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::SetDimension(int dimension)
{
  /*if (dimension > 0 && dimension <= mMaxDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;*/
  dim = dimension;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetDimension() const
{
  return BigProblem->GetDimension();
}

// ------------------------------------------------------------------------------------------------
int MyProblem::Initialize()
{
  /*BigProblem = new TGrishaginProblem(10);
  std::cout << "Grishagin Problem 1" << std::endl;
  std::cout << "GrishaginProblem ( 0.5, 0.5 ) = " << BigProblem->ComputeFunction({ 0.5, 0.5 }) << std::endl;
  std::cout << "GrishaginProblem Derivatives ( 0.5, 0.5 ) = {" <<
    BigProblem->ComputeFunctionDerivatives({ 0.5, 0.5 })[0] << ", " <<
    BigProblem->ComputeFunctionDerivatives({ 0.5, 0.5 })[1] << "}" << std::endl;
  std::cout << std::endl;
  return IProblem::OK;*/

  ///Копипаст из гришагинпроблем
  //if (mDimension == mSupportedDimension && !mIsInitialized)
  //{

  



    int funcNumber = 0;

   /* BigProblem = new TShekelProblem(funcNumber);
    mIsInitialized = true;

    return IProblem::OK;*/


    int funcClass = 0;
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
    //pugi::xml_node config = doc.child("config");
    if (config.child("function_class"))
    {
      try
      {
        funcClass = std::stoi(config.child("function_class").child_value());
        if (funcClass < 1 || funcClass > 10)
          throw std::invalid_argument("");
      }
      catch (std::invalid_argument except)
      {
        return IProblem::ERROR;
      }
    }
    else
      return IProblem::ERROR;
    if (funcClass == 1) {
      BigProblem = new TShekelProblem(funcNumber);
      mIsInitialized = true;

    }
    if (funcClass == 2) {
      BigProblem = new THillProblem(funcNumber);
      mIsInitialized = true;

    }


    ///Раскоментить
 /*   mFunction = vagrisfunction::GrishaginFunction();
    mFunction.SetFunctionNumber(funcNumber);
    mIsInitialized = true;*/

    return IProblem::OK;
  //}
 // else
    //return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
void MyProblem::GetBounds(double* lower, double *upper)
{
  std::vector<double> l(dim);
  std::vector<double> u(dim);
  BigProblem->GetBounds(l, u);
  for (int i = 0; i < dim; i++) {
    lower[i] = l[i];
    upper[i] = u[i];
  }
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetOptimumValue(double& value) const
{
  value = BigProblem->GetOptimumValue();
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetOptimumPoint(double* point) const
{
  std::vector<double> x = BigProblem->GetOptimumPoint();
  for (int i = 0; i < dim; i++)
    point[i] = x[i];

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double MyProblem::CalculateFunctionals(const double* x, int fNumber)
{
  std::vector<double> f;
  for (int i = 0; i < dim; i++)
    f.push_back(x[i]);

  return BigProblem->ComputeFunction(f);
}

// ------------------------------------------------------------------------------------------------
MyProblem::~MyProblem()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new MyProblem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------