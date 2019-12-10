#ifndef __GRISHAGINPROBLEM_H__
#define __GRISHAGINPROBLEM_H__

#include "../GlobalSearch/problem_interface.h"
#include "grishagin_function.h"

class TGrishaginProblem : public IProblem
{
protected:

  int mDimension;
  bool mIsInitialized;
  static const int mSupportedDimension = 2;
  vagrisfunction::GrishaginFunction mFunction;
  std::string mConfigPath;

public:

  TGrishaginProblem();

  virtual int SetConfigPath(const std::string& configPath);
  virtual int SetDimension(int dimension);
  virtual int GetDimension() const;
  virtual int Initialize();

  virtual void GetBounds(double* upper, double *lower);
  virtual int GetOptimumValue(double& value) const;
  virtual int GetOptimumPoint(double* x) const;

  virtual int GetNumberOfFunctions() const;
  virtual int GetNumberOfConstraints() const;
  virtual int GetNumberOfCriterions() const;

  virtual double CalculateFunctionals(const double* x, int fNumber);

  ~TGrishaginProblem();
};

extern "C" LIB_EXPORT_API IProblem* create();
extern "C" LIB_EXPORT_API void destroy(IProblem* ptr);

#endif
// - end of file ----------------------------------------------------------------------------------
