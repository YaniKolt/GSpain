#ifndef __GKLSPROBLEM_H__
#define __GKLSPROBLEM_H__

#include "../GlobalSearch/problem_interface.h"
#include "gkls_function.h"

class TGKLSProblem : public IProblem
{
protected:

  int mDimension;
  bool mIsInitialized;
  static const int mMaxDimension = 50;
  gklsfunction::GKLSFunction *mPFunction;
  std::string mConfigPath;

public:

  TGKLSProblem();

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

  ~TGKLSProblem();
};

extern "C" LIB_EXPORT_API IProblem* create();
extern "C" LIB_EXPORT_API void destroy(IProblem* ptr);

#endif
// - end of file ----------------------------------------------------------------------------------
