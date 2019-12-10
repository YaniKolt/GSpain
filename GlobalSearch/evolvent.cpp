
#include "evolvent.h"

#include <stdlib.h>
#include <string.h>
#include <omp.h>

//double Perebor(int maxTrial, int demension, IProblem* problem)
//{
//  int size;
//  int rank;
//  int delta;
//  int ost;
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  double min, sum;
//
//  double* lower_bounds = new double[demension];
//  double* upper_bounds = new double[demension];
//  problem->GetBounds(lower_bounds, upper_bounds);
//  double shag = lower_bounds[0];
//  double* masX = new double[maxTrial + 1];
//  masX[0] = shag;
//
//  //Подсчитаем координаты X у точек, в которых потом будем считать значение функции
//  for (int i = 1; i < maxTrial + 1; i++)
//  {
//    if (shag >= upper_bounds[0] + 0.1)
//      break;
//    shag += (upper_bounds[0] - lower_bounds[0]) / maxTrial;
//    masX[i] = shag;
//  }
//
//  int flag = 0;
//  int tmpsize = 0;
//  ost = (maxTrial + 1) % size;
//  if ((maxTrial + 1) <= size)
//  {
//    tmpsize = size;
//    size = 1;
//    ost = 0;
//    flag = -1;
//  }
//
//  delta = (maxTrial + 1) / size;
//
//  if (ost == 0)
//  {
//    if (rank == 0)
//    {
//      for (int i = delta; i < (maxTrial + 1) - delta; i += delta)
//      {
//        MPI_Send(&masX[i], delta, MPI_DOUBLE, i / delta, 0, MPI_COMM_WORLD);
//      }
//    }
//  }
//  else
//  {
//    if (rank == 0)
//    {
//      for (int i = delta + ost; i <= (maxTrial + 1) - delta; i += delta)
//      {
//        MPI_Send(&masX[i], delta, MPI_DOUBLE, (i - ost) / delta, 0, MPI_COMM_WORLD);
//      }
//    }
//  }
//
//  MPI_Status status;
//  double* b = new double[delta];
//  for (int i = 0; i < delta; i++)
//    b[i] = 0;
//  double prom_min;
//  double* tmp = new double[demension];
//  for (int i = 0; i < demension; i++)
//    tmp[i] = 0;
//  if (rank == 0)
//  {
//    for (int i = 0; i < delta + ost; i++)
//    {
//      for (int j = 0; j < demension; j++)
//        tmp[j] = masX[i];
//      sum = problem->CalculateFunctionals(tmp, 0);
//      //printf("Now Func %lf\n", sum);
//
//      if (i == 0)
//        prom_min = sum;
//      else
//        if (sum < prom_min)
//          prom_min = sum;
//    }
//  }
//  else
//  {
//    if (flag == 0)
//    {
//      MPI_Recv(&b[0], delta, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//      for (int i = 0; i < delta; i++)
//      {
//        for (int j = 0; j < demension; j++)
//          tmp[j] = b[i];
//        sum = problem->CalculateFunctionals(tmp, 0);
//        //printf("Now Func %lf\n", sum);
//
//        if (i == 0)
//          prom_min = sum;
//        else
//          if (sum < prom_min)
//            prom_min = sum;
//      }
//    }
//  }
//
//  if (flag == 0)
//  {
//    MPI_Reduce(&prom_min, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
//  }
//  else
//  {
//    if (rank == 0)
//    {
//      return prom_min;
//    }
//    else
//    {
//      return -1;
//    }
//  }
//  
//  MPI_Barrier(MPI_COMM_WORLD);
//  delete[] masX;
//  delete[] tmp;
//  delete[] lower_bounds;
//  delete[] upper_bounds;
//  delete[] b;
//  return min;
//}



void CalculateNode(double is, int n, int *u, int *v, int *l, double nexpExtended)
  // вычисление вспомогательного центра u(s) и соответствующих ему v(s) и l(s)
  // calculate u=u[s], v=v[s], l=l[s] by is=s
{
  int n1, i, j, k1, k2, iq;
  double iff;
  double nexp;

  iq = 1;
  for (nexp = 1, i = 0; i < n; nexp += nexp, i++);
  n1 = n - 1;
  *l = 0;
  if (is == 0)
  {
    *l = n1;
    for (i = 0; i < n; i++)
    {
      u[i] =- 1;
      v[i] =- 1;
    }
  }
  else if (is == (nexpExtended - 1.0))
  {
    *l = n1;
    u[0] = 1;
    v[0] = 1;
    for (i = 1; i < n; i++)
    {
      u[i] =- 1;
      v[i] =- 1;
    }
    v[n1] = 1;
  }
  else
  {
    iff = nexpExtended; 
    k1 =- 1;
    for (i = 0; i < n; i++)
    {
      iff = iff / 2;
      if (is >= iff)
      {
        if ((is == iff) && (is != 1.0))
        {
          *l = i;
          iq =- 1;
        }
        is -= iff;
        k2 = 1;
      }
      else
      {
        k2 =- 1;
        if ((is == (iff - 1.0)) && (is != 0))
        {
          *l = i;
          iq = 1;
        }
      }
      j =- k1 * k2;
      v[i] = j;
      u[i] = j;
      k1 = k2;
    }
    v[*l] = v[*l] * iq;
    v[n1] =- v[n1];
  }
}

// строим прообраз _y (принадлежащий гиперкубу по образу x (из отрезка 0..1)
// нижняя граница области поиска - A
// верхняя граница области поиска - B
// размерность - N
// плотность построения разверток - m
void GetImage(double x, double* _y, const double* A, const double* B, int N, int m)
{  // x ---> y
  if (N == 1)
  {
    _y[0] = x - 0.5;
  }
  else
  {
    double nexpExtended;
    int i = 0;
    for (nexpExtended = 1.0, i = 0; i < N; nexpExtended += nexpExtended, i++)
      ;
    int iu[MaxDim];
    int iv[MaxDim];
    int l;
    double d;
    int mn;
    double r;
    int iw[MaxDim];
    int it, j;
    double is;

    d = x;
    r = 0.5;
    it = 0;
    mn = m * N;
    for (i = 0; i < N; i++)
    {
      iw[i] = 1;
      _y[i] = 0.0;
    }
    for (j = 0; j < m; j++)
    {
      if (x == 1.0)
      {
        is = nexpExtended - 1.0;
        d = 0; // d = 0.0;
      }
      else
      {
        //Код из старой версии - уточнить работоспособность при N > 32
        d *= nexpExtended;
        is =(int)d;
        d  -= is;
      }
      CalculateNode(is, N, iu, iv, &l, nexpExtended);
      i = iu[0];
      iu[0] = iu[it];
      iu[it] = i;
      i = iv[0];
      iv[0] = iv[it];
      iv[it] = i;
      if (l == 0)
        l = it;
      else if (l == it)
        l = 0;
      r *= 0.5;
      it = l;
      for (i = 0; i < N; i++)
      {
        iu[i] *= iw[i];
        iw[i] *= -iv[i];
        _y[i] += r * iu[i];
      }
    }
    //return y; // it saves return value to y, so no need to call operator= again
  }
  //transform_P_to_D();
  int i;
  for (i = 0; i < N; i++)
    _y[i] = _y[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2;
}

void CalculateNumbr(double *s, int *u, int *v, int *l, double nexpExtended, int N)
  // calculate s(u)=is,l(u)=l,v(u)=iv by u=iu
{
  int i, k1, k2, l1;
  double is, iff;

  iff = nexpExtended;
  is = 0;
  k1 = -1;
  k2 = 0; 
  l1 = 0; 
  for (i = 0; i < N; i++)
  {
    iff = iff / 2;
    k2 =- k1 * u[i];
    v[i] = u[i];
    k1 = k2;
    if (k2 < 0)
      l1 = i;
    else
    {
      is += iff;
      *l = i;
    }
  }
  if (is == 0)
    *l = N - 1;
  else
  {
    v[N - 1] = -v[N - 1];
    if (is == (nexpExtended - 1))
      *l = N - 1;
    else
    {
      if (l1 == (N - 1))
        v[*l] = -v[*l];
      else
        *l = l1;
    }
  }
  *s = is;
}


// строим образ х (из отрезка 0..1) по прообразу  _y (принадлежащий гиперкубу)
// нижняя граница области поиска - A
// верхняя граница области поиска - B
// размерность - N
// плотность построения разверток - m
void GetInverseImage( double* _y, double& x, const double* A, const double* B, int N, int m)
{// y ---> x

  double* y = new double [N];
  memcpy(y, _y, N * sizeof(double));

  for (int i = 0; i < N; i++)
    y[i] = (y[i] - (A[i] + B[i]) / 2) / (B[i] - A[i]);


  double nexpExtended;
  int i = 0;
  for (nexpExtended = 1.0, i = 0; i < N; nexpExtended += nexpExtended, i++)
      ;
  int u[MaxDim], v[MaxDim];
  double r1;
  if (N == 1)
  {
    x = y[0] + 0.5;
    x;
  }

  double  r;
  int w[MaxDim];
  int j, it, l;
  double is;

  for (i = 0; i < N; i++)
    w[i] = 1;
  r = 0.5;
  r1 = 1;
  x = 0;
  it = 0;
  for (j = 0; j < m; j++)
  {
    r *= 0.5;
    for (i = 0; i < N; i++)
    {
      u[i] = (y[i] < 0) ? -1 : 1;
      y[i] -= r * u[i];
      u[i] *= w[i];
    }
    i = u[0];
    u[0] = u[it];
    u[it] = i;
    CalculateNumbr(&is, u, v, &l, nexpExtended, N);
    i = v[0];
    v[0] = v[it];
    v[it] = i;
    for (i = 0; i < N; i++)
      w[i] *= -v[i];
    if (l == 0)
      l = it;
    else
      if (l == it)
        l = 0;
    it = l;
    r1 = r1 / nexpExtended;
    x += r1 * is;
  } 

}