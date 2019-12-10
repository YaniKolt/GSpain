#include "evolvent.h"
#include "problem_manager.h"
#include "extended.h"
#include "math.h"
// #include "mpi.h"


#include <stdio.h>
#include <omp.h>

int minimI(double* m, int size) {
  double min = m[0]; 
  int in = 0;
    for (int i = 1; i < size; i++){
      if (m[i] < min) {
        min = m[i];
        in = i;
      }
    }
  return in;
}



//double* Perebor(int numOfPaint, int demension, IProblem* problem)
//{
//  int size;
//    int rank;
//    int segm;
//    int ost;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    segm = numOfPaint / size;
//    ost = numOfPaint % size;
//     
//  double* low_bo = new double[demension]; 
//  double* up_bo = new double[demension];
//  problem->GetBounds(low_bo, up_bo); //получаем нижнюю и верхнюю границу
//  double* delta = new double[demension];
//  for (int i = 0; i < demension; i++) {
//    delta[i] = ((up_bo[i] - low_bo[i]) / (numOfPaint - 1)); //ищем, с каким шагом мы идем по каждой координате
//  }
//  
//  if (rank == 0) {
//    double** mas = new double*[numOfPaint];
//    for (int i = 0; i < numOfPaint; i++)
//      mas[i] = new double[demension];
//    double* res = new double[numOfPaint];
//  }
//  double** lockmas = new double*[segm];
//  for (int i = 0; i < segm; i++)
//    lockmas[i] = new double[demension];
//
// for (int j = 0; j < demension; j++) { 
//    for (int i = 0; i < segm; i++) //ищем значения точек в одном процессе
//    {
//      lockmas[i][j] = low_bo[j] + delta[j] * (segm*rank+i); // если 3 процесса и 7 точек: 00 11 22 0
//    }
//    if (rank == 0) {
//      for (int i = 0; i < ost; i++) {
//        mas[(segm*size+1)][j]  = low_bo[j] + delta[j] * (segm*size + i);
//      }
//    }
//  }
//
// MPI_Gather(&lockmas, segm, MPI_DOUBLE, &mas, segm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
// double* lockres = new double[segm];
//  for (int i = 0; i < segm; i++) {
//    lockres[i] = problem->CalculateFunctionals(lockmas[i], 0);
//  }
//  if (rank == 0) {
//    for (int i = 0; i < ost; i++) {
//      res[segm*size +i] = problem->CalculateFunctionals(mas[segm*size+i], 0);
//    }
//  }
//
//  MPI_Gather(&lockres, segm, MPI_DOUBLE, &res, segm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//  if (rank == 0) {
//    int optim;
//    optim = minimI(res, numOfPaint);
//    return  mas[optim];
//  }
//}


void sort(double* &arr, double* &arr2, int size) {
  double temp; 
  double temp2;
  for (int i = 0; i < size - 1; i++) {
    for (int j = 0; j < size - i - 1; j++) {
      if (arr[j] > arr[j + 1]) {
        temp = arr[j];
        arr[j] = arr[j + 1];
        arr[j + 1] = temp;
        temp2 = arr2[j];
        arr2[j] = arr2[j + 1];
        arr2[j + 1] = temp2;
      }
    }
  }
}

int maxIn(double* m, int size) {
  double max = m[1];
  int in = 1;
  for (int i = 2; i < size; i++) {
    if (m[i] > max) {
      max = m[i];
      in = i;
    }
  }
  return in;
}

double maxZn(double* m, int size) {
  double max = m[1];
  for (int i = 2; i < size; i++) {
    if (m[i] > max) {
      max = m[i];
    }
  }
  return max;
}

int binSe(double* arr, int size, double key) {
  bool flag = false;
  int l = 0;
  int r = size; 
  int mid;

  while ((l <= r) && (flag != true)) {
    mid = (l + r) / 2;
    if (arr[mid] == key) flag = true; 
    if (arr[mid] > key) r = mid - 1;
    else l = mid + 1;
  }
  return mid;
}


double GlobSe(IProblem* problem, double e, int demension, int &iter,  int max) {
  double* low_bo = new double[demension];
  double* up_bo = new double[demension];
  double* tmpX = new double[demension];
  double r = 4.0;
  for (int i = 0; i < demension; i++) {
    tmpX[i] = 0;
  }
  problem->GetBounds(low_bo, up_bo);
 // double* masX = new double[(up_bo[0]-low_bo[0])/r]; //максимальное количество иттераций 
  //double* masX = new double[(up_bo[0] - low_bo[0]) / r];
  double* masX = new double[max];
  //double* masY = new double[(up_bo[0] - low_bo[0]) / r]; // храним значения предыдущих точек, чтобы не пересчитывать
  double* masY = new double[max];
  //double* muarr = new double[(up_bo[0] - low_bo[0]) / r];
  double* muarr = new double[max];
  muarr[0] = 0;
 // double* R = new double[(up_bo[0] - low_bo[0]) / r];
  double* R = new double[max];
  R[0] = 0;
  masX[0] = low_bo[0];
  masX[1] = up_bo[0];
  tmpX[0] = masX[0];
  masY[0] = problem->CalculateFunctionals(tmpX, 0);
  tmpX[0] = masX[1];
  masY[1] = problem->CalculateFunctionals(tmpX, 0);

  //problem->CalculateFunctionals(mas[segm*size + i], 0);
  int kk=2;
  int k = 2;
  double mu;
  int t;
  int indK;
  int signal = 0;
  double xk;
  int minkkk;
  if (masY[0] < masY[1])
    minkkk = 0;
  else
    minkkk = 1;
  double minY = masY[minkkk];
  double minX = masX[minkkk];
  while ((signal != 1) && (k < max+2)) {
   // for (k = 2; k < max + 2; k++) {
      //нужно сравнивать значения пришедшего по У со значением текущего минимума по У. Если от нового по У меньше, то приравниваем минимуму его. 
      // нужно сравнить по е текущий минимум. 
      sort(masX, masY, k);
      for (int i = 1; i < k; i++) {
        muarr[i] = ((abs(masY[i] - masY[i - 1])) / (masX[i] - masX[i - 1]));
      }
      mu = maxZn(muarr, k);
      for (int i = 1; i < k; i++) {
        R[i] = (masX[i] - masX[i - 1]) + (((masY[i] - masY[i - 1])*(masY[i] - masY[i - 1])) / (r*r*mu*mu*(masX[i] - masX[i - 1]))) - 2 * (masY[i] + masY[i - 1]) / (r*mu);
      }
      //R[0] = 2 * (masX[1] - masX[0]) - 4 * (masY[0]) / (r*mu);///Вот вообще не понятно, что вместо М на слайде. Пусть будет R*mu
      t = maxIn(R, k);
      masX[k] =( (masX[t] + masX[t - 1]) / 2 )- ((masY[t] - masY[t - 1]) / (2 * r*mu));
      tmpX[0] = masX[k];
      masY[k] = problem->CalculateFunctionals(tmpX, 0);

      if (masY[k] < minY) {
        if ((abs(masX[t] - masX[k]) < e) || (abs(masX[t-1] - masX[k]) < e)) {
          minX = masX[k];
          minY = masY[k];
          //printf("min  %d %lf %lf, \n  ", k, minX, minY);
          signal = 1;
        }
        else {
          minX = masX[k];
          minY = masY[k];
          //printf("%d %lf %lf, \n  ", k, minX, minY);


        }
      }
     // printf("%d %lf %lf, \n  ", k, masX[k], masY[k]);
      k++;
      kk++;
    //}
  }

  //minX = masX[k];
 // minY = masY[k];
  printf("Kol itter %d,  ", kk-2);
  iter = kk-2;
  return minX;
}



int  main(int argc, char* argv[])
{
  /// мэнаджер подключения задач
  TProblemManager manager;
  /// Число проведенных испытаний
  int trialCount = 0;
  /// Найдена ли точка глобального минимума
  bool isFindOptimalPoint = false;
  /// 
  //int iter = 0;
  /// Координаты точек испытания
  double** points = 0;
  /// Максимальное число испытаний
  int maxTrial = 10000;
  if (argc > 1)
  {
    /// Имя задачи
    std::string problemName = argv[1];
    /// Конфигурационный файл
    std::string problepConf = "";
    if (argc > 2)
    {
      problepConf = argv[2];
    }

    int err;
    /// Загружаем задачу
    if (manager.LoadProblemLibrary(problemName) != TProblemManager::OK_)
      return 0;
    /// Получаем задачу
    IProblem* problem = manager.GetProblem();
    /// Устанавливаем размерность
    err = problem->SetDimension(1);
    if (err != TProblemManager::OK_)
    {
      printf("Error SetDimension!\n");
      return 0;
    }
    /// Задаем файл с конфигурацией задачи
    err = problem->SetConfigPath(problepConf);
    if (err != TProblemManager::OK_)
    {
      printf("Error SetConfigPath!\n");
      return 0;
    }
    /// Инициализация
    err = problem->Initialize();
    if (err != TProblemManager::OK_)
    {
      printf("Error Initialize!\n");
      return 0;
    }
    /// Получаем размерность из задачи
    int dimension = problem->GetDimension();
    /// Граници области поиска
    double* lower_bounds = new double[dimension];
    double* upper_bounds = new double[dimension];
    /// Получаем Граници области
    problem->GetBounds(lower_bounds, upper_bounds);

    /// Координаты найденного  минимума
    double* minx = new double[dimension];
    /// Значение минимума
    double minf = -std::numeric_limits<double>::infinity();

    if (argc > 3)
    {
      points = new double*[maxTrial * 2];
    }



    //if (Extended::GetPrecision() == 0.01)
    //{
    //  if (parameters.m * (pTask->GetFreeN() - pTask->GetNumberOfDiscreteVariable()) <= 50)
    //  {
    //    Extended::SetTypeID(etDouble);
    //  }
    //  Extended::SetPrecision(1 / (::pow(2., parameters.m * (pTask->GetFreeN() -
    //    pTask->GetNumberOfDiscreteVariable()))));
    //}


    /////////////////////
    ///Вычисоения
    ////////////////////

    minx[0] = 0.5;
    minx[1] = 0;

    //    int n = 100;
    //
    //    Extended* ex = new Extended [n];
    //    Extended* ey = new Extended[n];
    //    Extended* ez = new Extended[n];
    //
    //    for (int i = 0; i < n; i++)
    //    {
    //      ex[i] = 1.0 / (i + 1.0);
    //      ey[i] = i * 100.0;
    //      ez[i] = (ex[i] + ey[i]) / 1000000.0;
    //
    //      //printf("ez[%d] = %lf\n", i, ez[i].toDouble());
    //    }
    //
    //    int t = 4;
    //#pragma omp parallel for num_threads(t)
    //    for (int i = 0; i < t; i++)
    //    {
    //      int thread = omp_get_thread_num();
    //      printf("%d\n", thread);
    //    }
    //
    //    Extended* exp = new Extended[n];
    //    Extended* eyp = new Extended[n];
    //    Extended* ezp = new Extended[n];
    //
    //#pragma omp parallel for num_threads(t)
    //    for (int i = 0; i < n; i++)
    //    {
    //      exp[i] = 1.0 / (i + 1.0);
    //      eyp[i] = i * 100.0;
    //      ezp[i] = (exp[i] + eyp[i]) / 1000000.0;
    //
    //      //printf("ez[%d] = %lf\n", i, ez[i].toDouble());
    //    }
    //
    //    for (int i = 0; i < n; i++)
    //    {
    //      if ((ex[i] != exp[i]) || (ey[i] != eyp[i]) || (ez[i] != ezp[i]))
    //        printf("ez[%d] = %lf\tezp[%d] = %lf\n", i, ez[i].toDouble(), i, ezp[i].toDouble());
    //    }
    //
    //
    minf = problem->CalculateFunctionals(minx, 0);

    double x = 0.5;
    GetInverseImage(minx, x, lower_bounds, upper_bounds, dimension, 10);
    GetImage(x, minx, lower_bounds, upper_bounds, dimension, 10);

    trialCount++;

    if (argc > 3)
    {
      points[0] = new double[2];

      points[0][0] = minx[0];
      points[0][1] = minx[1];
    }

    ////////////////////
    ///
    ///////////////////

    for (int i = 0; i < dimension; i++)
      printf("x[%d] = %lf\n", i, minx[i]);
    printf("min = %lf\n", minf);
    printf("Point count = %d\n\n", trialCount);

    /// Известная точка глобального минимума
    double* BestTrialy = new double[dimension];
    if (problem->GetOptimumPoint(BestTrialy) == TProblemManager::OK_)
    {
      for (int j = 0; j < dimension; j++)
        printf("Optimal x[%d] = %lf \n", j, BestTrialy[j]);
    }

    // Проверяем попадание в окрестность глобального оптимума
    double Epsilon = 0.01;
    isFindOptimalPoint = true;
    for (int i = 0; i < dimension; i++)
    {
      double fabsx = fabs(BestTrialy[i] - minx[i]);
      double fm = Epsilon * (upper_bounds[i] - lower_bounds[i]);
      if (fabsx > fm)
      {
        isFindOptimalPoint = false;
        break;
      }
    }

    if (isFindOptimalPoint)
    {
      printf("\nGlobal optimum FOUND!\n\n");
    }

    if (argc > 3)
    {
      /// Имя файла с координатами точек
      std::string pointLogName = argv[3];
      FILE* pointLog = fopen(pointLogName.c_str(), "w");
      /// Печатаем точки испытания
      fprintf(pointLog, "%d\n", trialCount);
      for (int i = 0; i < trialCount; i++)
      {
        for (int j = 0; j < dimension; j++)
          fprintf(pointLog, "%lf ", points[i][j]);
        fprintf(pointLog, "\n");
      }

      /// Печатаем найденную точку
      for (int j = 0; j < dimension; j++)
        fprintf(pointLog, "%lf ", minx[j]);
      fprintf(pointLog, "\n");

      /// Известная точка глобального минимума
      if (problem->GetOptimumPoint(BestTrialy) == TProblemManager::OK_)
      {
        for (int j = 0; j < dimension - 1; j++)
          fprintf(pointLog, "%lf ", BestTrialy[j]);
        fprintf(pointLog, "%lf", BestTrialy[dimension - 1]);
      }

      fclose(pointLog);
      delete[] points;
    }

    int dem = 1;
    int iter = 106;
    int max = 100000;
    double minGS = 10;
    minGS = GlobSe(problem, 0.000001, dem, iter, max);
    printf("AGP %lf, \n  ", minGS); 
   // printf("%lf,  ", iter);

    //int dem = 1;
    //double* pointis = new double[dem];
    //pointis = Perebor(2, dem, problem);
    //double* low = new double[dem];
    //double* up = new double[dem];
    //problem->GetBounds(low, up); //получаем нижнюю и верхнюю границу
    //for (int i = 0; i < dem; i++) {
    //  printf("%lf,  ", up[i]);
    //}
    //printf("%lf", pointis[0]); //как сделать вывод массива?
    //printf("\n\n ----- Result Perebor: ----- \n");
    //for (int i = 0; i < dem; i++) {
    //  printf("%lf,  ", pointis[i]);
    //}
  }
  return 0;
}