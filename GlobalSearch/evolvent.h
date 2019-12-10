#ifndef EVOLVENT_H
#define  EVOLVENT_H

#define MaxDim 50

// строим прообраз _y (принадлежащий гиперкубу по образу x (из отрезка 0..1)
// нижняя граница области поиска - A
// верхняя граница области поиска - B
// размерность - N
// плотность построения разверток - m
void GetImage(double x, double* _y, const double* A, const double* B, int N, int m);
// строим образ х (из отрезка 0..1) по прообразу  _y (принадлежащий гиперкубу)
// нижняя граница области поиска - A
// верхняя граница области поиска - B
// размерность - N
// плотность построения разверток - m
void GetInverseImage( double* _y, double& x, const double* A, const double* B, int N, int m);

#endif //EVOLVENT_H