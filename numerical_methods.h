#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H
#include <iostream>
#include "constants.h"

class numerical_methods
{
public:
    numerical_methods();
};
class gauss_seidel
{
public:
    int gauss_seidel_method();
private:

};

class thomas_algorithm
{
public:
    thomas_algorithm(int );
    ~thomas_algorithm();
    void A_tridiagonalmatrix_cons(double dia_ele, double upper_dia_ele, double lower_dia_ele);
    void thomas_algorithm_solver(double *x, double *B);
    int LU_Decomposition();

private:
//    int LU_Decomposition(T *A, T *L, T *U, int);
protected:
    static double *A,*L,*U,*y;static int n;
};

#endif // NUMERICAL_METHODS_H
