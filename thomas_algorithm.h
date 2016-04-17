#ifndef THOMAS_ALGORITHM_H
#define THOMAS_ALGORITHM_H

template<class T>
class thomas_algorithm
{
public:
    thomas_algorithm(int );
    ~thomas_algorithm();
    int thomas_algorithm_solver(T *x, T *B);
    int LU_Decomposition(T *A);

private:
//    int LU_Decomposition(T *A, T *L, T *U, int);
protected:
    static T *L,*U,*y;static int n;
};
#include "thomas_algorithm.cpp"
#endif // THOMAS_ALGORITHM_H
