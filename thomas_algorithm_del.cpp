#include "thomas_algorithm.h"
#include <iostream>
using namespace std;

/* (n by n diagonally dominant matrix) Gauss seidel solver
 *                    System
 *            A       *     x     =    B
 *         L * U      *     x     =    B
 *
 *  (a11 a12  0  ...     0)      (x1)        (b1)
 *  (a21 a22 a23 ...     0)      (x2)        (b2)
 *  (    a32 a33 a34     0)       .            .
 *  (                     )  *    .      =     .
 *  (                     )       .            .
 *  (            ann-1 ann)      (xn)        (bn)
*/

template <class T>
thomas_algorithm<T>::thomas_algorithm()
{


}
template<class T>
int thomas_algorithm<T>::thomas_algorithm_solver()
{
    T* L, *U, *y, *B, *x;
    int n =4;
    T *A=new T[n*n];
    for(int i=0;i<n;i++)
    {
        A[i*(n+1)]=2.04;
        if(i<n-1)
        {
            A[(i+1)*n+i]=-1;
            A[(i*n)+1+i]=-1;
        }

    }
    L=new T[n*n];U=new T[n*n];y=new T[n];B=new T[n], x=new T[n];


    B[0]=40.8;B[1]=0.8;B[2]=0.8;B[3]=200.8;


    //this->template LU_Decomposition<T>(A,L,U, n);

    y[0]=B[0];
    for(int i=0;i<n;i++)
    {
        y[i]=B[i]-L[(n*i)+(i-1)]*y[i-1];
    }
    x[n-1]=y[n-1]/U[n*n-1];
    for(int i=n-2;i>=0;i--)
    {
        x[i] = (y[i]-U[(n+1)*i+1]*x[i+1])/U[(n+1)*i];
    }

    for(int i=0;i<n;i++)
        cout<<x[i]<<endl;
    delete []L; delete[]U;
    delete []y; delete[]B; delete[]x;

    return 1;
}
template <class T>
int thomas_algorithm<T>::LU_Decomposition(T *A, T *L, T *U, int n)
{

    U[0]=A[0];U[1]=A[1];L[0]=1;
    for(int i=1;i<n;i++)
    {
        L[(n*i)+(i-1)]=A[(n*i)+(i-1)]/U[(i-1)*(n+1)];       //4,9,14,  ...
        L[i*(n+1)] = 1;                                     //0,5,10,15...
        U[i*(n+1)]=A[i*(n+1)]- L[(n*i)+(i-1)]*U[n*(i-1)+i]; //0,5,10,15...
        if(i<n-1)
            U[(n*i)+i+1]=A[(n*i+i+1)];                      //1,6,11...
    }





return 1;



}
