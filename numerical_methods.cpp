#include "numerical_methods.h"
#include <iostream>
#include <math.h>
#include <time.h>
#define tol 1e-5
using namespace std;

numerical_methods::numerical_methods()
{
}



int gauss_seidel::gauss_seidel_method()
{
    /* (n by n diagonally dominant matrix) Gauss seidel solver
     *                    System
     *            A       *     x     =    B
     *        Amatrix     *  Xmatrix  = Bmatrix
     *
     *  (a11 a12 ... a1n)      (x1)        (b1)
     *  (a21 a22 ... a2n)      (x2)        (b2)
     *  ( .   .   .   . )       .            .
     *  ( .   .   .   . )  *    .      =     .
     *  ( .   .   .   . )       .            .
     *  (an1 an2 ... ann)      (xn)        (bn)
   */
    int i=0,n=3;double sum1=0.0, iter=0;


    double *Amatrix, *Xmatrix, *Bmatrix, temp=0.0, time_taken=0.0; bool residual = true;
    clock_t tStart;

    Xmatrix  = new double[n];
    Amatrix  = new double[n*n];
    Bmatrix  = new double[n];

    for(i=0;i<n;i++)
        Xmatrix[i]=0;

    Amatrix[0] = 4;Amatrix[1] = 1;Amatrix[2] = -1;
    Amatrix[3] = 2;Amatrix[4] = 7;Amatrix[5] = 1;
    Amatrix[6] = 1;Amatrix[7] = -3;Amatrix[8] = 12;

    Bmatrix[0]=3;Bmatrix[1]=19;Bmatrix[2]=31;

    tStart =  clock();

    do
    {
        for(int m=0;m<n;m++)
        {
            for(i=0;i<n;i++)
            {
                if(i!=m)
                {
                    sum1+=Amatrix[(m*n)+i]*Xmatrix[i];
                }

            }

            temp = Xmatrix[m];
            Xmatrix[m]=1/Amatrix[(n+1)*m]*(Bmatrix[m] - sum1 );

            sum1=0.0;

            if(fabs(temp - Xmatrix[m]) < tol)
                residual = false;
            iter++;

        }

    }while(residual);

    time_taken = (double)(clock() - tStart)/double(CLOCKS_PER_SEC);

    for(i=0;i<n;i++)
        cout<<"Xmatrix["<<i<<"]"<<" = "<<Xmatrix[i]<<endl;
    cout<<"time taken (sec) = "<<time_taken<<", iterations = "<<iter<<endl;

    delete []Xmatrix; delete []Bmatrix; delete []Amatrix;

    return 1;

}

double *thomas_algorithm::L;
double *thomas_algorithm::U;
double *thomas_algorithm::y;
int thomas_algorithm::n;

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


thomas_algorithm::thomas_algorithm(int size)
{
    L=new double[size*size]();U=new double[size*size]();y=new double[size]();
     n = size;


}

thomas_algorithm::~thomas_algorithm()
{
    delete []L;delete []U; delete []y;
    L=NULL;U=NULL;y=NULL;
}


void thomas_algorithm::thomas_algorithm_solver(double *U_n, double *B)
{


    y[0]=B[0];

    for(int i=1;i<n;i++)
    {

        y[i]=B[i]-L[(n*i)+(i-1)]*y[i-1];

    }

    U_n[n-1]=y[n-1]/U[n*n-1];

    for(int i=n-2;i>=0;i--)
    {

        U_n[i] = (y[i]-U[(n+1)*i+1]*U_n[i+1])/U[(n+1)*i];

    }
}

int thomas_algorithm::LU_Decomposition(double *A)
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


