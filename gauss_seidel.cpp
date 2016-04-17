#include <iostream>
#include "gauss_seidel.h"
#include <math.h>
#include <time.h>
#define tol 1e-5
using namespace std;


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
