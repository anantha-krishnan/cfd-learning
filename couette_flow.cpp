#include "couette_flow.h"


couette_flow::couette_flow(int div)
{
    length=1.0;
    del_y=length/div;
    matrix_size=div+1-2;
    B=new double[matrix_size];U_n=new double[matrix_size];
}
couette_flow::~couette_flow()
{
    delete []B; delete []U_n;
}

void couette_flow::construction(double del_t,double Re)
{
    //double a, b;E=e;
    if(del_t==0)
        del_t = E*Re*del_y*del_y;
     //a = -E/2;b=1+E;


}


void couette_flow::set_BC_cond(double a, double b)
{
    initial_first=a; inital_last=b;
}

void couette_flow::set_initial_cond(double a)
{
    for(int i=0;i<matrix_size;i++)
        U_n[i]=a;
}

void couette_flow::set_constantBmatrix()
{
    //second point after BC1 => B[2]=(1-E)*initial_first
    B[0]=(1-E)*U_n[0] + (E/2)*(U_n[1] + initial_first) - (-E/2 * initial_first);
    for(int i=1;i<matrix_size-1;i++)
        B[i]=(1-E)*U_n[i] + (E/2)*(U_n[i+1]+U_n[i-1]);
    B[matrix_size-1]=(1-E)*U_n[matrix_size-1] + (E/2)*(inital_last + U_n[matrix_size-2]) - (-E/2*inital_last);

}

void couette_flow::set_U(double *updated_U)
{
    U_n=updated_U;
}

void couette_flow::couette_flow_solver()
{
    thomas_algorithm *h = new thomas_algorithm(matrix_size);
    h->A_tridiagonalmatrix_cons(1+E,-E/2,-E/2);

    h->LU_Decomposition();

    for(int i=0;i<1361;i++)
    {
        h->thomas_algorithm_solver(U_n, B);
        set_constantBmatrix();


    }
    cout<<initial_first<<endl;
    for(int i=0;i<matrix_size;i++)
        cout<<U_n[i]<<endl;
    cout<<inital_last<<endl;
    delete h;



}

