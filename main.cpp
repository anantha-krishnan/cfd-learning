#include <iostream>

#include "numerical_methods.h"
#include "couette_flow.h"
#include "simple.h"

using namespace std;

//#define extern int tol 1e-5

int main()
{
    /*gauss_seidel * h =new gauss_seidel();
    h->gauss_seidel_method();
    delete h;

    thomas_algorithm<double> *h = new thomas_algorithm<double>(4);
    h->thomas_algorithm_solver();
    delete h;*/

    /*couette_flow *h=new couette_flow(20);
    h->set_BC_cond(0,1);
    h->set_initial_cond(0);
    h->construction(0,5000,1);
    h->set_constantBmatrix();
    h->couette_flow_solver();
    delete h;*/

   SIMPLE2D *h = new SIMPLE2D();
    h->startsimple();
    delete h;

    return 0;
}

