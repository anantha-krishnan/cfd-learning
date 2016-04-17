#ifndef COUETTE_FLOW_H
#define COUETTE_FLOW_H
//#include <thomas_algorithm.h>
#include <numerical_methods.h>

using namespace std;

class couette_flow
{
public:
    couette_flow(int );
    ~couette_flow();
    void set_BC_cond(double , double );
    void set_initial_cond(double);
    void set_constantBmatrix();
    void set_U(double *);
    void couette_flow_solver();
    void construction(double , double , double );
    double *U_n,*A,*B,length, Re,del_y,del_t, E, initial_first, inital_last;
    int matrix_size, div;




private:






};

#endif // COUETTE_FLOW_H
