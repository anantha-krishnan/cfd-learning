#ifndef SIMPLE_H
#define SIMPLE_H

#include "numerical_methods.h"
#include "structured_domain.h"
#include "constants.h"
#include <math.h>
#include <cmath>
#include <iostream>
using namespace std;

class SIMPLE1D:public structured_domain1D
{
public:
    SIMPLE1D():structured_domain1D(){}



protected:


};

class SIMPLE2D:public structured_domain2D
{
public:
    SIMPLE2D():structured_domain2D(){}

    void startsimple();
    void getvalues();
    void BCs(double,double,double,double,double,double,double);
    void prediction();
    void correction();
    void correctedpressure();
    void case_setup(double l, double h, double t, double , double);
    void calculate_A_star(double *,int i, int j);
    void calculate_B_star(double *, int ,int);
    void calculate_u_dash_star();
    void calculate_v_dash_star();
    virtual ~SIMPLE2D();
protected:
    double delx, dely, delt,rho, mu;
};

#endif // SIMPLE_H
