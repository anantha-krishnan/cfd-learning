#include "simple.h"

void SIMPLE2D::initialization(double p_star)
{
    //guess value at all pressure points p*
    for(int j=0;j<m_rows;j++)
    {
        for(int i=0;i<n_columns;i++)
        {
            p_star_2D[i][j]=p_star;p_prime_2D[i][j]=p_star;
            rhou_vel_2D[i][j]=0;rhov_vel_2D[i][j]=0;

        }
        rhou_vel_2D[n_columns][j]=0;
        rhov_vel_2D[n_columns][j]=0;
        rhov_vel_2D[n_columns+1][j]=0;
    }

    for(int i=0;i<n_columns+1;i++)
        rhov_vel_2D[n_columns-1+i][m_rows]=0;


}

void SIMPLE2D::BCs(double u_upper, double v_up, double u_bottom, double v_bottom, double pprime_in, double v_in, double pprime_out)
{
    /*p runs from 1 to n, 1 to m
      u runs from 1 to n+1, 1 to m
      v runs from 1 to n+2, 1 to m+1
      */

    for(int i=0;i<m_rows;i++)
    {
        p_prime_2D[n_columns-1][i]=pprime_out;
        p_prime_2D[0][i]=pprime_in;
        rhov_vel_2D[0][i]=v_in;
    }
    rhov_vel_2D[0][m_rows]=v_in;
    for(int i=0;i<n_columns+1;i++)
    {
        rhov_vel_2D[i][m_rows]=v_up;
        rhov_vel_2D[i][0]=v_bottom;
    }
    for(int i=0;i<n_columns;i++)
    {
        rhou_vel_2D[i][0]=u_bottom;
        rhou_vel_2D[i][m_rows]=u_upper;
    }
}
void SIMPLE2D::prediction()
{
    double A ;
    for(int i=1;i<n_columns-1;i++)
        for(int j=1;j<m_rows-1;j++)
        {
            rhou_vel_2D[i+1][j]=rhou_vel_2D[i+1][j]; + calculate_A_star(i,j);
        }
}
void SIMPLE2D::case_setup(double l, double h, double t, double density, double viscosity)
{
    delx=l/(structured_domain2D::n_columns-1);
    dely=h/(structured_domain2D::m_rows-1);
    delt=t; rho=density; mu=viscosity;
}
void SIMPLE2D::calculate_u_dash_star()
{

}
double SIMPLE2D::calculate_A_star(int i, int j)
{
    double A;
    A = (rho * u_star_2D[i+2][j]*u_star_2D[i+2][j] - rho * u_star_2D[i][j]*u_star_2D[i][j])/(2*delx) + (rho * u_star_2D[i+1][j+1] * ((v_star_2D[i+1][j+1]+v_star_2D[i+2][j+1])/2)-(rho * u_star_2D[i+1][j-1] * (v_star_2D[i+1][j]+v_star_2D[i+2][j])/2))/(2*dely) + mu * ((u_star_2D[i+2][j]- (2*u_star_2D[i+1][j])+u_star_2D[i][j])/(delx*delx) + (u_star_2D[i+1][j+1]- (2*u_star_2D[i+1][j])+u_star_2D[i+1][j-1])/dely*dely );
    return (A);
}
