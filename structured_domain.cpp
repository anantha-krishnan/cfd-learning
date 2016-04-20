#include "structured_domain.h"

//structured_domain1D::structured_domain1D()
//{
//    dimension = 0;
//    points=0;
//    rhou_vel_1D=0,p_1D=0,p_star_1D=0, p_prime_1D=0;


//}
void structured_domain1D::initialise_memory(int discretepoints)
{
    dimension = 1; points=discretepoints;
    rhou_vel_1D=new double[points];p_1D=new double[points];p_star_1D=new double[points];p_prime_1D=new double[points];

}

//structured_domain2D::structured_domain2D()
//{
//    dimension =2; m_rows=0;n_columns=0;
//    rhou_vel_2D=0, rhov_vel_2D=0, p_2D=0, p_star_2D=0, p_prime_2D=0;
//}

void structured_domain2D::initialise_memory(int xpoints, int ypoints)
{
    dimension=2;n_columns=xpoints,m_rows=ypoints;
    rhou_vel_2D=new double*[n_columns+1];rhov_vel_2D=new double*[n_columns+2];
    u_star_2D=new double*[n_columns+1];v_star_2D=new double*[n_columns+2];
    p_2D=new double*[n_columns];p_star_2D=new double*[n_columns];p_prime_2D=new double*[n_columns];
    for(int i=0;i<n_columns;i++)
    {        
        p_2D[i]= new double[m_rows];
        p_star_2D[i]=new double[m_rows];
        p_prime_2D[i]=new double[m_rows];
    }
    for(int i=0;i<n_columns+1;i++)
    {
        u_star_2D[i]=new double[m_rows];
        rhou_vel_2D[i]=new double[m_rows];
    }
    for(int i=0;i<n_columns+2;i++)
    {
        v_star_2D[i]=new double[m_rows+1];
        rhov_vel_2D[i]=new double[m_rows+1];
    }
    for(int j=0;j<m_rows;j++)
    {
        for(int i=0;i<n_columns;i++)
        {
            p_star_2D[i][j]=0;p_prime_2D[i][j]=0;
        }
    }
    for(int j=0;j<m_rows;j++)
    {
        for(int i=0;i<n_columns+1;i++)
        {
            u_star_2D[i][j]=0;rhou_vel_2D[i][j]=0;
        }
    }
    for(int j=0;j<m_rows+1;j++)
    {
        for(int i=0;i<n_columns+2;i++)
        {
            v_star_2D[i][j]=0;
            rhov_vel_2D[i][j]=0;
        }
    }
}
structured_domain2D::~structured_domain2D()
{
    for(int i=0;i<n_columns;i++)
    {
        delete []p_2D[i];
        delete []p_star_2D[i];
        delete []p_prime_2D[i];
    }

    delete []p_2D;
    delete []p_star_2D;
    delete []p_prime_2D;

    for(int i=0;i<n_columns+1;i++)
    {
        delete []u_star_2D[i];
        delete []rhou_vel_2D[i];
    }
    delete []u_star_2D;
    delete []rhou_vel_2D;

    for(int i=0;i<n_columns+2;i++)
    {
        delete []v_star_2D[i];
        delete []rhov_vel_2D[i];
    }

    delete []v_star_2D;
    delete []rhov_vel_2D;





}
