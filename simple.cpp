#include "simple.h"
SIMPLE2D::~SIMPLE2D()
{

}

void SIMPLE2D::getvalues()
{
    double u_upper,  v_up,  u_bottom,  v_bottom,  pprime_in,  v_in,  pprime_out;
//    cout<<"enter upper u vel,upper v vel, bottom u,bottom v, v inlet , xpoints, y points"<<endl;
//    cin>>u_upper>>v_up>>u_bottom>>v_bottom>>v_in>>x>>y;
    pprime_in=0;pprime_out=0;u_upper=1;u_bottom=0;v_up=0;v_bottom=0;v_in=0;
//    cout<<"x = "<<x<<"y = "<<y<<endl;
    initialise_memory(21,11);
    BCs(u_upper,v_up,u_bottom,v_bottom,pprime_in,v_in,pprime_out);
    double vis;
//    cout<<"enter l, h, delt, density, viscosity"<<endl;
//    cin>>l>>h>>t>>den>>vis;
    vis=3.7e-7;
    case_setup(0.5,0.01,0.001,0.002733,vis);

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

    }
    for(int i=0;i<m_rows+1;i++)
    {
    v_star_2D[0][i]=v_in;
    }
    for(int i=0;i<n_columns+2;i++)
    {
        v_star_2D[i][m_rows]=v_up;
        v_star_2D[i][0]=v_bottom;
    }
    for(int i=0;i<n_columns+1;i++)
    {
        u_star_2D[i][0]=u_bottom;
        u_star_2D[i][m_rows-1]=u_upper;
    }
}

void SIMPLE2D::case_setup(double l, double h, double t, double density, double viscosity)
{
    delx=l/(structured_domain2D::n_columns-1);
    dely=h/(structured_domain2D::m_rows-1);
    delt=t; rho=density; mu=viscosity;
    v_star_2D[14][4]=0.5;
    rhov_vel_2D[14][4]=v_star_2D[14][4]*rho;
}

void SIMPLE2D::calculate_A_star(double *A,int i, int j)
{
    //double A;
    /*double t1,t2,t3,t4;
    t1=rho * u_star_2D[i+1][j]*u_star_2D[i+1][j] ;
    t2=rho * u_star_2D[i-1][j]*u_star_2D[i-1][j];
    t3=rho * u_star_2D[i][j+1];
    t4=v_star_2D[i][j]+v_star_2D[i+1][j]/2;
    */
    *A = (rho * u_star_2D[i+1][j]*u_star_2D[i+1][j] - rho * u_star_2D[i-1][j]*u_star_2D[i-1][j])/(2*delx) + (rho * u_star_2D[i][j+1] * ((v_star_2D[i+1][j+1]+v_star_2D[i][j+1])/2)-(rho * u_star_2D[i][j-1] * (v_star_2D[i+1][j-1]+v_star_2D[i][j-1])/2))/(2*dely) + mu * ((u_star_2D[i+1][j]- (2*u_star_2D[i][j])+u_star_2D[i-1][j])/(delx*delx) + (u_star_2D[i][j+1]- (2*u_star_2D[i][j])+u_star_2D[i][j-1])/(dely*dely) );
    //return (A);
}
void SIMPLE2D::calculate_B_star(double *B, int i, int j)
{
    *B= ((rho*v_star_2D[i+1][j]*(u_star_2D[i][j]+u_star_2D[i][j-1])/2)-(rho*v_star_2D[i-1][j]*(u_star_2D[i-1][j-1]+u_star_2D[i-1][j])/2))/(2*delx) +((rho*v_star_2D[i][j+1]*v_star_2D[i][j+1])-(rho*v_star_2D[i][j-1]*v_star_2D[i][j-1]))/(2*dely) + mu * ((v_star_2D[i+1][j]-2*v_star_2D[i][j]+v_star_2D[i-1][j])/(delx*delx) + (v_star_2D[i][j+1]-2*v_star_2D[i][j]+v_star_2D[i][j-1])/(dely*dely)) ;
}
void SIMPLE2D::prediction()
{
    double A,B ;
    for(int i=1;i<n_columns;i++)//u is calculated from one to n+1 columns
        for(int j=1;j<m_rows-1;j++)//one to m rows
        {
            calculate_A_star(&A,i,j);
            rhou_vel_2D[i][j]=rhou_vel_2D[i][j] + A * delt - ((p_star_2D[i][j]-p_star_2D[i-1][j])* delt/delx);
            u_star_2D[i][j]=rhou_vel_2D[i][j]/rho;
        }
    for(int i=1;i<n_columns+1;i++)//v is calculated from one to n+2 columns
        for(int j=1;j<m_rows;j++)//one to m+1 rows
        {
            calculate_B_star(&B, i, j);
            rhov_vel_2D[i][j]=rhov_vel_2D[i][j] + B *delt - ((p_star_2D[i-1][j]-p_star_2D[i-1][j-1])*delt/dely);
            v_star_2D[i][j]=rhov_vel_2D[i][j]/rho;
        }

    for(int i=0;i<m_rows;i++)
    {
        u_star_2D[0][i]=u_star_2D[1][i];
        u_star_2D[n_columns][i]=u_star_2D[n_columns-1][i];
        v_star_2D[n_columns+1][i]=v_star_2D[n_columns][i];
    }
        v_star_2D[n_columns+1][m_rows]=v_star_2D[n_columns][m_rows];

}
void SIMPLE2D::correction()
{
    //    clock_t tStart =  clock();
    double b = -delt/(delx*delx); double c = -delt/(dely*dely);
    double a = 2*(b+c);
    double d /*,temp*/;/*bool residue=true*/;int iter=0,m=0;

    thomas_algorithm *h=new thomas_algorithm(m_rows-2);
    h->A_tridiagonalmatrix_cons(a,c,c);h->LU_Decomposition();


    double *B=new double[m_rows-2]();double *P=new double[m_rows-2]();

    //do
    //for(m=0;m<300;m++)
    {
        for(int i=1;i<n_columns-1;i++)
        {

            for(int j=1;j<m_rows-1;j++)
            {
                d=(rhou_vel_2D[i+1][j] - rhou_vel_2D[i][j])/delx + (rhov_vel_2D[i+1][j+1] - rhov_vel_2D[i+1][j])/dely;

                B[j-1]=(d + b*p_prime_2D[i+1][j] + b*p_prime_2D[i-1][j]);
                if(j==1)
                    B[j-1]+= c*p_prime_2D[i][j-1];
                if(j==m_rows-2)
                    B[j-1]+= c*p_prime_2D[i][j+1];
                B[j-1]=-1/a*B[j-1];

                P[j-1]=p_prime_2D[i][j];
            }

            h->thomas_algorithm_solver(P,B);

            for(int j=1;j<m_rows-1;j++)
                p_prime_2D[i][j]=P[j-1];
            /*temp = p_prime_2D[i][j];
            //cout<<"temp = "<<temp<<endl;


            p_prime_2D[i][j]=-(b*p_prime_2D[i+1][j] + b*p_prime_2D[i-1][j] + c*p_prime_2D[i][j+1] + c*p_prime_2D[i][j-1] + d)/a ;
            if(fabs(temp - p_prime_2D[i][j]) < TOL)
                residue = false;
            else
                residue = true;*/
            iter++;
        }
    }//while(residue);
    //cout<<"over "<<endl;
    //cout<<"m "<<m<<endl;
    delete h;
    delete []B;
    delete []P;

}


    //    time_taken = (double)(clock() - tStart)/double(CLOCKS_PER_SEC);

    //    for(i=0;i<n;i++)
    //        cout<<"Xmatrix["<<i<<"]"<<" = "<<Xmatrix[i]<<endl;
    //    cout<<"time taken (sec) = "<<time_taken<<", iterations = "<<iter<<endl;



void SIMPLE2D::correctedpressure()
{

    for(int i=1;i<n_columns-1;i++)
        for(int j=1;j<m_rows-1;j++)
            p_star_2D[i][j]=p_star_2D[i][j] + underrelaxationfactor * p_prime_2D[i][j];


    //for(int i=0;i<m_rows+1;i++)
    //cout<<" vel u "<<v_star_2D[n_columns-2][i]<<endl;
    //cout<<v_star_2D[14][4]<<endl;



}
void SIMPLE2D::startsimple()
{
    getvalues();

    for(int i=0;i<300;i++)
    {
        prediction();
        correction();
        correctedpressure();
    }

    for(int i=0;i<m_rows;i++)
        cout<<u_star_2D[14][i]<<endl;
}
