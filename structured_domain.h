#ifndef STRUCTURED_DOMAIN_H
#define STRUCTURED_DOMAIN_H

class structured_domain
{
public:
    structured_domain(int d):dimension(d){}
    //virtual void initialise_memory()=0;
protected:
    int dimension;
};

class structured_domain1D:public structured_domain
{
public:
    structured_domain1D():structured_domain(1){}



protected:
    virtual void initialise_memory(int);
    int points;
    double *rhou_vel_1D,*p_1D,*p_star_1D, *p_prime_1D;

};

class structured_domain2D:structured_domain
{
public:
    structured_domain2D():structured_domain(2){}
    virtual void initialise_memory(int, int);


protected:
    virtual ~structured_domain2D();//virtual required so as to call the destructors in a sequence. first derived class then this base class. protected to ensure that base class pointer alone is not called diectly thru derived class object- will result in memory leak thru derived class not getting destroyed.
    int m_rows, n_columns;
    double **rhou_vel_2D, **rhov_vel_2D, **p_2D, **p_star_2D, **p_prime_2D, **u_star_2D, **v_star_2D;
};


#endif // STRUCTURED_DOMAIN_H
