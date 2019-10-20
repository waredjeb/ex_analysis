#include "prng.h"
#include "xorshiro.h"
#include "xorshiroGen.h"
#include "simulation.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <algorithm>
#include "TRandom3.h"
#include <ctime>

simulation::simulation(int N, double sm_fact_n,double lambda_n,double diameter_n, double I0_n)

{
    N_ev = N;
    sm_fact = sm_fact_n;
    lambda = lambda_n;
    diameter = diameter_n;
    I0 = I0_n;

}

simulation::~simulation(){
    //delete clean;
    //delete smeared;
}

double simulation::diffraction_fun(double theta){

    double k = 2*M_PI/lambda;
    double alpha = k*diameter*sin(theta);
    double J_x = gsl_sf_bessel_J1(alpha);

    return (I0*(J_x)/alpha*2)*((J_x)/alpha*2);

}

double* simulation::diffraction_tc(double xmin,double xmax,
                                   double ymin,double ymax){

        clean = new double[N_ev];
        int i = 0;

        while(i < N_ev){
            double x = uniform.rand(xmin,xmax);
            double y = uniform.rand(ymin,ymax);

            //std::cout << "y: " << y << " diff(x)" << diffraction_fun(x) << std::endl;


            if(y < diffraction_fun(x)){

                
                
                
                clean[i] = x;
                i++;
            }
        }
        return clean;
}

double* simulation::diffraction_smeared(double xmin,double xmax, int bins){
    
    double std = sm_fact*(xmax-xmin);
    smeared = new double[N_ev];
    
        for(long int i = 0; i < N_ev; i++){

            double x = clean[i];
            smeared[i] = x + dist.gauss_bm(0,std);


        }
            return smeared;
}
