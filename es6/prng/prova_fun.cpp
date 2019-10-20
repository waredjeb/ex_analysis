#include <iostream>
#include "simulation.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <vector>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

int main(int argc, char** argv){

    //simulation exp;
    int N = 3;
    //int N = 50000;
    double diameter = 0.75e-3;
    double lambda = 400e-9;
    double c = 1;


    double i = -3.11;
    int k = 0;
    double* x_axis = new double[62];
    double* y_axis = new double[62];

    simulation exp(N,c,lambda,diameter);

    while(i<N){
        
        std::cout << k << std::endl;
        y_axis[k]= exp.diffraction_fun(i);;
        x_axis[k]=i;
        i +=0.1;
        k+=1;
        std::cout << "x: " << i << "  y: " << exp.diffraction_fun(i)*pow(10,6) << std::endl;               
    }

    TCanvas* c1 = new TCanvas("c","c",1000,1000);
    TGraph * g = new TGraph(62,x_axis,y_axis);
    g->Draw("APC");
    c1->Draw("");
    c1->Print("./try.png","png");




    return 0;
}