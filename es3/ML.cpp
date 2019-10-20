#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include <iostream>
#include "ML.h"

using namespace RooFit ;

/* double to_fit(double *x, double *par){ //to be passed to the fitting routine
    return par[3]*(3./(4.*M_PI))*(0.5*(1.-par[0]) + (0.5)*(3.*par[0]-1)*cos(x[0])*cos(x[0]) - par[1]*sin(x[0])*sin(x[0])*cos(2.*x[1])- sqrt(2.)*par[2]*sin(2.*x[0])*cos(x[1]));
}

int main()
{
    
    RooRealVar x("x","x",-20,20) ;
    
    // Model (intentional strong correlations)
    RooRealVar mean("mean","mean of g1 and g2",0) ;
    RooRealVar sigma_g1("sigma_g1","width of g1",3) ;
    RooGaussian g1("g1","g1",x,mean,sigma_g1) ;
    
    RooRealVar sigma_g2("sigma_g2","width of g2",4,3.0,6.0) ;
    RooGaussian g2("g2","g2",x,mean,sigma_g2) ;
    
    RooRealVar frac("frac","frac",0.5,0.0,1.0) ;
    RooAddPdf model("model","model",RooArgList(g1,g2),frac) ;
    
    // Generate 1000 events
    RooDataSet* data = model.generate(x,1000) ;
    
    for(int i = 0; i <1000; i++){
        data->get(i)->Print("v");
        std::cout << std::endl;
    }
    
    
    return 0;
}
*/

void roofit_methods::ml_unbinned_fit(){

    RooRealVar  theta("theta","theta",0,2*M_PI) ;
    RooRealVar  phi("phi","phi", 0, 2*M_PI) ;

    RooDataSet data("data","data",RooArgSet(theta, phi)) ;

    for (int i=0 ; i < vectors.theta.size() ; i++) {
        
        theta = vectors.theta[i];
        phi = vectors.phi[i];
        data.add(RooArgSet(theta,phi)) ;
    }

    std::cout<< "Filling complete" << "\n";
    
    return;
}
