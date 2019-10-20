#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include <iostream>
#include "TF2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include <cmath>

using namespace RooFit ;

 
struct estimates{ //struct to contain result of fit
    double a_est;
    double b_est;
    double g_est;
};

struct mc_generated{
    std::vector<double> theta;
    std::vector<double> phi;
};

void ml_unbinned_fit(mc_generated a){
    
    RooRealVar  theta("theta","theta",0,2*M_PI) ;
    RooRealVar  phi("phi","phi", 0, 2*M_PI) ;
    
    RooDataSet data("data","data",RooArgSet(theta, phi)) ;
    
    for (int i=0 ; i < a.theta.size() ; i++) {
        
        theta = a.theta[i];
        phi = a.phi[i];
        data.add(RooArgSet(theta,phi)) ;
    }
    
    std::cout<< "Filling complete" << "\n";
    RooRealVar alpha("alpha", "alpha", 0.65, 0,1);
    RooRealVar beta("beta", "beta", 0.06, -.5, .5 );
    RooRealVar gamma("gamma", "gamma", -0.18, -1, 1);
    RooRealVar scale("scale", "scale", 5, -10, 10);

    RooWorkspace* ws = new RooWorkspace("ws") ;
    ws->factory("EXPR::mypdf('scale*(3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(theta)*cos(theta) - b*sin(theta)*sin(theta)*cos(2.*phi)- sqrt(2.)*gamma*sin(2.*theta)*cos(phi))',theta, phi, alpha, beta, gamma, scale)");
    //RooGenericPdf genpdf("gen_pdf", "gen_pdf", " scale*(3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(theta)*cos(theta) - b*sin(theta)*sin(theta)*cos(2.*phi)- sqrt(2.)*gamma*sin(2.*theta)*cos(phi))", RooArgSet(theta, phi, alpha, beta, gamma, scale));
    
    //Generating toy MC:
    //RooDataSet* data_roofit = genpdf.generate(RooArgSet(theta,phi),50000);
    //TH2F* h_pdf = new TH2F("h2", "h2", 100, 0, 2*M_PI, 100, 0, 2*M_PI);

    //for(int i=0; i < 50000; i++){
        //data_roofit->get(i)->Print("v");
    //}   


    return;
}


double to_fit(double *x, double *par){ //to be passed to the fitting routine
  return par[3]*(3./(4.*M_PI))*(0.5*(1.-par[0]) + (0.5)*(3.*par[0]-1)*cos(x[0])*cos(x[0]) - par[1]*sin(x[0])*sin(x[0])*cos(2.*x[1])- sqrt(2.)*par[2]*sin(2.*x[0])*cos(x[1]));
}

double evaluate(double t, double p){ // this evaluates the distribution at theta and phi

    double a = 0.65;
    double b = 0.06;
    double g = -0.18;

    return (3./(4.*M_PI))*(0.5*(1.-a) + (0.5)*(3.*a-1)*cos(t)*cos(t) - b*sin(t)*sin(t)*cos(2.*p)- sqrt(2.)*g*sin(2.*t)*cos(p));
}

estimates monte_carlo(int measures, int seed ){
    
    std::vector<double> theta_vec;
    std::vector<double> phi_vec;

    TRandom3* gen = new TRandom3(seed);
    TH2F* h = new TH2F("h", "h", 300, 0, 2*M_PI, 300, 0, 2*M_PI);

    for(int i = 0; i < measures; i++){ //rejection method to build distribution

        bool it = true; 

        while(it){

            double theta = gen->Uniform(0,2*M_PI);
            double phi = gen->Uniform(0, 2*M_PI);
            double y = gen->Uniform(0,1);

            if(y < evaluate(theta, phi)){
                theta_vec.push_back(theta);
                phi_vec.push_back(phi);
                h->Fill(theta, phi);
                it = false;
            }
        }
    }

    mc_generated simu;
    simu.theta = theta_vec;
    simu.phi = phi_vec;

    ml_unbinned_fit(simu);

    const int npar = 4;
    double dist_params[npar] = {0.65, 0.06, -0.18, 1.};
    dist_params[3] = h->Integral("width") / 4.; // "#scaling" estimation
    std::cout << "Par 3" << dist_params[3] << "\n";

    TF2* fit = new TF2("fit", to_fit, 0, 2*M_PI, 0, 2*M_PI, npar);
    fit->SetParameters(dist_params);
    fit->SetParNames("#alpha", "#beta", "#gamma", "#scaling");

    
    h->Fit(fit, "LL");
    estimates est;
    std::cout<< fit->GetParameter(0) << "\n";
    est.a_est = fit->GetParameter(0);
    est.b_est = fit->GetParameter(1);
    est.g_est = fit->GetParameter(2);

    //TCanvas* c = new TCanvas("c", "c", 1000,1000, 1000,800);
    //h->Draw("colz");
    //fit->Draw("same cont1");
    //c->Draw();
    //c->SaveAs("./miao.png", "png");

    return est;
}

int main(){

    srand(time(NULL));

    int experiments = 1;
    int measures = 50000;

    //this is to check CLT on parameters.

    TH1F* h_a = new TH1F("#alpha parameter distribution", "#alpha parameter distribution", 100, 0.63, 0.67);
    h_a->SetLineWidth(2);
    h_a->GetXaxis()->SetTitle("#alpha");
    h_a->GetYaxis()->SetTitle("counts");

    TH1F* h_b = new TH1F("#beta parameter distribution", "#beta parameter distribution", 100, 0.05, 0.07);
    h_b->SetLineWidth(2);
    h_b->GetXaxis()->SetTitle("#beta");
    h_b->GetYaxis()->SetTitle("counts");

    TH1F* h_g = new TH1F("#gamma parameter distribution", "#gamma parameter distribution", 100, -0.17, -0.19);
    h_g->SetLineWidth(2);
    h_g->GetXaxis()->SetTitle("#gamma");
    h_g->GetYaxis()->SetTitle("counts");

    TF1* fit_a = new TF1("fit_a", "gaus", 0.63, 0.67);
    TF1* fit_b = new TF1("fit_b", "gaus", 0.05, 0.07);
    TF1* fit_g = new TF1("fit_g", "gaus", -0.17, -0.19);

    for(int i = 0; i < experiments; i++){
        if(i%10 == 0 ){
            std::cout << i << "\n";
        }
        int seed = rand();
        estimates est =  monte_carlo(measures, seed );
        h_a->Fill(est.a_est);
        h_b->Fill(est.b_est);
        h_g->Fill(est.g_est);

    }

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    TCanvas*c = new TCanvas("c", "c", 1000,1000, 1000, 800);
    h_a->Fit(fit_a);
    h_a->Draw("hist");
    fit_a->Draw("same");
    c->Draw();
    c->SaveAs("./a.png", "png");

    TCanvas*c1 = new TCanvas("c1", "c1", 1000,1000, 1000, 800);
    h_b->Fit(fit_b);
    h_b->Draw("hist");
    fit_b->Draw("same");
    c1->Draw();
    c1->SaveAs("./b.png", "png");

    TCanvas*c2 = new TCanvas("c2", "c2", 1000,1000, 1000, 800);
    h_g->Fit(fit_g);
    h_g->Draw("hist");
    fit_g->Draw("same");
    c2->Draw();
    c2->SaveAs("./g.png", "png");
    
    

    return 0;;
    
}
