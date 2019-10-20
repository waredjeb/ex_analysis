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
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "RooClassFactory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooMCStudy.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include <cmath>
#include "string.h"
#include <iostream>

using namespace RooFit ;

double myPDF(double *x, double *par){ //to define TF2 obj
  return par[3]*(3./(4.*M_PI))*(0.5*(1.-par[0]) + (0.5)*(3.*par[0]-1)*cos(x[0])*cos(x[0]) - par[1]*sin(x[0])*sin(x[0])*cos(2.*x[1])- sqrt(2.)*par[2]*sin(2.*x[0])*cos(x[1]));
}


int main(){
    
    //defining Real variables according to RooFit
    //Angles physical defined in [0-2*pi]
    RooRealVar  t("t","t",0,2*M_PI) ;
    RooRealVar  p("p","p", 0, 2*M_PI) ;
    
    //Defining Parameters according to RooFit.
    //Initial value is the 3rd argument, 4th and 5th arg
    //are the range of existance of the parameter.
    //RooRealVar alpha("alpha", "alpha", 0.65, -1,1);
    RooRealVar alpha("alpha", "alpha", 0.65, 0.62, 0.66);

    //RooRealVar beta("beta", "beta", 0.06, -1, 1);
    RooRealVar beta("beta", "beta", 0.06, 0.04, 0.075);
    RooRealVar gamma("gamma", "gamma", -0.18, -.2, -0.16);
    RooRealVar scale("scale", "scale", 18, 1, 25.);
    
    //ClassFactory generates a pdf instance and normalize it.
    RooAbsPdf* genpdf = RooClassFactory::makePdfInstance("GenPdf","scale*(3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(t)*cos(t) - beta*sin(t)*sin(t)*cos(2.*p)- sqrt(2.)*gamma*sin(2.*t)*cos(p))",RooArgSet(t,p,alpha, beta, gamma, scale)) ;
    
    // Generate a toy MC dataset from the interpreted p.d.f
    //data2 will contain unbinned data.
    RooDataSet* data = genpdf->generate(RooArgSet(t,p),50000) ;
    
    //visualizing the events generated through a ROOT TH2F.
    TH2F* dh2 = new TH2F("h2", "h2", 100, 0, 2*M_PI, 100, 0, 2*M_PI);
    data->fillHistogram(dh2, RooArgList(t,p));
    TCanvas*c = new TCanvas("c", "c", 1000,1000,1000,800);
    dh2->SetTitle("Binned #theta, #phi distribution");
    dh2->GetXaxis()->SetTitle("#theta [rad]");
    dh2->GetYaxis()->SetTitle("#phi [rad]");
    dh2->Draw("colz");
    c->Draw();
    c->SaveAs("./generated_distribution.png", "png");

    //Unbinned maximum likelihood fit:
    genpdf->fitTo(*data) ;

    //Retireve statistics:
    auto alpha_est = alpha.getVal();
    auto beta_est = beta.getVal();
    auto gamma_est = gamma.getVal();
    auto scale_est = scale.getVal();

    //defining a TF2 object with estimated max LL parameters
    const int npar = 4;
    double dist_params[npar] = {alpha_est, beta_est, gamma_est, scale_est};
    TF2* myPdf = new TF2("myPdf", myPDF, 0, 2*M_PI, 0, 2*M_PI, npar);
    myPdf->SetParameters(dist_params);
    myPdf->SetParNames("#alpha", "#beta", "#gamma", "#scaling");

    gStyle->SetOptStat(1111);
    TCanvas* from_ML = new TCanvas("from_ML", "from_ML", 1000,1000,1000,800);
    dh2->SetTitle(Form("Unbinned ML fit: #alpha = %.2f, #beta = %.2f, #gamma = %.2f",alpha_est, beta_est, gamma_est));
    dh2->GetXaxis()->SetTitle("#theta [rad]");
    dh2->GetYaxis()->SetTitle("#phi [rad]");
    dh2->Draw("colz");
    myPdf->Draw("same cont1");
    from_ML->Draw();
    from_ML->SaveAs("./result_ML_unbinned.pdf", "pdf");

    TCanvas* from_ML_2 = new TCanvas("from_ML_2", "from_ML_2", 1000,1000,1000,800);
    dh2->SetTitle(Form("Unbinned ML fit: #alpha = %.2f, #beta = %.2f, #gamma = %.2f",alpha_est, beta_est, gamma_est));
    dh2->GetXaxis()->SetTitle("#theta [rad]");
    dh2->GetYaxis()->SetTitle("#phi [rad]");
    dh2->Draw("lego");
    myPdf->Draw("same cont1");
    from_ML_2->Draw();
    from_ML_2->SaveAs("./result_ML_unbinned_2.pdf", "pdf");


    //obtaining Likelihhod from model.
    RooAbsReal* nll = genpdf->createNLL(*data, NumCPU(4));

    //instantiating minimizer for likelihood
    RooMinimizer m(*nll) ;
    m.migrad();

    //plotting likelihood for parameter alpha
    RooPlot* frame1 = alpha.frame(Title("Likelihood #alpha")) ;
    nll->plotOn(frame1,ShiftToZero()) ;

    //plotting likelihood for parameter beta
    RooPlot* frame2 = beta.frame(Title("Likelihood #beta")) ;
    nll->plotOn(frame2,ShiftToZero()) ;

    //plotting likelihood for parameter gamma
    RooPlot* frame3= gamma.frame(Title("Likelihood #gamma")) ;
    nll->plotOn(frame3,ShiftToZero()) ;
    
    TCanvas*c1 = new TCanvas("c1", "c1", 1000, 1000, 1000, 800);
    frame1->Draw();
    c1->Draw();
    c1->SaveAs("./likelihood_alpha.png", "png");

    TCanvas*c2 = new TCanvas("c2", "c2", 1000, 1000, 1000, 800);
    frame2->Draw();
    c2->Draw();
    c2->SaveAs("./likelihood_beta.png", "png");

    TCanvas*c3 = new TCanvas("c3", "c3", 1000, 1000, 1000, 800);
    frame3->Draw();
    c3->Draw();
    c3->SaveAs("./likelihood_gamma.png", "png");


    //plotting Minuit contour plot.
    RooPlot* frame_cont = m.contour(gamma, alpha, 1, 2, 3);
    
    TCanvas* ct = new TCanvas("ct", "ct", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont->GetYaxis()->SetTitleOffset(2);
    frame_cont->GetYaxis()->SetTitle("#gamma");
    frame_cont->GetXaxis()->SetTitle("#beta");
    frame_cont->Draw();
    ct->Draw();
    ct->SaveAs("./cont_ab.png", "png");


    //Creating many montecarlo events with this setup
    RooMCStudy mcs(*genpdf, RooArgSet(t,p));
    mcs.generateAndFit(1,50000) ;  // generate & fit 100 samples of 500 events each

    TCanvas* mc_can = new TCanvas("mc_can", "mc_can", 1000,1000,1000,800);
    RooPlot* mc_frame1 =  mcs.plotParam(alpha) ; // Draw distribution of parameter
    mc_frame1->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./alpha_dist.png", "png");

    RooPlot* mc_frame2 =  mcs.plotError(alpha) ; // Draw distribution of parameter error
    mc_frame2->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./alpha_error_dist.png", "png");

    RooPlot* mc_frame3 =  mcs.plotPull(alpha) ; // Draw distribution of parameter pull
    mc_frame3->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./alpha_pull.png", "png");

    //#######BINNING DATA AND CHI2 FIT#############

    t.setBins(100);
    p.setBins(100);
    RooDataHist data_hist("data_hist","binned version of data",RooArgSet(t,p),*data) ;
    RooChi2Var chi2("chi2","chi2",*genpdf,data_hist,DataError(RooAbsData::SumW2)) ;
    RooMinimizer n(chi2) ;
    n.migrad();

    auto alpha_est_chi = alpha.getVal();
    auto beta_est_chi = beta.getVal();
    auto gamma_est_chi = gamma.getVal();
    auto scale_est_chi = scale.getVal();

    //defining a TF2 object with estimated chi2 parameters
    const int npar_chi = 4;
    double dist_params_chi[npar_chi] = {alpha_est_chi, beta_est_chi, gamma_est_chi, scale_est_chi};
    TF2* myPdf_chi = new TF2("myPdf_chi", myPDF, 0, 2*M_PI, 0, 2*M_PI, npar_chi);
    myPdf_chi->SetParameters(dist_params_chi);
    myPdf_chi->SetParNames("#alpha", "#beta", "#gamma", "#scaling");

    //defining TH2F to plot 
    TH2F* dh2_chi = new TH2F("h2_chi", "h2_chi", 100, 0, 2*M_PI, 100, 0, 2*M_PI);
    data->fillHistogram(dh2_chi, RooArgList(t,p));

    //plotting chi2 results
    TCanvas* from_chi2 = new TCanvas("from_ML", "from_ML", 1000,1000,1000,800);
    dh2_chi->SetTitle(Form("Binned #chi^{2} fit: #alpha = %.2f, #beta = %.2f, #gamma = %.2f",alpha_est_chi, beta_est_chi, gamma_est_chi));
    dh2_chi->GetXaxis()->SetTitle("#theta [rad]");
    dh2_chi->GetYaxis()->SetTitle("#phi [rad]");
    dh2_chi->Fit(myPdf_chi);
    dh2_chi->Draw("lego");
    myPdf->Draw("same cont1");
    from_chi2->Draw();
    from_chi2->SaveAs("./result_CHI2_binned.pdf", "pdf");
    


    genpdf->chi2FitTo(data_hist) ;








//Plotting contours
    /*
    RooPlot* frame_cont3 = n.contour(gamma, beta, 1, 2, 3);
    RooPlot* frame_cont2 = m.contour(alpha, beta, 1, 2, 3);

    frame_cont->SetTitle("Minuit contour plot #alpha #gamma");
    frame_cont2->SetTitle("Minuit contour plot #alpha #beta");
    frame_cont3->SetTitle("Minuit contour plot #gamma #beta");

    TCanvas* ct = new TCanvas("ct", "ct", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont->GetYaxis()->SetTitleOffset(2);
    frame_cont->GetYaxis()->SetTitle("#gamma");
    frame_cont->GetXaxis()->SetTitle("#alpha");
    frame_cont->Draw();
    ct->Draw();
    ct->SaveAs("./cont_ag.png", "png");

    TCanvas* ct2 = new TCanvas("ct2", "ct2", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont2->GetYaxis()->SetTitleOffset(2);
    frame_cont2->GetYaxis()->SetTitle("#beta");
    frame_cont2->GetXaxis()->SetTitle("#alpha");
    frame_cont2->Draw();
    ct2->Draw();
    ct2->SaveAs("./cont_ab.png", "png");

    TCanvas* ct3 = new TCanvas("ct3", "ct3", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont3->GetYaxis()->SetTitleOffset(2);
    frame_cont3->GetYaxis()->SetTitle("#beta");
    frame_cont3->GetXaxis()->SetTitle("#alpha");
    frame_cont3->Draw();
    ct3->Draw();
    ct3->SaveAs("./cont_bg.png", "png");
    
    */
    std::cout<<"Over" << "\n";

    
    return 0;

}
