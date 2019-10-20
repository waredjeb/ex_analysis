//c++ -o prova prova.cpp -L/usr/local/lib -lgsl -lgslcblas -lm
//c++ prng.cpp xorshiro.cpp xorshiroGen.cpp  prova.cpp -L/usr/local/lib -lgsl -lgslcblas -lm `root-config --libs --glibs --cflags` -o prova
#include <cmath>
#include <iostream>
#include "xorshiroGen.h"
#include <algorithm>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "simulation.h"
#include <string>
#include <sstream>

void line(){

    std::cout << "######################" << std::endl;
}

int main(int argc, char** argv){
 
//Simulation parameters
    int N = 50000;
    double diameter = 0.75e-6;
    double lambda = 400e-9;
    double c = 0.03;
    double I0 = 1;
    double ka = (2*M_PI/lambda*diameter);

    //Histo parameter
    double xmax = ka;
    double xmin = -1*ka;
    double ymax = 10;
    double ymin = 0;
    double bins = sqrt(N);

    double std = c * (xmax - xmin);

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    xorshiroGen rndxor;
    xorshiro uniform;

    simulation experiment(N,c,lambda,diameter,I0);

    TH1D* h = new TH1D("histo","histo",bins,xmin,xmax);
    TH1D* h_sm = new TH1D("histo_sm","histo_sm",bins,xmin,xmax);
    TCanvas* c1 = new TCanvas("c1","c1",1000,1000);

    double inf_x = -M_PI/2;
    double max_x = M_PI/2;
    double* clean = experiment.diffraction_tc(inf_x,max_x,ymin,ymax);   
    double* smeared = experiment.diffraction_smeared(inf_x,max_x,bins); 
    //double* diffraction_smeared(long int N,int bins,double xmin, double xmax,double sm_fact);

    RooUnfoldResponse response(bins,xmin,xmax); //Unfolding

    for(int i = 0; i < N; i++){

        double x = clean[i]; 
        
        h->Fill(ka*sin(clean[i]));
        h_sm->Fill(ka*sin(smeared[i]));
        response.Fill(ka*sin(smeared[i]),ka*sin(clean[i]));     
    }
    //last par, iterations
    RooUnfoldBayes unfold(&response, h, 3);
    RooUnfoldBinByBin unfold_bbb(&response,h);
    RooUnfoldSvd unfold_svd(&response,h,5);

    TH1D* h_unfold =(TH1D*)unfold.Hreco();
    TH1D* h_unfold_bbb = (TH1D*)unfold_bbb.Hreco();
    TH1D* h_unfold_svd = (TH1D*)unfold_svd.Hreco();
    
   TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);

   leg->AddEntry(h, "Ideal-data");
   leg->AddEntry(h_sm, "Smeared-data");
   leg->AddEntry(h_unfold, "Bayes unfold");
   leg->AddEntry(h_unfold_bbb, "BinByBin unfold");
   leg->AddEntry(h_unfold_svd, "SVD unfold");
   //Saving the results
   h->SetLineColor(kBlue);
   h_sm->SetLineColor(kRed);
   //h_sm->GetXaxis()->SetRangeUser(0,5);
  
   h->SetTitle("");
   h->SetStats(0000);
   h_unfold->GetYaxis()->SetTitle("I(#theta)");
   h_unfold->GetYaxis()->SetTitleOffset(1.1);
   h_unfold->GetXaxis()->SetTitle("kasin(#theta)");
   
   
   //h->Scale(I0/h->GetMaximum());
   //h_sm->Scale(I0/h_sm->GetMaximum());

   //h->DrawNormalized("",I0);
   //h_sm->DrawNormalized("SAME",I0);
   //h->Write("clean");
   
   h_unfold->SetLineColor(kGreen);   
   h_unfold->DrawNormalized("");
   h->DrawNormalized("same");
   h_sm->DrawNormalized("same");

   h_unfold_bbb->SetLineColor(kPink);
   h_unfold_bbb->DrawNormalized("SAME");

   h_unfold_svd->SetLineColor(kOrange);
   h_unfold_svd->DrawNormalized("SAME");
   //h_sm->Write("smeared");
   leg->Draw("SAME");
   c1->Draw("");

    std::ostringstream string_temp;
    string_temp << "_" << std << "-" 
        << lambda << "-" << diameter;
    std::string name_file = string_temp.str();

   std::cout << "\n\n Results simulation" << std::endl;
   line();
   line();
   std::cout << "Smearing factor: " << std << "\n" << "Lambda: " << lambda << "\n" << "Diameter: " << diameter << std::endl; 
   line();
   line();
   //c1->SaveAs(("/home/wahid/Wahid/University/analisi_statistica/esercizi/es6/results/unfold"+name_file+".root").c_str(),"root");
   c1->SaveAs("/home/wahid/Wahid/University/analisi_statistica/esercizi/es6/results/unfold_bayes_bbb_svd.pdf","pdf");
   //data->Close();
   return 0;
     

     


}



 
