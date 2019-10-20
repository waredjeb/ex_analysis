#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "prng.h"
#include "xorshiro.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"



//funzione gaussiana 2d
double gaussian2d(double x, double y, double mux,double muy, double sx, double sy,double rho){


    double t_rho = 1 - rho*rho;
    double norm = 1./(2*M_PI * sx*sy*sqrt(t_rho));
    double t_mux = (x - mux)/sx;
    double t_muy = (y - muy)/sy;
    double exponent = (-1./(2*t_rho))*(t_mux*t_mux + t_muy*t_muy - 2*rho*(t_mux*t_muy));

    return norm*exp(exponent);

    }

//genero numeri casuali con algoritmo xorshiro128p
xorshiro random_number; //Nota istanza da fare fuori dalla funzione!

double rand_uni(double min = 0., double max = 1.){
    
    return random_number.rand(min,max);

}
//struct punti x e y
struct xy{
    float x;
    float y;
};

//genero coppie di numeri che seguono distribuzione gauss2d con try and catch
//Nota ritorno una struct!!! La funzione deve essere la struct, non un tipo!
xy gauss_tc(double mux,double muy, double sx, double sy,double rho,
                double xmin, double xmax, double ymin, double ymax){
    
    xy rndxy; 
    //double rndx;
    //double rndy;
    
    int i = 0;
    while(1){

        //rndx = rand_uni(xmin,xmax);
        //rndy = rand_uni(ymin,ymax);
        rndxy.x = rand_uni(xmin,xmax);
        rndxy.y = rand_uni(ymin,ymax);

        float rndz = rand_uni();
        //tc
        float rndGauss = gaussian2d(rndxy.x,rndxy.y,mux,muy,sx,sy,rho);


        //std::cout << "fxy:" << rndGauss << "\t rndz: " << rndz << std::endl;
        //std::cout << "Sono qui" << std::endl;
        if(rndGauss > rndz){ //try and catch

            return rndxy; //c++ non permette ritorno di più funzioni, mi serve una struct
        }        
        
    }
    
}



int main(int argc, char** argv){

    //Signal gaussian params
    double mux_sig = 0.;
    double muy_sig = 0.;
    double sx_sig = 0.3;
    double sy_sig = 0.3;
    double rho_sig = 0.5;

    double xmin = -10;
    double xmax = 10;
    double ymin = -10;
    double ymax = 10;

    int N_sig = 10000;

    TCanvas* c1 = new TCanvas("c1","c1",700,700);
    TH2F* gauss_sig = new TH2F("gauss_sig","gauss_sig", 100, xmin, xmax,100,ymin,ymax);
    TH2F* gauss_bkg = new TH2F("gauss_bkg","gauss_bkg",100,xmin,xmax,100,ymin,ymax);
    
    xy xy_sig; //Struct per il segnale
    xy xy_bkg; //Struct per il fondo


    float xentry_sig;
    float yentry_sig;
    float xentry_bkg;
    float yentry_bkg;

    //TMVA Legge dati dai tree, ne creo uno con dentro le variabili x e y per segnale e fondo.
    TTree* tree_sig = new TTree("Signal","Signal");
    tree_sig->Branch("x", &xentry_sig);
    tree_sig->Branch("y", &yentry_sig);

    TTree* tree_bkg = new TTree("Background","Background");
    tree_bkg->Branch("x",&xentry_bkg);
    tree_bkg->Branch("y",&yentry_bkg);

    TFile* outfile = new TFile("dataset.root", "RECREATE");


    for(int i = 0; i < N_sig; i++){

        xy_sig = gauss_tc(mux_sig,muy_sig,sx_sig,sy_sig,rho_sig,xmin,xmax,ymin,ymax);

        xentry_sig = xy_sig.x;
        yentry_sig = xy_sig.y;

        gauss_sig->Fill(xy_sig.x,xy_sig.y);
        tree_sig->Fill();
    }

    //Background gaussian params
    double mux_bkg = 4.;
    double muy_bkg = 4.;
    double sx_bkg = 1;
    double sy_bkg = 1;
    double rho_bkg = 0.4;

    int N_bkg = 10000;  
   

    for(int i = 0; i < N_bkg; i++){

        xy_bkg = gauss_tc(mux_bkg,muy_bkg,sx_bkg,sy_bkg,rho_bkg,xmin,xmax,ymin,ymax);

        xentry_bkg = xy_bkg.x;
        yentry_bkg = xy_bkg.y;

        gauss_bkg->Fill(xy_bkg.x,xy_bkg.y);
        tree_bkg->Fill();
    }
    
    gauss_sig->SetMarkerColor(kRed);
    gauss_bkg->SetMarkerColor(kBlue);
    gauss_sig->SetStats(0000);
    gauss_sig->Draw();
    gauss_bkg->Draw("SAME");
    c1->Print("./analysis/dataset.png","png");
    
    tree_sig->Write();
    tree_bkg->Write();
    outfile->Close();

 return 0;   
}
