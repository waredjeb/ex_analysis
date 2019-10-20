void analysis(){

    std::fstream uni_genFile;
    uni_genFile.open("../data/un_xorshiro128p.txt");

    TCanvas * c1 = new TCanvas("xorshiro128","uni_gen", 700,700);

    double min = 0;
    double max = 1.;
    int bins = 100;

    TH1D * uni_Histo = new TH1D("uni_gen","uni_gen", bins,min,max);
    
    double entry;

    while(uni_genFile >> entry){
        uni_Histo -> Fill(entry);
    }

    uni_Histo->GetYaxis()->SetRangeUser(0,1400);
    uni_Histo->SetTitle("Uniform distribution - xorshiro1280plus");
    uni_Histo->GetXaxis()->SetTitle("x");
    uni_Histo->GetYaxis()->SetTitleOffset(1.5);
    uni_Histo->GetYaxis()->SetTitle("counts");
    uni_Histo->SetLineWidth(1);
    //uni_Histo->SetFillStyle(3001);
    uni_Histo -> Draw();
    
    uni_Histo->SetFillColorAlpha(kTeal+10, 0.4);
    c1->Print("un_xorshiro128p.pdf","pdf");
}
