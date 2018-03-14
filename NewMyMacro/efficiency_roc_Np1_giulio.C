#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TFile.h"
#include "TROOT.h"
#include <TLatex.h>



void setHistoStyle(TH1F *histo, std::string title_name, std::string Yaxis_name) {
    
    histo->GetYaxis()->SetTitleOffset(.8);
    histo->GetXaxis()->SetTitleOffset(0.5);
    histo->SetStats(0);
    histo->SetTitleFont(42,"x");
    histo->SetTitleFont(42,"y");
    histo->SetTitleSize(0.05, "XYZ");
    histo->SetTitle(title_name.c_str());
    histo->SetYTitle(Yaxis_name.c_str());
    histo->SetXTitle("");

            
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->GetXaxis()->SetLabelSize(0.05);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleOffset(0.8);
}



void computeSensitivityAndError(TH1D* hist_BDT_S, TH1D* hist_BDT_B, TH1D* hist_sensitivity, int seed_idx, Double_t totalSensitivity_thisVariable[], Double_t totalSensitivityError_thisVariable[]) {
    
    float totalSensitivitySquared_thisVariable = 0;
    float totalSensitivityErrorSquared_thisVariable = 0;
    
                for (int n = 1; n <= hist_BDT_S->GetXaxis()->GetNbins(); n++) {

                
                    double SignalScan = hist_BDT_S->GetBinContent(n);
                    double BackgroundScan = hist_BDT_B->GetBinContent(n);
                    double Signal_RelativeError = hist_BDT_S->GetBinError(n)/SignalScan;
                    double BKG_RelativeError = hist_BDT_B->GetBinError(n)/BackgroundScan;

                    totalSensitivitySquared_thisVariable +=  ((BackgroundScan > 0.0001) && (SignalScan > 0.0001)) ? SignalScan*SignalScan/(BackgroundScan) : 0.;
                    hist_sensitivity->SetBinContent(n, ((BackgroundScan > 0.0001) && (SignalScan > 0.0001)) ? SignalScan/sqrt(BackgroundScan) : 0.); 
                
                    float SensitivityRelativeErrorSquared = (BackgroundScan > 0.00001 && SignalScan > 0.00001) ? (4*Signal_RelativeError*Signal_RelativeError + BKG_RelativeError*BKG_RelativeError) : 0.;
                    totalSensitivityErrorSquared_thisVariable += SensitivityRelativeErrorSquared*hist_sensitivity->GetBinContent(n)*hist_sensitivity->GetBinContent(n);
                    
                    std::cout << "totalSensitivityErrorSquared   "  << hist_BDT_B->GetBinError(n) <<  "  \t  "  << BackgroundScan <<  "  \t  "  << BKG_RelativeError <<  "  \t  "  << SensitivityRelativeErrorSquared <<  "  \t  "  <<hist_sensitivity->GetBinContent(n)*hist_sensitivity->GetBinContent(n) <<  "  \t  "  << totalSensitivityErrorSquared_thisVariable << std::endl;
                    
                }
                
                totalSensitivity_thisVariable[seed_idx]      = sqrt(totalSensitivitySquared_thisVariable);
                totalSensitivityError_thisVariable[seed_idx] = sqrt(totalSensitivityErrorSquared_thisVariable);
                
    
    
}




void  DrawCorrelationMatrix(TH2F * CorrelationMatrixS, TH2F * CorrelationMatrixB, std::string fileName) {
    
    CorrelationMatrixS->SetTitle("signal");
    CorrelationMatrixB->SetTitle("background");
    
    TCanvas * canv = new TCanvas("canvas", "", 1200, 600);
    canv->Divide(2,1);
    
    canv->cd(1);
    CorrelationMatrixS->Draw("colz, text");
    
    canv->cd(2);
    CorrelationMatrixB->Draw("colz, text");
    
    canv->Print(("figure/efficiency_figure/CorrelationMatrix/CorrelationMatrix_"+fileName+".png").c_str()); 
}



void DrawAndSave_bkg_and_signal(TH1D* hist_BDT_S, TH1D* hist_BDT_B, std::string Sample) {
    
    
        TCanvas * canv = new TCanvas("canvas", "", 800, 600);    

        gStyle->SetOptStat(0000);

                
        std::string nameFig_Super = hist_BDT_S->GetTitle();
        std::string nameFig_S     = hist_BDT_S->GetTitle();
        std::string nameFig_B     = hist_BDT_B->GetTitle();
        nameFig_Super = "figure/efficiency_figure/BDToutput/sensitivity"+Sample+"/"+nameFig_S+"_"+Sample+"_SuperImposed";
        nameFig_S = "figure/efficiency_figure/BDToutput/singleDistribution/"+Sample+"/"+nameFig_S+"_"+Sample+"_S";
        nameFig_B = "figure/efficiency_figure/BDToutput/singleDistribution/"+Sample+"/"+nameFig_B+"_"+Sample+"_B";
            
        gPad->SetLogy();
        
        hist_BDT_S->SetLineWidth(2);
        hist_BDT_B->SetLineWidth(2);
        
        hist_BDT_S->SetLineColor(kRed);
        hist_BDT_S->SetFillColorAlpha(kRed, 0.2);
        hist_BDT_S->SetFillStyle(3004);


        hist_BDT_B->SetLineColor(kBlue);
        hist_BDT_B->SetFillColorAlpha(kBlue, 0.2);
        hist_BDT_B->SetFillStyle(3005);
                
                
        hist_BDT_S->Draw("hist");
        canv->Print((nameFig_S+".png").c_str());
        
        
        hist_BDT_B->Draw("hist");
        canv->Print((nameFig_B+".png").c_str());   
        
        
        TLegend *leg = new TLegend(0.55,0.7,0.89,0.85);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->AddEntry(hist_BDT_S,"Signal","f");
        leg->AddEntry(hist_BDT_B,"Background","f");
        
        hist_BDT_S->GetXaxis()->SetTitle("tanh^{-1}((BDT response + 1)/2)");
        hist_BDT_S->GetYaxis()->SetTitle("Entries");  
        
        hist_BDT_S->GetYaxis()->SetRangeUser(0.01, 100000); 
//         hist_BDT_B->Write();
//         hist_BDT_S->Write();
        hist_BDT_S->Draw("hist");
        hist_BDT_B->Draw("hist same");
        leg->Draw("same");

        canv->Print((nameFig_Super+".png").c_str());        
    
}
            
   
   
   
   
void DrawSensitivity(TH1D* hist_sensitivity, float sensitivityToPrint, std::string Sample) {
   
    
        std::string nameFig     = hist_sensitivity->GetTitle();
        nameFig = "figure/efficiency_figure/BDToutput/sensitivity"+Sample+"/"+nameFig+"_Sensitivity";
        
        


        float totalSensitivity = 0.;
        for (int n = 1; n <= hist_sensitivity->GetXaxis()->GetNbins(); ++n) {
            totalSensitivity += hist_sensitivity->GetBinContent(n)*hist_sensitivity->GetBinContent(n);
//             std::cout << "bin " << n << " sensitivity   " << hist_sensitivity->GetBinContent(n) << " \t Sigma:   " << totalSensitivity <<"signal  "<< hist_BDT_S->GetBinContent(n)<<" "<<"bkg "<<hist_BDT_B->GetBinContent(n)<< std::endl;
        }

        if(abs(sqrt(totalSensitivity) - sensitivityToPrint) > 0.001)  std::cout << "ERROR: TWO DIFFERENT SENSIBILITY" << std::endl;
        std::cout << "sensitivityToPrint    "  << sensitivityToPrint << std::endl;
        
        std::ostringstream sigmaString;
        sigmaString<<std::setprecision(3)<<sensitivityToPrint;
        TLatex* tex = new TLatex(0.50,0.80,("Tot sigma = " + sigmaString.str() ).c_str());
        tex->SetNDC();
        tex->SetTextAlign(35);
        tex->SetTextFont(42);
        tex->SetTextSize(0.06);
        tex->SetLineWidth(2);




        hist_sensitivity->SetLineColor(kOrange+4);
        hist_sensitivity->SetLineWidth(2);
        gPad->SetGridy();
        
        hist_sensitivity->GetXaxis()->SetTitle("tanh^{-1}((BDT response + 1)/2)");
        hist_sensitivity->GetYaxis()->SetTitle("bin per bin sensitivity");  

        
        hist_sensitivity->SetLineWidth(3);

//         double titleOffset  = 1.2;
//         double titleSize    = 0.06;
//         double labelOffset  = 0.01;
//         double labelSize    = 0.04;
// 
//         hist_sensitivity->GetYaxis()->SetTitleOffset(titleOffset);
//         hist_sensitivity->GetXaxis()->SetTitleOffset(titleOffset);
//         hist_sensitivity->GetYaxis()->SetLabelOffset(labelOffset);
//         hist_sensitivity->GetXaxis()->SetLabelOffset(labelOffset);
//         hist_sensitivity->GetYaxis()->SetTitleSize(titleSize);
//         hist_sensitivity->GetXaxis()->SetTitleSize(titleSize);
//         hist_sensitivity->GetYaxis()->SetLabelSize(labelSize*1.5);
//         hist_sensitivity->GetXaxis()->SetLabelSize(labelSize);
        
        TCanvas *canv=new TCanvas("c","",800, 600); 
        canv->cd();
        hist_sensitivity->Draw();
        tex->Draw();
        canv->Print((nameFig+".png").c_str());
        
        
        
}
     
     
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////  FLOAT --> STRING //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string makeString(float x, int precision = 7) {
    
    std::ostringstream x_s;
    x_s << std::setprecision(precision) << x; 
    std::string x_s_string(x_s.str());
    return x_s_string;
}
    

    

void drawTestAndTrainSuperimposed(TH1D* htest,TH1D* htrain,std::string figname, float kolmog, float prob_chiq){
    TCanvas *c=new TCanvas("c","",1000,800);
    c->cd();
    gPad->SetLogy();
    
    
    htest->SetLineWidth(2);
    htrain->SetLineWidth(2);

    htest->SetLineColor(2);
    htest->Draw("E0");
    htrain->Draw("same,E0");
    
    float Xbottom = 0.2;
    float Ybottom = 0.2;
    float Xtop = 0.3;
    float Ytop = 0.3;
    TLegend *leg = new TLegend(Xbottom,Ybottom,Xtop,Ytop);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);

    leg->AddEntry(htest,"test","l");
    leg->AddEntry(htrain,"train","l");
    leg->Draw();
            
    
    std::string kolmog_string         = makeString(kolmog);
    std::string prob_chiq_string      = makeString(prob_chiq);

            
            
    float yt2 = 0.5;
    float xt2 = 0.7;
    TLatex *t2 =new TLatex(xt2,yt2,("KS= "+kolmog_string).c_str());
    t2->SetNDC();
    t2->SetTextAlign(22);
    t2->SetTextSize(0.05);
    t2->Draw();
            

    TLatex *t3 =new TLatex(xt2,yt2+0.2,( "Prob(Chi2)="+prob_chiq_string).c_str());
    t3->SetNDC();
    t3->SetTextAlign(22);
    t3->SetTextSize(0.05);
    t3->Draw();
    c->Print(("figure/efficiency_figure/BDToutput/TestAndTrainComparison/"+figname).c_str());
}



    
    
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////  INTEGRA I BIN CON MENO DI 0.5     NON USATA //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void integrate_when_bkg_is_less_than_05(TH1D * hist_BDT_S, TH1D * hist_BDT_B) {
    
    
        bool oltreBKG = false;
        int lastGood=0;
        for(int m=1; m <= 40;++m){
            if( hist_BDT_B->GetBinContent(m)<0.5 && oltreBKG == false) {
                oltreBKG = true;
                lastGood=m-1;
                cout<<"m="<<m<<endl;
            }
            if(oltreBKG){
                hist_BDT_S->SetBinContent(lastGood,hist_BDT_S->GetBinContent(lastGood)+ hist_BDT_S->GetBinContent(m));
                hist_BDT_S->SetBinContent(m,0);
                hist_BDT_B->SetBinContent(lastGood,hist_BDT_B->GetBinContent(lastGood)+hist_BDT_B->GetBinContent(m));
                hist_BDT_B->SetBinContent(m,0);
            }
        } 
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////  DISEGNA IL TEST N+1 E LE PROBABILITA /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void drawGraph(TCanvas * canv, TH1F *frame, TGraph *gr, std::string fig_name, bool isSensitivity, int n_variables, float significance_nomore) {

        canv->cd();
        canv->SetBottomMargin(.3);
        canv->SetLeftMargin(.1);
        gPad->SetLogy(0);
        gPad->SetGridx();
        gPad->SetGridy();
        frame->Draw();
                
        
        gr->SetMarkerStyle(21);
        gr->SetLineWidth(2);
        gr->Draw("PL");
        
        
        
        if (isSensitivity) {
            
        TLine *line2 = new TLine(0., significance_nomore ,(Double_t) n_variables, significance_nomore);
        line2->SetLineStyle(2);
        line2->SetLineColor(2);
        line2->SetLineWidth(2);
        line2->Draw("Lsame");
        }
        
        
        gr->SaveAs(("figure/efficiency_figure/"+fig_name+".root").c_str());
        canv->Print(("figure/efficiency_figure/"+fig_name+".png").c_str());
                
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////  TROVA IL BIN DOWN LIMIT //////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int  FindBinDown(TH1D *hist_BDT_binning_B, int binLimitUp, int minNumberOfEventPerBin, int binMinNumber) {
        
    int binLimitDown = 0.;
        for(int n = binLimitUp-binMinNumber; n > 0; n--) {
            if (hist_BDT_binning_B->Integral(n+1, binLimitUp) >= minNumberOfEventPerBin ) {
                binLimitDown = n;
                break;
            }
             
        }
        
    return binLimitDown;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////  TROVA IL BINNING DEGLI ISTOGRAMMI ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void findBinning(float  binning_BDT[], int Nbins, int Nbins_binning, int max, int minNumberOfEventPerBin, int binMinNumber, TH1D * hist_BDT_binning_B) {

    
    for(int n = 0; n<Nbins; n++)  binning_BDT[n]=(n-Nbins)*0.000001;   
                
                
    int binLimitUp   = Nbins_binning;
    int binLimitDown = Nbins_binning;

    
    
    int binLimitDown10= Nbins_binning;
    int binLimitDown11= Nbins_binning;



    for(int n = Nbins-1; n > 0 && binLimitDown>0; n--) {
        binning_BDT[n] = (1.*binLimitDown*max)/Nbins_binning;//binLimitDown*max/Nbins_binning;
        binLimitUp = binLimitDown;
        binLimitDown10 = FindBinDown(hist_BDT_binning_B, binLimitUp, minNumberOfEventPerBin, binMinNumber);
        binLimitDown11 = FindBinDown(hist_BDT_binning_B, binLimitUp, minNumberOfEventPerBin+1, binMinNumber);
        binLimitDown = (binLimitDown10+binLimitDown11)/2;
                    
    }
                    
                    
//     for(int n = Nbins-1; n > 0 && binLimitDown>0; n--) {
// //         binning_BDT[n] = (1.*binLimitDown*max)/Nbins_binning;
//         binning_BDT[n] = binLimitDown;
//         binLimitUp = binLimitDown;
//         binLimitDown = FindBinDown(hist_BDT_binning_B, binLimitUp, minNumberOfEventPerBin, binMinNumber);
//     }
    binning_BDT[0] = 0.;
    
    for(int n = 0; n < Nbins; n++) std::cout << binning_BDT[n] << ", \t";
    std::cout  << std::endl;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////  FUNZIONI PER LA DISTRIBUZIONE DELLE SENSIBILITA CON DIVERSI SEED ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1F * makeSensitivityhisto(std::string histoName, Double_t totalSensitivity[], int numberOfFile, float limitUp, float limitDown, int numberOfBin) {
    
    TH1F * histo = new TH1F (histoName.c_str(), histoName.c_str(), numberOfBin, limitUp, limitDown);   
    
    for(int n = 0; n < numberOfFile; n++)  histo->Fill(totalSensitivity[n]);
    
    setHistoStyle(histo, "sensitivity distribution", "Entries");
    histo->SetLineWidth(2);
    histo->SetLineColor(4);
    return histo;
}








float findMean(Double_t totalSensitivity[], int numberOfFile) {
    float mean=0;
    for(int n = 0; n < numberOfFile; n++)  mean += totalSensitivity[n];
    return mean/numberOfFile;
}

float findRMS(Double_t totalSensitivity[], int numberOfFile) {
    float RMS=0;
    float mean=findMean(totalSensitivity, numberOfFile);
    
    for(int n = 0; n < numberOfFile; n++)  RMS += (totalSensitivity[n] - mean)*(totalSensitivity[n] - mean);
    return sqrt(RMS/(numberOfFile-1)) ;
}




void drawSinglePlot(Double_t arrayToDraw[], int numberOfFile, std::string histoTitle, std::string figureName, float limitUp, float limitDown, int numberOfBin) {
    
    TH1F * histoDiff = makeSensitivityhisto(("histoDiff"+figureName).c_str(), arrayToDraw, numberOfFile, limitUp, limitDown, numberOfBin);
    histoDiff->SetTitle(histoTitle.c_str());

    float mean = findMean(arrayToDraw, numberOfFile);
    float RMS  = findRMS(arrayToDraw, numberOfFile);

    
    TCanvas *canvDiff = new TCanvas(("canvDiff"+figureName).c_str(), ("canvDiff"+figureName).c_str(), 800, 600);
    canvDiff->cd();
    
 
    TLatex* tex = new TLatex(0.40,0.80,("mean = " + makeString(mean, 3) + " +- " + makeString(RMS/sqrt(numberOfFile), 2)  ).c_str());
    tex->SetNDC();
    tex->SetTextAlign(35);
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    

    
    histoDiff->Draw();
    tex->Draw();
    canvDiff->Print(("figure/efficiency_figure/" + figureName +".png").c_str());
    
}






void drawAndSaveSensitivities(Double_t totalSensitivity[], Double_t totalSensitivity_train[], int numberOfFile, std::string legend, std::string legend_train, std::string histoTitle, std::string figureName, float limitUp, float limitDown, int numberOfBin) {


    
    TH1F * histo = makeSensitivityhisto("histo", totalSensitivity, numberOfFile, limitUp, limitDown, numberOfBin);
    TH1F * histo_train = makeSensitivityhisto("histo_train", totalSensitivity_train, numberOfFile, limitUp, limitDown, numberOfBin);
    histo->SetTitle(histoTitle.c_str());
    histo_train->SetTitle(histoTitle.c_str());
    
    histo_train->SetLineStyle(2);
    
    float mean_test  = findMean(totalSensitivity, numberOfFile);
    float mean_train = findMean(totalSensitivity_train, numberOfFile);
    
    float RMS_test   = findRMS(totalSensitivity, numberOfFile);
    float RMS_train  = findRMS(totalSensitivity_train, numberOfFile);
    
    
    float Xbottom = 0.1;
    float Ybottom = 0.7;
    float Xtop = 0.5;
    float Ytop = 0.9;
    TLegend *myLegend=new TLegend(Xbottom, Ybottom, Xtop, Ytop, "");
    myLegend->SetBorderSize(0);
    myLegend->SetFillColor(0);
    myLegend->SetFillStyle(0);
    myLegend->SetTextFont(42);
    myLegend->SetTextSize(0.04);
    myLegend->SetBorderSize(0);                  // without border
    myLegend->SetFillColorAlpha(1,0);       
    myLegend->AddEntry(histo,(legend + "mean = " + makeString(mean_test, 3) + " +- " + makeString(RMS_test, 2)).c_str(),"l");    
    myLegend->AddEntry(histo_train,(legend_train + "mean = " + makeString(mean_train, 3) + " +- " + makeString(RMS_train, 2)).c_str(),"l");    
    
    
    TCanvas *canvS = new TCanvas("canvS", "canvS", 800, 600);
    canvS->cd();
    histo->Draw();
    histo_train->Draw("same");
    myLegend->Draw();
    canvS->Print(("figure/efficiency_figure/testNp1/" + figureName +".png").c_str());
    

    
}







void efficiency_roc_Np1_giulio(){
//	gROOT->ProcessLine(".x setTDRStyle.C");

            TCanvas *canv = new TCanvas("canv", "canv", 1400, 700);
            canv->SetBottomMargin(.15);
            canv->SetLeftMargin(.15);
            gPad->SetLogy();
            


//                 const int n_variables =38;
//         std::string variables_names[n_variables]={"Inv_mass", "energytot", "W_mass_virtual1", "W_mass_virtual2", "qgl_1q", "qgl_2q", "thetastarW2", "thetastarW1", "theta1", "qq_pt", "theta2", "W_Pt_virtual1", "W_Pt_virtual2", "ll_eta", "EWKHTsoft", "DeltaEtaQQ", "diffMassWWH", "ll_pt", "Jet3_pt", "ll_zstar", "met_pt", "softLeadingJet_pt", "btagCMVA", "cosThetaStarJet", "WWmass", "impulsoZ", "deltaMRel", "randomVariable", "cosThetaPlane", "softActivityEWK_njets2", "softActivityEWK_njets5", "softActivityEWK_njets10", "W_eta_virtual1", "W_eta_virtual2", "E_parton1", "E_parton2", "deltaM", "nomore"};
            
            

            
            const int primary_variables_number=6;
            std::string variables_names_array_primary[primary_variables_number]={"ll_mass","Mqq", "RptHard","ll_zstar", "softActivityEWK_njets10", "ll_pt"};
            
//             const int primary_variables_number=8;
//             std::string variables_names_array_primary[primary_variables_number]={"ll_mass", "Mqq", "RptHard","ll_zstar","softActivityEWK_njets5","ll_pt","W_mass_virtual2","W_Pt_virtual1"};
            
            
            const int max_variables_number=39;
            std::string all_variables_names[max_variables_number]={"Inv_mass","ll_mass","energytot","W_mass_virtual1","W_mass_virtual2","qgl_1q","qgl_2q" ,"thetastarW2","thetastarW1","theta1", "qq_pt","theta2","W_Pt_virtual1","W_Pt_virtual2","Mqq", "RptHard", "ll_eta", "EWKHTsoft", "DeltaEtaQQ" ,"diffMassWWH", "ll_pt","Jet3_pt","ll_zstar","met_pt","softLeadingJet_pt","btagCMVA", "cosThetaStarJet","WWmass","impulsoZ", "deltaMRel", "cosThetaPlane","softActivityEWK_njets2","softActivityEWK_njets5","softActivityEWK_njets10","W_eta_virtual1","W_eta_virtual2","E_parton1","E_parton2","deltaM"};
            

            
            
//             const int max_variables_number=17;
//             std::string all_variables_names[max_variables_number]={"Inv_mass","ll_mass","energytot","W_mass_virtual1","W_mass_virtual2","qgl_1q","qgl_2q" ,"W_Pt_virtual1","W_Pt_virtual2","Mqq", "RptHard", "ll_eta","randomVariable", "cosThetaPlane","softActivityEWK_njets2","softActivityEWK_njets5","deltaM"};
            
//             const int max_variables_number=14;
//             std::string all_variables_names[max_variables_number]={"ll_mass", "Mqq", "RptHard","ll_zstar","softActivityEWK_njets5", "softActivityEWK_njets10","ll_pt","W_mass_virtual2","W_Pt_virtual1", "energytot","met_pt","qgl_1q","impulsoZ","randomVariable"};
            
//             const int max_variables_number=6;
//             std::string all_variables_names[max_variables_number]={"ll_mass","Mqq", "RptHard","ll_zstar", "softActivityEWK_njets10","randomVariable"};
               

//             const int max_variables_number=101;
// //             const int max_variables_number=10;
//             std::string all_variables_names[max_variables_number];
//             all_variables_names[0]="ll_mass";
//             all_variables_names[1]="Mqq";
//             all_variables_names[2]="RptHard";
                
                
                
//             for(int n = 0; n < max_variables_number-3; n++){
//                 
//                 std::ostringstream n_s;
//                 n_s << n+2; 
//                 std::string n_s_string(n_s.str());
//                 all_variables_names[n+3] = "qgl_1q_" + n_s_string + "SeedGen_";
//             }
            
            
//             const int n_variables =max_variables_number - primary_variables_number;
            const int n_variables =max_variables_number - primary_variables_number +1;
            std::string variables_names[n_variables];
            variables_names[n_variables-1] = "nomore";
            

            
            int variable_idx = 0;
            for(int i=0; i<max_variables_number;i++){
                bool damettere=true;
                for(int j=0; j<primary_variables_number;j++){
                    if(all_variables_names[i].compare(variables_names_array_primary[j]) == 0)
                        damettere=false;
                }
                if(damettere==true) {
                variables_names[variable_idx] = all_variables_names[i];
                variable_idx++;
                }
            }
            

        std::string variable_used_in_BDT = "";
        for(int i=0; i<primary_variables_number-1;i++)  
            variable_used_in_BDT = variable_used_in_BDT + variables_names_array_primary[i] + ", ";
        variable_used_in_BDT = variable_used_in_BDT + variables_names_array_primary[primary_variables_number-1];
        

        
	std::string file_names[n_variables];
	Double_t sensitivity_ex[n_variables];
	Double_t sensitivity_ey[n_variables];
        Double_t sensitivity_train_ey[n_variables];
        
	Double_t frame2_axisx[n_variables];
	Double_t frame2_axisx_train[n_variables];
	Double_t Chi_S_axisx[n_variables];
	Double_t Chi_B_axisx[n_variables];
	Double_t totalProbS_axisx[n_variables];
	Double_t totalProbB_axisx[n_variables];
	for (int i=0;i<n_variables;i++){
            sensitivity_ex[i] = 0.;
            sensitivity_ey[i] = 0.3;
            sensitivity_train_ey[i] = 0.3;
            frame2_axisx[i] = 0.5+i;
            frame2_axisx_train[i] = 0.6+i;
            Chi_S_axisx[i] = 0.5+i;
            Chi_B_axisx[i] = 0.5+i;
            totalProbS_axisx[i] = 0.5+i;
            totalProbB_axisx[i] = 0.5+i;
            file_names[i] = variables_names[i];
	}


	std::string end = "muaxis2jet2q";
        Double_t totalSensitivityErrorSquared[n_variables];
	Double_t totalSensitivity_train[n_variables];
	Double_t totalSensitivity[n_variables];
	Double_t totalChiS[n_variables];
        Double_t totalChiB[n_variables];
        Double_t totalProbS[n_variables];
        Double_t totalProbB[n_variables];
        Double_t totalKS_S[n_variables];
        Double_t totalKS_B[n_variables];


        
        std::string trainingOption = "";
//         trainingOption = "nominal";
//         trainingOption = "noVarTransform";
//         trainingOption = "sinVarTransform";
//         trainingOption = "maxDepth";

        std::string nameFileToAdd = "_Trees150_nodeSize7_maxDepth6";
        nameFileToAdd = "";
        
        
        
        for (int i=0;i<n_variables;i++) std::cout << "file_names[" << i << "]  \t" << file_names[i] << std::endl;
// 	for (int i=0;i<n_variables;i++){
//             file_names[i] = "/scratch/mandorli/Hmumu/training/CMSSW_8_0_28/src/training/output/TMVA_main_v25_Np1_"+file_names[i];
//             file_names[i] = file_names[i]+trainingOption+nameFileToAdd+end+".root";
// //             file_names[i] = file_names[i]+trainingOption+end+".root";
// 	}
	
	

        const int Nbins = 10;          //Bins are Nbins -1
        float max = 4.;
        
        int Nbins_binning=10000;
        int minNumberOfEventPerBin = 5;
        float binWidth = 0.15;
        int binMinNumber = (int) ( binWidth*((1.*Nbins_binning)/(1.*max)));
        


        trainingOption = trainingOption + "_"+std::to_string(primary_variables_number)+"To"+std::to_string(primary_variables_number+1)+"_"+std::to_string(minNumberOfEventPerBin)+"ev";    // "_4To5_5ev"
	trainingOption = trainingOption + nameFileToAdd;
        
        
	for (int current_file=0;current_file<n_variables;current_file++){
            
            
                const int numberOfSeedToRead = 19;
            
                Double_t totalSensitivity_thisVariable[numberOfSeedToRead];
                Double_t totalSensitivity_thisVariable_train[numberOfSeedToRead];
                
                Double_t totalSensitivityErrorSquared_thisVariable[numberOfSeedToRead];
                Double_t totalSensitivityErrorSquared_thisVariable_train[numberOfSeedToRead];
                Double_t totalSensitivityError_thisVariable[numberOfSeedToRead];
                Double_t totalSensitivityError_thisVariable_train[numberOfSeedToRead];

                Double_t totalChiS_thisVariable[numberOfSeedToRead];
                Double_t totalChiB_thisVariable[numberOfSeedToRead];
                Double_t totalProbS_thisVariable[numberOfSeedToRead];
                Double_t totalProbB_thisVariable[numberOfSeedToRead];
                Double_t totalKS_S_thisVariable[numberOfSeedToRead];
                Double_t totalKS_B_thisVariable[numberOfSeedToRead];

                
            
            
            
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++){  
                
                totalSensitivity_thisVariable[seed_idx] = 0;
                totalSensitivity_thisVariable_train[seed_idx] = 0;
                
                totalSensitivityErrorSquared_thisVariable[seed_idx] = 0;
                totalSensitivityErrorSquared_thisVariable_train[seed_idx] = 0;
                totalSensitivityError_thisVariable[seed_idx] = 0;
                totalSensitivityError_thisVariable_train[seed_idx] = 0;
                
                totalChiS_thisVariable[seed_idx] = 0;
                totalChiB_thisVariable[seed_idx] = 0;
                totalProbS_thisVariable[seed_idx] = 0;
                totalProbB_thisVariable[seed_idx] = 0;
                totalKS_S_thisVariable[seed_idx] = 0;
                totalKS_B_thisVariable[seed_idx] = 0;
                
                
                int seed = 61 + seed_idx;

                std::ostringstream seed_s;
                seed_s  << seed; 
                std::string seed_string(seed_s.str());
    
                std::string fileToRead = file_names[current_file];
                fileToRead = "/scratch/mandorli/Hmumu/training/CMSSW_8_0_28/src/training/output/TMVA_main_v25_Np1_"+fileToRead+nameFileToAdd+"_"+seed_string+"SeedGen_"+end+".root";
            

            
                TFile *file = new TFile(fileToRead.c_str(), "read");
                TFile *file2= new TFile(fileToRead.c_str(), "read");
                
                TTree *tree     = (TTree*)file->Get("TestTree");
                TTree *train_tree  = (TTree*)file->Get("TrainTree");
                std::cout << "fileToRead  " << fileToRead << std::endl;
                

                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////     CORRELAZIONI      ///////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                if (seed_idx == 0) {
                    TH2F * CorrelationMatrixS = (TH2F*) ((TH2F*)file->Get("CorrelationMatrixS"))->Clone(("CorrelationMatrixS_"+variables_names[current_file]).c_str());
                    TH2F * CorrelationMatrixB = (TH2F*) ((TH2F*)file->Get("CorrelationMatrixB"))->Clone(("CorrelationMatrixB_"+variables_names[current_file]).c_str());
                    DrawCorrelationMatrix(CorrelationMatrixS, CorrelationMatrixB, variables_names[current_file]);
                }
                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////   FINE CORRELAZIONI   ///////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                
                
                
                
                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////    TROVO IL BINNING   ///////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                TH1D *hist_BDT_binning_B        = new TH1D ("hist_BDT_binning_B", variables_names[current_file].c_str(), Nbins_binning, 0., max);
                TH1D *hist_BDT_binning_B_train  = new TH1D ("hist_BDT_binning_B_train", variables_names[current_file].c_str(), Nbins_binning, 0., max);
                
                tree->Draw("atanh((BDTG+1.)/2.)>>hist_BDT_binning_B", "(classID ==1 )");
                train_tree->Draw("atanh((BDTG+1.)/2.)>>hist_BDT_binning_B_train", "(classID ==1 )");


                float binning_BDT[Nbins];
                float binning_BDT_train[Nbins];

                
                
                findBinning(binning_BDT,        Nbins, Nbins_binning, max, minNumberOfEventPerBin, binMinNumber, hist_BDT_binning_B);
                findBinning(binning_BDT_train,  Nbins, Nbins_binning, max, minNumberOfEventPerBin, binMinNumber, hist_BDT_binning_B_train);
                
//                  std::cout << "spero sia uguale" << std::endl;
//                 for(int n = 0; n < Nbins; n++) std::cout << binning_BDT[n] << ", \t";
//                 std::cout  << std::endl;
//                 for(int n = 0; n < Nbins; n++) std::cout << binning_BDT_train[n] << ", \t";
//                 std::cout  << std::endl;
//                 
//                 int binLimitUp   = Nbins_binning;
//                 int binLimitDown = Nbins_binning;
//                 int binLimitUp_train   = Nbins_binning;
//                 int binLimitDown_train = Nbins_binning;
//                 
//                 int binLimitDown10= Nbins_binning;
//                 int binLimitDown11= Nbins_binning;
//                 int binLimitDown10_train= Nbins_binning;
//                 int binLimitDown11_train= Nbins_binning;
// 
// 
// 
//                 for(int n = Nbins-1; n > 0 && binLimitDown>0; n--) {
//                     binning_BDT[n] = binLimitDown*max/Nbins_binning;
//                     binLimitUp = binLimitDown;
//                     binLimitDown10 = FindBinDown(hist_BDT_binning_B, binLimitUp, minNumberOfEventPerBin, binMinNumber);
//                     binLimitDown11 = FindBinDown(hist_BDT_binning_B, binLimitUp, minNumberOfEventPerBin+1, binMinNumber);
//                     binLimitDown = (binLimitDown10+binLimitDown11)/2;
//                     
//                     binning_BDT_train[n] = binLimitDown_train*max/Nbins_binning;
//                     binLimitUp_train = binLimitDown_train;
//                     binLimitDown10_train = FindBinDown(hist_BDT_binning_B_train, binLimitUp_train, minNumberOfEventPerBin, binMinNumber);
//                     binLimitDown11_train = FindBinDown(hist_BDT_binning_B_train, binLimitUp_train, minNumberOfEventPerBin+1, binMinNumber);
//                     binLimitDown_train = (binLimitDown10_train+binLimitDown11_train)/2;
//                     
//                 }
//                 
//                 if(binning_BDT[1] > 0.) binning_BDT[0] = 0.;
//                 if(binning_BDT_train[1] > 0.) binning_BDT_train[0] = 0.;
//                 
//                 for(int n = 0; n < Nbins; n++) std::cout << binning_BDT[n] << ", \t";
//                 std::cout  << std::endl;
//                 for(int n = 0; n < Nbins; n++) std::cout << binning_BDT_train[n] << ", \t";
//                 std::cout  << std::endl;
                

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        

                
                
//                 TH1D *hist_BDT_S = new TH1D (("hist_BDT_"+variables_names[current_file]+"_S").c_str(), (variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_BDT_B = new TH1D (("hist_BDT_"+variables_names[current_file]+"_B").c_str(), (variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_sensitivity = new TH1D ("hist_sensitivity",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_sensitivity_train = new TH1D ("hist_sensitivity_train",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_BDT_S_train_testBinning = new TH1D ("hist_BDT_S_train_testBinning",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_BDT_B_train_testBinning = new TH1D ("hist_BDT_B_train_testBinning",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_BDT_S_train = new TH1D ("hist_BDT_S_train",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_BDT_B_train = new TH1D ("hist_BDT_B_train",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), 40, 0., 4);
//                 TH1D *hist_sum=new TH1D ("hist_sum","Sum_signal+bkg", 40, 0., 4);
                
                TH1D *hist_BDT_S = new TH1D (("hist_BDT_"+variables_names[current_file]+"_S").c_str(), (variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT);
                TH1D *hist_BDT_B = new TH1D (("hist_BDT_"+variables_names[current_file]+"_B").c_str(), (variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT);
                TH1D *hist_BDT_S_train_testBinning = new TH1D ("hist_BDT_S_train_testBinning",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT);
                TH1D *hist_BDT_B_train_testBinning = new TH1D ("hist_BDT_B_train_testBinning",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT);
                TH1D *hist_BDT_S_train = new TH1D ("hist_BDT_S_train",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT_train);
                TH1D *hist_BDT_B_train = new TH1D ("hist_BDT_B_train",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT_train);
                TH1D *hist_sum=new TH1D ("hist_sum","Sum_signal+bkg",Nbins-1, binning_BDT);
                TH1D *hist_sensitivity = new TH1D ("hist_sensitivity",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT);
                TH1D *hist_sensitivity_train = new TH1D ("hist_sensitivity_train",(variables_names[current_file]+"_"+seed_string+"Seed").c_str(), Nbins-1, binning_BDT_train);
                
                hist_BDT_S->Sumw2();
                hist_BDT_B->Sumw2();
                hist_BDT_S_train->Sumw2();
                hist_BDT_B_train->Sumw2();
                


            

                
                tree->Draw(("atanh((BDTG+1.)/2.)>>hist_BDT_"+variables_names[current_file]+"_S").c_str(), "(classID ==0 )");
                tree->Draw(("atanh((BDTG+1.)/2.)>>hist_BDT_"+variables_names[current_file]+"_B").c_str(), "(classID ==1 )");
                train_tree->Draw("atanh((BDTG+1.)/2.)>>hist_BDT_S_train", "(classID ==0 )");
                train_tree->Draw("atanh((BDTG+1.)/2.)>>hist_BDT_B_train", "(classID ==1 )");
                train_tree->Draw("atanh((BDTG+1.)/2.)>>hist_BDT_S_train_testBinning", "(classID ==0 )");
                train_tree->Draw("atanh((BDTG+1.)/2.)>>hist_BDT_B_train_testBinning", "(classID ==1 )");

      
//////////////////////////////////0.5 non significa niente perche sta prima di Scale(). nel codice vecchio stava qui (per errore)////////////////////////////////////////////////////////////////////
/*                bool oltreBKG = false;
                int lastGood=0;
                for(int m=1; m <= 40;++m){
                if( hist_BDT_B->GetBinContent(m)<0.5 && oltreBKG == false) {
                    oltreBKG = true;
                    lastGood=m-1;
                    cout<<"m="<<m<<endl;
                    }
                if(oltreBKG){
                    hist_BDT_S->SetBinContent(lastGood,hist_BDT_S->GetBinContent(lastGood)+ hist_BDT_S->GetBinContent(m));
                    hist_BDT_S->SetBinContent(m,0);
                    hist_BDT_B->SetBinContent(lastGood,hist_BDT_B->GetBinContent(lastGood)+hist_BDT_B->GetBinContent(m));
                    hist_BDT_B->SetBinContent(m,0);
                }                
        } 

                
                 oltreBKG = false;
                 lastGood=0;
                for(int m=1; m <= 40;++m){
                if( hist_BDT_B_train->GetBinContent(m)<0.5 && oltreBKG == false) {
                    oltreBKG = true;
                    lastGood=m-1;
                    cout<<"m="<<m<<endl;
                    }
                if(oltreBKG){
                    hist_BDT_S_train->SetBinContent(lastGood,hist_BDT_S_train->GetBinContent(lastGood)+ hist_BDT_S_train->GetBinContent(m));
                    hist_BDT_S_train->SetBinContent(m,0);
                    hist_BDT_B_train->SetBinContent(lastGood,hist_BDT_B_train->GetBinContent(lastGood)+hist_BDT_B_train->GetBinContent(m));
                    hist_BDT_B_train->SetBinContent(m,0);
                }                
        } */ 
                
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
                
                
                
                
                    
                hist_BDT_S->Scale(9.2/hist_BDT_S->Integral());
                hist_BDT_B->Scale(8610./hist_BDT_B->Integral());
                hist_BDT_S_train->Scale(9.2/hist_BDT_S_train->Integral());
                hist_BDT_B_train->Scale(8610./hist_BDT_B_train->Integral());
                hist_BDT_S_train_testBinning->Scale(9.2/hist_BDT_S_train_testBinning->Integral());
                hist_BDT_B_train_testBinning->Scale(8610./hist_BDT_B_train_testBinning->Integral());
                

                              


        
        
                
                DrawAndSave_bkg_and_signal(hist_BDT_S,       hist_BDT_B,        "Test");
                DrawAndSave_bkg_and_signal(hist_BDT_S_train, hist_BDT_B_train,  "Train");
                
            
                totalSensitivityErrorSquared_thisVariable[seed_idx] = 0.;
                totalSensitivity_thisVariable[seed_idx] = 0.;
                totalSensitivity_thisVariable_train[seed_idx] = 0.;
                
                
//                 computeSensitivityAndError(hist_BDT_S, hist_BDT_B, hist_sensitivity, seed_idx, totalSensitivity_thisVariable, totalSensitivityError_thisVariable);
                
                
                std::cout << "totalSensitivity_thisVariable  prima:    " << totalSensitivity_thisVariable[seed_idx] << std::endl;
                
                for (int n = 1; n <= hist_BDT_S->GetXaxis()->GetNbins(); n++) {
    //                 double SignalScan = 0.;
    //                 double BackgroundScan = 0.;
                
                    double SignalScan = hist_BDT_S->GetBinContent(n);
                    double BackgroundScan = hist_BDT_B->GetBinContent(n);
                    double Signal_RelativeError = hist_BDT_S->GetBinError(n)/SignalScan;
                    double BKG_RelativeError = hist_BDT_B->GetBinError(n)/BackgroundScan;

                    totalSensitivity_thisVariable[seed_idx] +=  ((BackgroundScan > 0.0001) && (SignalScan > 0.0001)) ? SignalScan*SignalScan/(BackgroundScan) : 0.;
                    hist_sensitivity->SetBinContent(n, ((BackgroundScan > 0.0001) && (SignalScan > 0.0001)) ? SignalScan/sqrt(BackgroundScan) : 0.); 
                
                    float SensitivityRelativeErrorSquared = (BackgroundScan > 0.00001 && SignalScan > 0.00001) ? (4*Signal_RelativeError*Signal_RelativeError + BKG_RelativeError*BKG_RelativeError) : 0.;
                    totalSensitivityErrorSquared_thisVariable[seed_idx] += SensitivityRelativeErrorSquared*hist_sensitivity->GetBinContent(n)*hist_sensitivity->GetBinContent(n);
                    
                    std::cout << "totalSensitivityErrorSquared   "  << hist_BDT_B->GetBinError(n) <<  "  \t  "  << BackgroundScan <<  "  \t  "  << BKG_RelativeError <<  "  \t  "  << SensitivityRelativeErrorSquared <<  "  \t  "  <<hist_sensitivity->GetBinContent(n)*hist_sensitivity->GetBinContent(n) <<  "  \t  "  << totalSensitivityErrorSquared_thisVariable[seed_idx] << std::endl;
                    
                    double SignalScan_train = hist_BDT_S_train->GetBinContent(n);
                    double BackgroundScan_train = hist_BDT_B_train->GetBinContent(n);
                    totalSensitivity_thisVariable_train[seed_idx] +=  ((BackgroundScan_train > 0.0001) && (SignalScan_train > 0.0001)) ? SignalScan_train*SignalScan_train/(BackgroundScan_train) : 0.;
                    hist_sensitivity_train->SetBinContent(n, ((BackgroundScan_train > 0.0001) && (SignalScan_train > 0.0001)) ? SignalScan_train/sqrt(BackgroundScan_train) : 0.);
                    
//                     float SensitivityRelativeErrorSquared_train = (BackgroundScan > 0.00001 && SignalScan > 0.00001) ? (4*Signal_RelativeError*Signal_RelativeError + BKG_RelativeError*BKG_RelativeError) : 0.;
//                     totalSensitivityErrorSquared_thisVariable_train[seed_idx] += SensitivityRelativeErrorSquared_train*hist_sensitivity_train->GetBinContent(n)*hist_sensitivity_train->GetBinContent(n);
                }
                
                totalSensitivity_thisVariable[seed_idx] = sqrt(totalSensitivity_thisVariable[seed_idx]);
                totalSensitivity_thisVariable_train[seed_idx] = sqrt(totalSensitivity_thisVariable_train[seed_idx]);

                std::cout << "totalSensitivity_thisVariable  dopo:     " << totalSensitivity_thisVariable[seed_idx] << std::endl;
                std::cout << "totalSensitivityError_thisVariable  prima:    " << totalSensitivityError_thisVariable[seed_idx] << std::endl;

                totalSensitivityError_thisVariable[seed_idx] = sqrt(totalSensitivityErrorSquared_thisVariable[seed_idx]);
                std::cout << "totalSensitivityError_thisVariable  prima:    " << totalSensitivityError_thisVariable[seed_idx] << std::endl;


                DrawSensitivity(hist_sensitivity, totalSensitivity_thisVariable[seed_idx], "Test");
                DrawSensitivity(hist_sensitivity_train, totalSensitivity_thisVariable_train[seed_idx], "Train");

                
                
                double X_S= hist_BDT_S_train_testBinning->Chi2Test(hist_BDT_S,"WWCHI2/NDF");
                double X_B= hist_BDT_B_train_testBinning->Chi2Test(hist_BDT_B,"WWCHI2/NDF");
                double prob_S=TMath::Prob(hist_BDT_S_train_testBinning->Chi2Test(hist_BDT_S,"WWCHI2"),hist_BDT_S_train_testBinning->GetNbinsX()-1);
                double prob_B=TMath::Prob(hist_BDT_B_train_testBinning->Chi2Test(hist_BDT_B,"WWCHI2"),hist_BDT_B_train_testBinning->GetNbinsX()-1);
                double K_s= hist_BDT_S_train_testBinning->KolmogorovTest(hist_BDT_S);
                double K_b= hist_BDT_B_train_testBinning->KolmogorovTest(hist_BDT_B);



        

        
        
                drawTestAndTrainSuperimposed(hist_BDT_B, hist_BDT_B_train_testBinning,"BDT_bkg"+variables_names[current_file]+"_"+seed_string+"Seed"+".png",K_b,prob_B);
                drawTestAndTrainSuperimposed(hist_BDT_S, hist_BDT_S_train_testBinning,"BDT_signal"+variables_names[current_file]+"_"+seed_string+"Seed"+".png",K_s,prob_S);

                
                totalChiS_thisVariable[seed_idx]=X_S;
                totalChiB_thisVariable[seed_idx]=X_B;
                totalProbS_thisVariable[seed_idx]=prob_S;
                totalProbB_thisVariable[seed_idx]=prob_B;
                totalKS_S_thisVariable[seed_idx]=K_s;
                totalKS_B_thisVariable[seed_idx]=K_b;

                
            }
        

        
        
        totalSensitivity[current_file]              = findMean(totalSensitivity_thisVariable, numberOfSeedToRead);
        totalSensitivity_train[current_file]        = findMean(totalSensitivity_thisVariable_train, numberOfSeedToRead);
        totalChiS[current_file]                     = findMean(totalChiS_thisVariable, numberOfSeedToRead);
        totalChiB[current_file]                     = findMean(totalChiB_thisVariable, numberOfSeedToRead);
        totalProbS[current_file]                    = findMean(totalProbS_thisVariable, numberOfSeedToRead);
        totalProbB[current_file]                    = findMean(totalProbB_thisVariable, numberOfSeedToRead);
        totalKS_S[current_file]                     = findMean(totalKS_S_thisVariable, numberOfSeedToRead);
        totalKS_B[current_file]                     = findMean(totalKS_B_thisVariable, numberOfSeedToRead);
        
        
        
        
        drawSinglePlot(totalSensitivity_thisVariable, numberOfSeedToRead, "Sensitivity test", "BDToutput/sensitivityTest/seedDistribution/sensitivity_" + variables_names[current_file] + "_Test", totalSensitivity[current_file]*0.5, totalSensitivity[current_file]*1.4, 40);
        drawSinglePlot(totalSensitivity_thisVariable_train, numberOfSeedToRead, "Sensitivity train", "BDToutput/sensitivityTrain/seedDistribution/sensitivity_" + variables_names[current_file] + "_Train", totalSensitivity_train[current_file]*0.5, totalSensitivity_train[current_file]*1.4, 40);
        
        drawSinglePlot(totalChiS_thisVariable, numberOfSeedToRead, "Sensitivity test", "Chi2_and_KS_prob/Chi2/Chi2_" + variables_names[current_file] + "_S", 0., 1., 40);
        drawSinglePlot(totalChiB_thisVariable, numberOfSeedToRead, "Sensitivity train", "Chi2_and_KS_prob/Chi2/Chi2_" + variables_names[current_file] + "_B", 0., 1., 40);
        drawSinglePlot(totalProbS_thisVariable, numberOfSeedToRead, "Chi2 probability signal", "Chi2_and_KS_prob/Chi2Probability/Prob_" + variables_names[current_file] + "_S", 0., 1., 40);
        drawSinglePlot(totalProbB_thisVariable, numberOfSeedToRead, "Chi2 probability bkg", "Chi2_and_KS_prob/Chi2Probability/Prob_" + variables_names[current_file] + "_B", 0., 1., 40);
        drawSinglePlot(totalKS_S_thisVariable, numberOfSeedToRead, "KS signal", "Chi2_and_KS_prob/KS/KS_" + variables_names[current_file] + "_S", 0., 1., 40);
        drawSinglePlot(totalKS_B_thisVariable, numberOfSeedToRead, "KS bkg", "Chi2_and_KS_prob/KS/KS_" + variables_names[current_file] + "_B", 0., 1., 40);
        


            
        
        sensitivity_ey[current_file] = findRMS(totalSensitivity_thisVariable, numberOfSeedToRead);
        sensitivity_train_ey[current_file] = findRMS(totalSensitivity_thisVariable, numberOfSeedToRead);
            
            
            
        Double_t sensitivityDifference_thisVariable_test_and_train[numberOfSeedToRead];
        for(int n=0; n < numberOfSeedToRead; n++) {
            sensitivityDifference_thisVariable_test_and_train[n] = totalSensitivity_thisVariable_train[n] - totalSensitivity_thisVariable[n];
        }
        drawSinglePlot(sensitivityDifference_thisVariable_test_and_train, numberOfSeedToRead, "sensitivity difference train - test", "BDToutput/sensitivityDifference/sensitivity_" + variables_names[current_file] + "_difference", -0.2, 0.2, 20);
            
            
            
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++)  std::cout << totalChiS_thisVariable[seed_idx] << " \t ";
                std::cout <<  std::endl;
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++)  std::cout << totalChiB_thisVariable[seed_idx] << " \t ";
                std::cout <<  std::endl;
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++)  std::cout << totalProbS_thisVariable[seed_idx] << " \t ";
                std::cout <<  std::endl;
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++)  std::cout << totalProbB_thisVariable[seed_idx] << " \t ";
                std::cout <<  std::endl;
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++)  std::cout << totalKS_S_thisVariable[seed_idx] << " \t ";
                std::cout <<  std::endl;
            for (int seed_idx=0;seed_idx<numberOfSeedToRead;seed_idx++)  std::cout << totalKS_B_thisVariable[seed_idx] << " \t ";
                std::cout <<  std::endl;
             
            
                
        }

        
        
        
        
        gPad->SetGridx();
        gPad->SetGridy();



        gStyle->SetOptStat(0000);
        TH1F *frame2 = new TH1F("frame2","",n_variables,0.,n_variables);
        TH1F *hChi_S = new TH1F("hChi_S","",n_variables,0.,n_variables);
        TH1F *hChi_B = new TH1F("hChi_B","",n_variables,0.,n_variables);
        TH1F *hProb_S= new TH1F("hProb_S","",n_variables,0.,n_variables);
        TH1F *hProb_B= new TH1F("hProb_B","",n_variables,0.,n_variables);
        TH1F *hKS_S= new TH1F("hKS_S","",n_variables,0.,n_variables);
        TH1F *hKS_B= new TH1F("hKS_B","",n_variables,0.,n_variables);
        
        
        frame2->SetMinimum(totalSensitivity[n_variables-1]*0.);
        frame2->SetMaximum(totalSensitivity[n_variables-1]*2.);



        
        
        for (int i=0;i<n_variables;i++){
            frame2->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            hChi_S->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            hChi_B->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            hProb_S->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            hProb_B->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            hKS_S->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            hKS_B->GetXaxis()->SetBinLabel(i+1,variables_names[i].c_str());
            
            
            
            frame2->GetXaxis()->SetLabelSize(0.02);
            hChi_S->GetXaxis()->SetLabelSize(0.02);
            hChi_B->GetXaxis()->SetLabelSize(0.02);
            hProb_S->GetXaxis()->SetLabelSize(0.02);
            hProb_B->GetXaxis()->SetLabelSize(0.02);
            hKS_S->GetXaxis()->SetLabelSize(0.02);
            hKS_B->GetXaxis()->SetLabelSize(0.02);
                
        }
        
        frame2->GetXaxis()->LabelsOption("v");
        
        setHistoStyle(frame2,    variable_used_in_BDT.c_str(),                    "Sensitivity");
        setHistoStyle(hChi_S,    "Signal Chi^2",        "Chi2");
        setHistoStyle(hChi_B,    "Background Chi^2",    "Chi2");
        setHistoStyle(hProb_S,   "Signal Prob",         "Probability (Chi2)");
        setHistoStyle(hProb_B,   "Background Prob",     "Probability (Chi2)");
        setHistoStyle(hKS_S,     "Signal KS",         "KS");
        setHistoStyle(hKS_B,     "Background KS",     "KS");
        

        
        

        TGraphErrors *gr_train    = new TGraphErrors(n_variables,frame2_axisx_train,totalSensitivity_train, sensitivity_ex, sensitivity_train_ey);
        TGraphErrors *gr          = new TGraphErrors(n_variables,frame2_axisx,totalSensitivity, sensitivity_ex, sensitivity_ey);
        TGraph *gr_ProbS    = new TGraph(n_variables,frame2_axisx,totalProbS);
        TGraph *gr_ProbB    = new TGraph(n_variables,frame2_axisx,totalProbB);
        TGraph *gr_ChiS     = new TGraph(n_variables,frame2_axisx,totalChiS);
        TGraph *gr_ChiB     = new TGraph(n_variables,frame2_axisx,totalChiB);
        TGraph *gr_KS_S     = new TGraph(n_variables,frame2_axisx,totalKS_S);
        TGraph *gr_KS_B     = new TGraph(n_variables,frame2_axisx,totalKS_B);
        

        gr_train->SetMarkerColor(4);
        gr_train->SetMarkerStyle(21);
        gr->SetMarkerColor(4);
        gr->SetMarkerStyle(21);
        
        drawGraph(canv, frame2,  gr, "totalSensitivity"+trainingOption,             true, n_variables, totalSensitivity[n_variables-1]);        
        drawGraph(canv, frame2,  gr_train, "totalSensitivity_train"+trainingOption, true, n_variables, totalSensitivity_train[n_variables-1]);        
        drawGraph(canv, hProb_S, gr_ProbS, "Chi2_and_KS_prob/probability_signal",            false, 0, 0);
        drawGraph(canv, hProb_B, gr_ProbB, "Chi2_and_KS_prob/probability_background",        false, 0, 0);
        drawGraph(canv, hChi_S,  gr_ChiS, "Chi2_and_KS_prob/chi2_signal",                    false, 0, 0);
        drawGraph(canv, hChi_B,  gr_ChiB, "Chi2_and_KS_prob/chi2_background",                false, 0, 0);
        drawGraph(canv, hKS_S,   gr_KS_S, "Chi2_and_KS_prob/KS_signal",                      false, 0, 0);
        drawGraph(canv, hKS_B,   gr_KS_B, "Chi2_and_KS_prob/KS_background",                  false, 0, 0);



//         drawAndSaveSensitivities(totalSensitivity, totalSensitivity_train, n_variables, "test:  ", "train: ", "sensitivity distributions",  "sensitivityDistribution", 0.45, 0.7, 40);
//         drawAndSaveSensitivities(totalKS_S, totalKS_B, n_variables, "KS signal:  ", "KS bkg:    ", "KS distributions", "KSDistribution", 0.9, 1., 40);
//         drawAndSaveSensitivities(totalProbB, totalProbS, n_variables, "Chi2 signal:  ", "Chi2 bkg:    ", "Chi2 distributions", "ChiDistribution", 0.9, 1., 40);
//        
//         
//         Double_t sensitivityDifference_test_and_train[n_variables];
//         for(int n=0; n < n_variables; n++) {
//             sensitivityDifference_test_and_train[n] = totalSensitivity_train[n] - totalSensitivity[n];
//         }
//         
//         drawSinglePlot(sensitivityDifference_test_and_train, n_variables, "sensitivity difference train - test", "sensitivityDistribution_difference", -0.2, 0.2, 40);

    
}





