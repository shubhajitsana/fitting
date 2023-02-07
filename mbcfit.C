//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooBifurGauss.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "TH1.h"
#include "TH2.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "TChain.h"
#include<cmath>
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooNLLVar.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace RooFit ;
using namespace std;
//int main(){

void mbcfit(){
    /*******************Fit Variables***********************************/
    RooRealVar mbc("mbc","M_{bc} (GeV)",5.27,5.29);
    /**defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT**/
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(mbc));
    /*******************Input root file**********************************/
    TChain* chain=new TChain();
    chain->Add("/home/sana/ssana/fitting/charged.root/tree");

    Double_t  de3, md03, mbc3, r23, kid3,pid3,sig,cont_prob;
    Int_t run;
    Int_t nevt3=(int)chain->GetEntries();

    // chain->SetBranchAddress("isSignal",&sig);
    chain->SetBranchAddress("deltaE",&de3);
    chain->SetBranchAddress("Mbc",&mbc3);
    chain->SetBranchAddress("D0_bar_InvM",&md03);
    chain->SetBranchAddress("ContProb",&cont_prob);
    chain->SetBranchAddress("Kp_PID_bin_kaon",&kid3);
    // chain->SetBranchAddress("R2",&r23);
    // chain->SetBranchAddress("pi_PID_bin_pion",&pid3);
    // chain->SetBranchAddress("__run__",&run);
    

    //Loading data 
    Double_t counter =0;
    for(int l=0;l<nevt3;l++) {
      chain->GetEntry(l);
      mbc.setVal(mbc3);
      if(md03>1.84 && md03<1.89 && mbc3>5.27 && mbc3 < 5.29 && de3 < 0.1 && de3 > -0.1 && kid3 > 0.6 && cont_prob < 0.86){//&& r23 < 0.3  )// && (run <= 1702 || run >= 1835))
          data->add(RooArgSet(mbc));
          counter++;
      }
    }

    /*****************************Delta E Fit***********************/
    // --- Build Signal PDF ---
    RooRealVar mean1("mean1","mean of Gaussian-1",5.279145,5.27,5.29);
    RooRealVar mean2("mean2","mean of Gaussian-2",5.279260,5.27,5.29);
    RooRealVar sigma1("sigma1","sigma of Gaussian-1",0., 0., 1);	
    RooRealVar sigma2("sigma2","sigma of Gaussian-2",0.01191,0.00000001,1);

    RooGaussian sig1("sig1","Gaussian-1",mbc,mean1,sigma1);  
    RooGaussian sig2("sig2","Gaussian-2",mbc,mean2,sigma2);

    RooRealVar fsig_1("frac_gaussians", "signal fraction", 0.5,0.,1.);
    RooAddPdf twoGaussians("twoGaussians", "sum of two Gaussians ",RooArgList(sig1, sig2), RooArgList(fsig_1));

    // --- Build Argus background PDF ---
  RooRealVar argpar("argpar","argus shape parameter",-34.70,-100.,-1.) ;
  RooArgusBG bkg("argus","Argus PDF",mbc,RooConst(5.29),argpar) ; //mbc background

    //Initialization of parameter before adding two pdf
    // cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
    Double_t event_count = counter; 
    Double_t signal_count = counter*0.4;
    Double_t back_count = counter*0.6;
    RooRealVar n_sig("n_sig", "n_sig", signal_count, 0., event_count);//52000
    RooRealVar n_bkg("n_bkg", "n_bkg", back_count, 0., event_count);//95000
    RooAddPdf sum("sum","sum",RooArgList(twoGaussians,bkg),RooArgList(n_sig, n_bkg));//adding two pdf
    sum.fitTo(*data);
    /****************************FIT COMPLETE*************************************/

    /*********************Start Plotting and showing outpouts*****************/
    //Plotting fitted result
    RooPlot* deframe = mbc.frame(Title("Fitting M_{bc} of B^{#pm}"), Bins(550)) ;                          
    data->plotOn(deframe, Binning(550), DataError(RooAbsData::SumW2)) ;
    sum.plotOn(deframe, LineColor(kBlue)	, LineStyle(kSolid)) ;
    // sum.paramOn(deframe,data,"", 2, Format("NEU"),AutoPrecision(1)), 0.7, 0.9, 0.9); //"NELU",  Prints the fitted parameter on the canvas
    sum.paramOn(deframe,data,"", 2, "NU", 0.1, 0.3, 0.9); //"NELU",  Prints the fitted parameter on the canvas
    sum.plotOn(deframe,Components(sig1),LineColor(kGreen),LineStyle(kDashed)) ;
    sum.plotOn(deframe,Components(sig2),LineColor(kBlack),LineStyle(kDashed)) ;
    sum.plotOn(deframe,Components(twoGaussians),LineColor(kRed),LineStyle(kDashed));
    sum.plotOn(deframe,Components(bkg),LineColor(kMagenta),LineStyle(kDashed)) ;

    //Extract info. from fitting frame and showing
    cout<<"chisq of the fit is : "<<deframe->chiSquare()<<endl;//chi-square of the fit
    cout<<"chi-square/ndof : "<<deframe->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)
    RooHist* hpull = deframe->pullHist() ;
    RooPlot* frame3 = mbc.frame(Title("Pull Distribution")) ;
    hpull->SetFillColor(1);
    frame3->addPlotable(hpull,"X0B"); // "X0" is for errorbar; and "B" is for histograms



    TCanvas* c1 = new TCanvas("c1", "c1", 2550, 1500) ;
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();  
    deframe->Draw() ;
    // Adding legend
    // TLegend *legend1 = new TLegend(0.1,0.7,0.3,0.9);
    // TLegendEntry *entry = legend1->AddEntry("sig1","1st Gaussian pdf","l");
    // entry->SetLineColor(kGreen);
    // entry->SetLineStyle(kDashed);
    // entry = legend1->AddEntry("sig2","2nd Gaussian pdf","l");
    // entry->SetLineColor(kBlack);
    // entry->SetLineStyle(kDashed);
    // entry = legend1->AddEntry("twoGaussians","Combined signal pdf","l");
    // entry->SetLineColor(kRed);
    // entry->SetLineStyle(kDashed);
    // entry = legend1->AddEntry("bkg","bkg-Chebyshev","l");
    // entry->SetLineColor(kMagenta);
    // entry->SetLineStyle(kDashed);
    // entry = legend1->AddEntry("sum","Fitted pdf","l");
    // entry->SetLineColor(kBlue);
    // entry->SetLineStyle(kSolid);
    // legend1->Draw();

    c1->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->Draw();
    pad2->cd(); 
    frame3->Draw() ;
    frame3->SetLineStyle(9);
    frame3->GetYaxis()->SetNdivisions(505);
    frame3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}"); 
    frame3->GetXaxis()->SetTitle("M_{bc} (GeV)"); 
    frame3->GetXaxis()->SetTitleSize(0.13);
    frame3->GetYaxis()->SetTitleSize(0.15);
    frame3->GetXaxis()->SetLabelSize(0.120);
    frame3->GetYaxis()->SetLabelSize(0.120); 
    frame3->GetXaxis()->SetTitleOffset(0.90);      
    frame3->GetYaxis()->SetTitleOffset(0.22);       
    frame3->GetYaxis()->SetRangeUser(-10.0, 10.0);       
    frame3->GetYaxis()->SetLimits(-10.0, 10.0);       
    frame3->GetXaxis()->SetNdivisions(505);
    frame3->GetYaxis()->CenterTitle(true);
    frame3->GetXaxis()->CenterTitle(true);
    frame3->Draw("AXISSAME");

    c1->Print("mbc_plot/mbc_fit_c_2.png");
    cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
    cout<<"chisq of the fit is :"<<deframe->chiSquare()<<endl;//chi-square of the fit
    cout<<"chi-square/ndof :"<<deframe->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)
}