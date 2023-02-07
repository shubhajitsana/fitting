#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Fetching problem in paramOn


import ROOT

######################Fit Variables#################
deltae = ROOT.RooRealVar("deltae","#DeltaE (GeV)", -0.1, 0.1)
####defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT
data = ROOT.RooDataSet("data", "data", ROOT.RooArgSet(deltae))

# Loading input root file and creating new root file
inFile = ROOT.TFile.Open("/home/sana/ssana/fitting/charged.root")
inTree = inFile.Get('tree')
NEvent = inTree.GetEntries()

# Loading data 
counter =0
for iEvent in range(inTree.GetEntries()):
    inTree.GetEntry(iEvent)
    # signal = getattr(inTree, 'isSignal')
    o_de = getattr(inTree, 'deltaE')
    o_mbc = getattr(inTree, 'Mbc')
    o_md0 = getattr(inTree, 'D0_bar_InvM')
    # o_r2 = getattr(inTree, 'R2')
    # o_kid = getattr(inTree, 'Kp_PID_bin_kaon')
    # o_ = getattr(inTree, '')
    if((o_md0>1.84) and (o_md0<1.89) and (o_mbc>5.27) and (o_mbc < 5.29) and (o_de < 0.1) and (o_de > -0.1)): # and (o_kid > 0.6)(o_r2 < 0.3) and 
        data.add({o_de})
        counter += 1


# //Drawing casual Histogram(BINNED) before fitting
# binDataSet = ROOT.RooDataHist("binDataSet", "binDataSet", mbc, data)
# c1 = ROOT.TCanvas("c1", "", 1500, 1500)  
# xframe1 = mbc.frame(Title="prefit histooo", Bins=200)
# binDataSet.plotOn(xframe1, Binning=200, DataError="SumW2")
# xframe1.Draw()
# c1.SaveAs("py_prefit_histo.png")


#***********************************Mbc fit****************************
#--- Build Signal PDF ---
mean1 = ROOT.RooRealVar("mean1","mean of Gaussian-1",-0.001,-0.02,0.02)
mean2 = ROOT.RooRealVar("mean2","mean of Gaussian-2",0.002,-0.02,0.02)
sigma1 = ROOT.RooRealVar("sigma1","sigma of Gaussian-1",0., 0., 0.05)
sigma2 = ROOT.RooRealVar("sigma2","sigma of Gaussian-2",0.01191,0.00000001,0.05)

sig1 = ROOT.RooGaussian("sig1","Gaussian-1",deltae,mean1,sigma1)
sig2 = ROOT.RooGaussian("sig2","Gaussian-2",deltae,mean2,sigma2)

fsig_1 = ROOT.RooRealVar("fsig_1", "signal fraction", 0.4,0.,1.)
twoGaussians = ROOT.RooAddPdf("twoGaussians", "sum of two Gaussians ", [sig1, sig2], [fsig_1])
# // --- Build Argus background PDF ---
b1 = ROOT.RooRealVar("Chbyshv-prm", "Chbyshv-prm", -0.062, -10., 10.)
bkg = ROOT.RooChebychev("bkg","Background",deltae, [b1]) 
    
# //Initialization of parameter before adding two pdf
event_count = counter 
signal_count = counter*0.4
back_count = counter*0.6
n_sig = ROOT.RooRealVar("n_sig", "n_sig", signal_count, 0., event_count)#52000
n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", back_count, 0., event_count)#95000
sum = ROOT.RooAddPdf("sum","sum", [twoGaussians,bkg], [n_sig, n_bkg])#adding two pdf

sum.fitTo(data)#fitting

#   /*********************Start Plotting and showing outpouts*****************/
#   //Plotting fitted result
deframe = deltae.frame(Title="Fitting #DeltaE", Bins=200) 
data.plotOn(deframe, Binning=200, DataError="SumW2") 
sum.plotOn(deframe, LineColor="kBlue")#, LineStyle(kSolid)) 
#Fetching problem in paramOn
sum.paramOn(deframe,data, FillColor="kRed",  Layout=(0.65, 0.9, 0.9),Format="NU", ShowConstants=True)#Prints the fitted parameter on the canvas
sum.plotOn(deframe,Components="sig1",LineColor="kGreen",LineStyle="--") 
sum.plotOn(deframe,Components="sig2",LineColor="kBlack",LineStyle="--") 
sum.plotOn(deframe,Components="twoGaussians",LineColor="kRed",LineStyle="--")
sum.plotOn(deframe,Components="bkg",LineColor="kMagenta",LineStyle="--") 


# //Extract info. from fitting frame and showing
# chisq = frame.chiSquare()#extract chi2 value
print(f"chi2 of mbc fitting = {frame.chiSquare()}") #Printing chi2 value
print(f"chi-square/ndof : {frame.chiSquare(7)}") #Printing chi2 value
hpull = deframe.pullHist() 
frame3 = deltae.frame(Title="Pull Distribution") 
hpull.SetFillColor(1)
frame3.addPlotable(hpull,"X0B") #"X0" is for errorbar; and "B" is for histograms




# //Start drawing fitting results
c1 = ROOT.TCanvas("c1", "c1", 2550, 1500) 
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.Draw()             # Draw the upper pad: pad1
pad1.cd()  
deframe.Draw()
# legend1 = TLegend(0.1,0.7,0.3,0.9)
# TLegendEntry *entry = legend1.AddEntry("sig1","1st Gaussian pdf","l")
# entry.SetLineColor(kGreen)
# entry.SetLineStyle(kDashed)
# entry = legend1.AddEntry("sig2","2nd Gaussian pdf","l")
# entry.SetLineColor(kBlack)
# entry.SetLineStyle(kDashed)
# entry = legend1.AddEntry("twoGaussians","Combined signal pdf","l")
# entry.SetLineColor(kRed)
# entry.SetLineStyle(kDashed)
# entry = legend1.AddEntry("bkg","bkg-Chebyshev","l")
# entry.SetLineColor(kMagenta)
# entry.SetLineStyle(kDashed)
# entry = legend1.AddEntry("sum","Fitted pdf","l")
# entry.SetLineColor(kBlue)
# entry.SetLineStyle(kSolid)
# legend1.Draw()

c1.cd()          # Go back to the main canvas before defining pad2
pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
pad2.Draw()
pad2.cd() 
frame3.Draw() 
# frame3.SetLineStyle(9)
# frame3.GetYaxis().SetNdivisions(505)
# frame3.GetYaxis().SetTitle("#sqrt{#chi^{2}}");
# frame3.GetXaxis().SetTitle("M_{bc} (GeV)") 
# frame3.GetXaxis().SetTitleSize(0.13)
# frame3.GetYaxis().SetTitleSize(0.15)
# frame3.GetXaxis().SetLabelSize(0.120)
# frame3.GetYaxis().SetLabelSize(0.120) 
# frame3.GetXaxis().SetTitleOffset(0.90)      
# frame3.GetYaxis().SetTitleOffset(0.22)       
# frame3.GetYaxis().SetRangeUser(-10.0, 10.0)       
# frame3.GetYaxis().SetLimits(-10.0, 10.0)       
# frame3.GetXaxis().SetNdivisions(505)
# frame3.GetYaxis().CenterTitle(true)
# frame3.GetXaxis().CenterTitle(true)
# frame3.Draw("AXISSAME")

c2.SaveAs("py_mbcfit.pdf")