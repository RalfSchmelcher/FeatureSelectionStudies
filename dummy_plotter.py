import ROOT, os, subprocess, sys, optparse, datetime
from array import array
import math
import numpy as np
date = datetime.date.today().isoformat()

def histStyle(hist,xtitle,color):
    markerstyle = 20 
    #hist.SetMarkerStyle(markerstyle)
    #hist.SetMarkerColor(cmode)
    hist.SetLineColor  (color)
    hist.SetLineWidth  (2)
    hist.SetLineStyle  (1)
    hist.SetMarkerSize(0.8)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle('a.u.')
    #hist.Scale(1.0/hist.Integral())
    return hist


def gethist(tree,selcuts,vname,xtitle,nXbins,xmin,xmax, color = 1, Xbins=[]):
    if Xbins != []:
        hist1 = ROOT.TH1F("hist1" ,'',len(Xbins)-1,Xbins)
    else:
        hist1 = ROOT.TH1F("hist1" ,'',nXbins,xmin,xmax)
    if vname.startswith('d'):
        tree.Draw('abs({here})>>hist1'.format(here=vname), selcuts, 'goff')
    else :
        tree.Draw('{here}>>hist1'.format(here=vname), selcuts, 'goff')
    #tree.Draw('{here}>>hist1'.format(here=abs(vname) if vname.startswith('d') else vname),selcuts, 'goff')
    #tree.Draw('{here}>>hist1'.format(here=vname), '', 'goff')
    hist1 = histStyle(hist1,xtitle,color)
    hist1.SetDirectory(0)
    return hist1

def get2dhist(tree, selcuts, vnamex, vnamey, xtitle, ytitle, nXbins, nYbins, xmin, ymin, xmax, ymax, Xbins, Ybins):
    if Xbins!=[] and Ybins!=[]: # Xbins and Ybins arrays instead
        hist2d = ROOT.TH2F("hist2d",'',len(Xbins)-1,Xbins,len(Ybins)-1,Ybins)
        # print("1st option")
    elif Xbins!=[]: # Xbins array instead of xmin & xmax
        hist2d = ROOT.TH2F("hist2d",'',len(Xbins)-1,Xbins,nYbins,ymin,ymax)
        # print("2nd option")
    elif Ybins!=[]: # Ybins array instead of ymin & ymax
        hist2d = ROOT.TH2F("hist2d",'',nXbins,xmin,xmax,len(Ybins)-1,Ybins)
        # print("3rd option")
    else: # usual
        hist2d = ROOT.TH2F("hist2d" ,'',nXbins,xmin,xmax,nYbins,ymin,ymax)
        # print("4th option")
    tree.Draw("{herey}:{herex}>>hist2d".format(herex=vnamex,herey=vnamey), selcuts, "goff")
    hist2d.GetXaxis().SetTitle(xtitle)
    hist2d.GetYaxis().SetTitle(ytitle)
    print(hist2d.Integral())
    hist2d.SetDirectory(0)
    return hist2d

sel_str_Gen_matching = "((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1)) || ((Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1) && (Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0))"
sel_str_Gen_paper = "GenJetAK8_pt[0]>200 && GenJetAK8_pt[1]>200 && abs(GenJetAK8_eta[0])<2.5 && abs(GenJetAK8_eta[1])<2.5 && Gen_M_jj>1250 && GenJetAK8_mass[0]>30 && GenJetAK8_mass[1]>30"
sel_str_Gen_matching_paper = "(" + sel_str_Gen_paper + ") && (" + sel_str_Gen_matching + ")"
sel_str_Reco_matching = "((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1)) || ((Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1) && (Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0))"
sel_str_Reco_paper = "FatJet_pt[0]>200 && FatJet_pt[1]>200 && abs(FatJet_eta[0])<2.5 && abs(FatJet_eta[1])<2.5 && Reco_M_jj>1250 && FatJet_mass[0]>30 && FatJet_mass[1]>30"
sel_str_Reco_matching_paper = "(" + sel_str_Reco_paper + ") && (" + sel_str_Reco_matching + ")"

binsForX = array('f',[1000, 1500, 2000, 3000, 4000, 5000, 6000, 6500]) # +1 for last bin
binsForY = array('f',[ 30,  60,  90, 120, 150, 200, 250, 300, 400, 500, 550]) # +1 for last bin
binsForYP = array('f',[ 30,  60,  90, 120, 150, 200, 250, 300, 400, 500, 550]) # +1 for last bin
# min_valueYX = min(binsForY[:-1]) / max(binsForX[:-1])
# max_valueYX = max(binsForY[:-1]) / min(binsForX[:-1])
# num_binsYX = len(binsForX[:-1]) * len(binsForY[:-1])
# step_sizeYX = (max_valueYX - min_valueYX) / num_binsYX
# binning_valuesYX = np.arange(min_valueYX, max_valueYX + 2*step_sizeYX, step_sizeYX)
# print(binning_valuesYX)

logBinning_YX = np.logspace(np.log10(1/200),np.log10(0.55),24)


def overlayshapes(outdir,l_ex0):

    pdraw={
        #### vname : [vname, nXbins, xmin, xmax, sel, _info]
        # "nV" : ["nV",100,0,4,"",""],
        # "V_pt" : ["V_pt",100,0,4000,"",""],
        # "V_eta" : ["V_eta",100,-3.5,3.5,"",""],
        # "V_phi" : ["V_phi",100,-3.5,3.5,"",""],
        # "V_mass" : ["V_mass",100,0,200,"",""],
        # "V_pdgId" : ["V_pdgId",100,22,26,"",""],
        "M_X" : ["M_X",50,0,6500,"",""],
        "M_Y" : ["M_Y",50,0,500,"",""],
        "M_YP" : ["M_YP",50,0,500,"",""],
        "delta_R_qq_Y23" : ["delta_R_qq_Y23",100,-0.1,2,"",""],
        "delta_R_qq_YP25" : ["delta_R_qq_YP25",100,-0.1,2,"",""],
        "Gen_M_jj" : ["Gen_M_jj",326,0,7000,"","no selection"],
        # "Gen_M_jj_match" : ["Gen_M_jj",326,0,7000,sel_str_Gen_matching,'matching'],
        # "Gen_M_jj_paper" : ["Gen_M_jj",326,0,7000,sel_str_Gen_paper, "VV/VH paper"],
        "Gen_M_jj_match_paper" : ["Gen_M_jj",326,0,7000,sel_str_Gen_matching_paper,"matching & VV/VH paper"],
        "Reco_M_jj" : ["Reco_M_jj",326,0,7000,"","no selection"],
        # "Reco_M_jj_match" : ["Reco_M_jj",326,0,7000,sel_str_Reco_matching,'matching'],
        # "Reco_M_jj_paper" : ["Reco_M_jj",326,0,7000,sel_str_Reco_paper, "VV/VH paper"],
        "Reco_M_jj_match_paper" : ["Reco_M_jj",326,0,7000,sel_str_Reco_matching_paper,"matching & VV/VH paper"],
        "Gen_M_j0Y" : ["GenJetAK8_mass[0]",100,0,1000,"((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1))",''],
        "Gen_M_j1Y" : ["GenJetAK8_mass[1]",100,0,1000,"((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1))",''],
        "Gen_M_j0YP" : ["GenJetAK8_mass[0]",100,0,1000,"((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1))",''],
        "Gen_M_j1YP" : ["GenJetAK8_mass[1]",100,0,1000,"((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1))",''],
        "Reco_M_j0Y" : ["FatJet_mass[0]",100,0,1000,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1))",''],
        "Reco_M_j1Y" : ["FatJet_mass[1]",100,0,1000,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1))",''],
        "Reco_M_j0YP" : ["FatJet_mass[0]",100,0,1000,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1))",''],
        "Reco_M_j1YP" : ["FatJet_mass[1]",100,0,1000,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1))",''],
        # "Reco_M_jY" : ["FatJet_mass",326,0,7000,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1)) || ((Reco_matched_Y[1]==1 && Reco_matched_Y[0]==0) && (Reco_matched_YP[1]==0 && Reco_matched_YP[0]==1))",''],
        # "Reco_M_jYP" : ["FatJet_mass",326,0,7000,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1)) || ((Reco_matched_YP[1]==1 && Reco_matched_YP[0]==0) && (Reco_matched_Y[1]==0 && Reco_matched_Y[0]==1))",''],
        # "nFatJet" : ["nFatJet",50,0,5,"",""],
        # "FatJet_pt" : ["FatJet_pt",100,0,4000,"",""],
        # "FatJet_eta" : ["FatJet_eta",100,-3.5,3.5,"",""],
        # "FatJet_phi" : ["FatJet_phi",100,-3.5,3.5,"",""],
        # "FatJet_mass" : ["FatJet_mass",100,0,700,"",""],
        # "FatJet_msoftdrop" : ["FatJet_msoftdrop",100,0,700,"",""],
        # "FatJet_n2b1" : ["FatJet_n2b1",100,-2,5,"",""],
        # "FatJet_n3b1" : ["FatJet_n3b1",100,-2,5,"",""],
        # "FatJet_tau1" : ["FatJet_tau1",100,0,1,"",""],
        # "FatJet_tau2" : ["FatJet_tau2",100,0,1,"",""],
        # "FatJet_tau3" : ["FatJet_tau3",100,0,1,"",""],
        # "FatJet_tau4" : ["FatJet_tau4",100,0,1,"",""],
        # "nSubJet" : ["nSubJet",50,0,5,"",""],
        # "SubJet_mass" : ["SubJet_mass",100,0,700,"",""],
        # "SubJet_n2b1" : ["SubJet_n2b1",100,-2,5,"",""],
        # "SubJet_n3b1" : ["SubJet_n3b1",100,-2,5,"",""],
        # "SubJet_tau1" : ["SubJet_tau1",100,0,1,"",""],
        # "SubJet_tau2" : ["SubJet_tau2",100,0,1,"",""],
        # "SubJet_tau3" : ["SubJet_tau3",100,0,1,"",""],
        # "SubJet_tau4" : ["SubJet_tau4",100,0,1,"",""],
        "GenJet0_deltaR_Y[0]"  : ["GenJet0_deltaR_Y[0]",100,-0.1,2,"",""],
        "GenJet1_deltaR_Y[0]"  : ["GenJet1_deltaR_Y[0]",100,-0.1,2,"",""],
        "GenJet0_deltaR_YP[0]" : ["GenJet0_deltaR_YP[0]",100,-0.1,2,"",""],
        "GenJet1_deltaR_YP[0]" : ["GenJet1_deltaR_YP[0]",100,-0.1,2,"",""],
        "FatJet0_deltaR_Y[0]"  : ["FatJet0_deltaR_Y[0]",100,-0.1,2,"",""],
        "FatJet1_deltaR_Y[0]"  : ["FatJet1_deltaR_Y[0]",100,-0.1,2,"",""],
        "FatJet0_deltaR_YP[0]" : ["FatJet0_deltaR_YP[0]",100,-0.1,2,"",""],
        "FatJet1_deltaR_YP[0]" : ["FatJet1_deltaR_YP[0]",100,-0.1,2,"",""],
        "GenJet0_deltaR_Y[1]"  : ["GenJet0_deltaR_Y[1]",100,-0.1,2,"",""],
        "GenJet1_deltaR_Y[1]"  : ["GenJet1_deltaR_Y[1]",100,-0.1,2,"",""],
        "GenJet0_deltaR_YP[1]" : ["GenJet0_deltaR_YP[1]",100,-0.1,2,"",""],
        "GenJet1_deltaR_YP[1]" : ["GenJet1_deltaR_YP[1]",100,-0.1,2,"",""],
        "FatJet0_deltaR_Y[1]"  : ["FatJet0_deltaR_Y[1]",100,-0.1,2,"",""],
        "FatJet1_deltaR_Y[1]"  : ["FatJet1_deltaR_Y[1]",100,-0.1,2,"",""],
        "FatJet0_deltaR_YP[1]" : ["FatJet0_deltaR_YP[1]",100,-0.1,2,"",""],
        "FatJet1_deltaR_YP[1]" : ["FatJet1_deltaR_YP[1]",100,-0.1,2,"",""],

        # tree->Draw("GenJetAK8_mass","(Gen_matched_YP[Iteration$]==1) && (Gen_matched_Y[abs(Iteration$-1)]==1)")


        }  
    sel=''
    for ikey,ival in pdraw.items():
        print(ikey)
        h_ex0 = gethist(l_ex0,ival[4],ival[0],ikey,ival[1],ival[2],ival[3],ROOT.kGreen+2)
        ROOT.gROOT.SetBatch()
        ROOT.gStyle.SetOptStat(0)
        canv = ROOT.TCanvas('{here}'.format(here=ikey), 'bar', 600, 600)
        canv.cd()
        ymax=max(h_ex0.GetMaximum(),0)#h_ex1.GetMaximum())
        h_ex0.GetYaxis().SetRangeUser(0.,ymax+ymax/2.0)
        #h_ex0.GetXaxis().SetTitle(ikey)
        #h_ex0.GetYaxis().SetTitle("a.u.")
        #h_ex0.GetXaxis().SetRangeUser(-2.5,2.5)
        #h_ex0.GetYaxis().SetNdivisions(505)
        #h_ex0.GetXaxis().SetNdivisions(510)
        h_ex0.GetYaxis().SetTitleOffset(1.2)
        h_ex0.GetXaxis().SetTitleOffset(1.0)
        h_ex0.GetYaxis().SetTitleSize(0.045)
        h_ex0.GetXaxis().SetTitleSize(0.045)


        h_ex0.Draw('hist')

        leg = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
        leg.SetTextSize(0.045); leg.SetShadowColor(0);        leg.SetFillStyle(0)
        leg.SetTextFont(42);   leg.SetBorderSize(0);
        #leg.SetNColumns(3)
        leg.AddEntry(h_ex0, 'XToYYPrime'  , 'l')
        leg.AddEntry('NULL', ival[5],'')
        leg.AddEntry('NULL', str(h_ex0.GetEntries())+" entries",'')
        # leg.AddEntry(h_ex1, 'LL'  , 'l')
        # leg.AddEntry(h_ex2, 'TLnLT'  , 'l')
        leg.SetLineColor(ROOT.kWhite)
        leg.SetFillColor(ROOT.kWhite)
        leg.Draw('same')
        # if "Gen_M_jj" in ikey:
        #     q = __import__("functools").partial(__import__("os")._exit, 0)  # FIXME
        #     __import__("IPython").embed()  # FIXME

        
        lat = ROOT.TLatex()
        lat.SetNDC()
        #fout.WriteTObject(canv)
        canv.SaveAs('{od}/{mww}.pdf'.format(od=outdir,mww=ikey))
        canv.SaveAs('{od}/{mww}.png'.format(od=outdir,mww=ikey))

###############################
def overlay2dshapes(outdir,l_ex0):
    # outdir='genplots'

    # basedir='./'
    # outdir='{Here}{when}{here}'.format(Here=basedir,when=date,here='gen')
    # if outdir not in os.listdir(basedir):
    #     os.system('mkdir -p {od}'.format(od=outdir))

    # f_ex0 = ROOT.TFile(f'{fname}.root', 'READ'); l_ex0 = f_ex0.Get('tree')
    # f_ex1 = ROOT.TFile('/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/lhetools/LatinoTreesLHE/WpmWmp_LL_genlevel.root', 'READ'); l_ex1 = f_ex1.Get('tree')
    # f_ex2 = ROOT.TFile('/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/lhetools/LatinoTreesLHE/WpmWmp_TLnLT_genlevel.root', 'READ'); l_ex2 = f_ex2.Get('tree')
    pdraw2d = {
        #####  name : [vnamex, vnamey, xtitle, ytitle, nXbins, nYbins, xmin, ymin, xmax, ymax, sel, Xbins, Ybins]
        # "delta_R_qq_YP25_vs_V_pt": ["delta_R_qq_YP25","V_pt","delta_R_qq_YP25","V_pt",200,200,-0.1,0,2,4000,'V_pdgId==25',[],[]],
        # "delta_R_qq_Y23_vs_V_pt": ["delta_R_qq_Y23","V_pt","delta_R_qq_Y23","V_pt",200,200,-0.1,0,2,4000,'V_pdgId==23',[],[]],
        # "delta_R_qq_YP25_vs_V_eta": ["delta_R_qq_YP25","V_eta","delta_R_qq_YP25","V_eta",50,50,-0.1,-4,2,4,'V_pdgId==25',[],[]],
        # "delta_R_qq_Y23_vs_V_eta": ["delta_R_qq_Y23","V_eta","delta_R_qq_Y23","V_eta",50,50,-0.1,-4,2,4,'V_pdgId==23',[],[]],
        # "delta_R_qq_YP25_vs_V_phi": ["delta_R_qq_YP25","V_phi","delta_R_qq_YP25","V_phi",50,50,-0.1,-4,2,4,'V_pdgId==25',[],[]],
        # "delta_R_qq_Y23_vs_V_phi": ["delta_R_qq_Y23","V_phi","delta_R_qq_Y23","V_phi",50,50,-0.1,-4,2,4,'V_pdgId==23',[],[]],
        "M_X_vs_M_YP_delta_R_high" : ["M_X", "M_YP", "M_X", "M_YP", 50, 50, 0, 0, 6500, 500, "delta_R_qq_YP25>=0.8 && V_pdgId==25",binsForX,binsForYP],
        "M_X_vs_M_Y_delta_R_high" : ["M_X", "M_Y", "M_X", "M_Y", 50, 50, 0, 0, 6500, 500, "delta_R_qq_Y23>=0.8 && V_pdgId==23",binsForX,binsForY],
        "M_X_vs_M_YP_delta_R_low" : ["M_X", "M_YP", "M_X", "M_YP", 50, 50, 0, 0, 6500, 500, "delta_R_qq_YP25<0.8 && V_pdgId==25",binsForX,binsForYP],
        "M_X_vs_M_Y_delta_R_low" : ["M_X", "M_Y", "M_X", "M_Y", 50, 50, 0, 0, 6500, 500, "delta_R_qq_Y23<0.8 && V_pdgId==23",binsForX,binsForY],
        "M_X_vs_M_Y" : ["M_X", "M_Y", "M_X", "M_Y", 50, 50, 0, 0, 6500, 500, "",binsForX,binsForY],
        "M_X_vs_M_YP" : ["M_X", "M_YP", "M_X", "M_YP", 50, 50, 0, 0, 6500, 500, "",binsForX,binsForYP],
        "M_Y_vs_M_YP" : ["M_Y", "M_YP", "M_Y", "M_YP", 50, 50, 0, 0, 500, 500, "",binsForY,binsForYP],
        # "M_X_vs_M_Y+M_YP" : ["M_X", "M_Y+M_YP", "M_X", "M_Y+M_YP", 50, 50, 0, 0, 6500, 1000, "",binsForX,[]],
        # "delta_R_qq_Y23_vs_delta_M" : ["delta_R_qq_Y23", "M_X-(M_Y+M_YP)", "delta_R_qq_Y23", "delta M", 200, 200, 0, 0, 2, 6500, "",[],[]],
        # "delta_R_qq_YP25_vs_delta_M" : ["delta_R_qq_YP25", "M_X-(M_Y+M_YP)", "delta_R_qq_YP25", "delta M", 200, 200, 0, 0, 2, 6500, "",[],[]],
        "delta_R_qq_Y23_vs_M_Y_ratio_log" : ["4*M_Y/M_X", "delta_R_qq_Y23", "4*M_Y_ratio", "delta_R_qq_Y23", 500, 200, 0, 0, 0.55, 2, "V_pdgId==23",logBinning_YX,[]],
        "delta_R_qq_YP25_vs_M_YP_ratio_log" : ["4*M_YP/M_X", "delta_R_qq_YP25", "4*M_YP_ratio" ,"delta_R_qq_YP25", 100, 200, 0, 0, 0.55, 2, "V_pdgId==25",logBinning_YX,[]],

        # # "Gen_M_jj_vs_M_jY" : ["Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "M_jY", 50, 50, 0, 0, 6500, 1000, "((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1)) || ((Gen_matched_Y[1]==1 && Gen_matched_Y[0]==0) && (Gen_matched_YP[1]==0 && Gen_matched_YP[0]==1))"],
        # # "Gen_M_jj_vs_M_jYP" : ["Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "M_jYP", 50, 50, 0, 0, 6500, 1000, "((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1)) || ((Gen_matched_YP[1]==1 && Gen_matched_YP[0]==0) && (Gen_matched_Y[1]==0 && Gen_matched_Y[0]==1))"],

        # "Reco_M_jj_vs_M_YP_delta_R_high" : ["Reco_M_jj", "M_YP", "Reco_M_jj", "M_YP", 50, 50, 0, 0, 6500, 500, "delta_R_qq_YP25>=0.8 && V_pdgId==25",binsForX,binsForYP],
        
        "GenJet0_deltaR_Y[0]_vs_Gen_Mjj"  : ["GenJet0_deltaR_Y[0]","Gen_M_jj","GenJet0_deltaR_Y[0]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet1_deltaR_Y[0]_vs_Gen_Mjj"  : ["GenJet1_deltaR_Y[0]","Gen_M_jj","GenJet1_deltaR_Y[0]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet0_deltaR_YP[0]_vs_Gen_Mjj" : ["GenJet0_deltaR_YP[0]","Gen_M_jj","GenJet0_deltaR_YP[0]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet1_deltaR_YP[0]_vs_Gen_Mjj" : ["GenJet1_deltaR_YP[0]","Gen_M_jj","GenJet1_deltaR_YP[0]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet0_deltaR_Y[0]_vs_Reco_Mjj"  : ["FatJet0_deltaR_Y[0]","Reco_M_jj","FatJet0_deltaR_Y[0]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet1_deltaR_Y[0]_vs_Reco_Mjj"  : ["FatJet1_deltaR_Y[0]","Reco_M_jj","FatJet1_deltaR_Y[0]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet0_deltaR_YP[0]_vs_Reco_Mjj" : ["FatJet0_deltaR_YP[0]","Reco_M_jj","FatJet0_deltaR_YP[0]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet1_deltaR_YP[0]_vs_Reco_Mjj" : ["FatJet1_deltaR_YP[0]","Reco_M_jj","FatJet1_deltaR_YP[0]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet0_deltaR_Y[1]_vs_Gen_Mjj"  : ["GenJet0_deltaR_Y[1]","Gen_M_jj","GenJet0_deltaR_Y[1]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet1_deltaR_Y[1]_vs_Gen_Mjj"  : ["GenJet1_deltaR_Y[1]","Gen_M_jj","GenJet1_deltaR_Y[1]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet0_deltaR_YP[1]_vs_Gen_Mjj" : ["GenJet0_deltaR_YP[1]","Gen_M_jj","GenJet0_deltaR_YP[1]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "GenJet1_deltaR_YP[1]_vs_Gen_Mjj" : ["GenJet1_deltaR_YP[1]","Gen_M_jj","GenJet1_deltaR_YP[1]","Gen_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet0_deltaR_Y[1]_vs_Reco_Mjj"  : ["FatJet0_deltaR_Y[1]","Reco_M_jj","FatJet0_deltaR_Y[1]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet1_deltaR_Y[1]_vs_Reco_Mjj"  : ["FatJet1_deltaR_Y[1]","Reco_M_jj","FatJet1_deltaR_Y[1]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet0_deltaR_YP[1]_vs_Reco_Mjj" : ["FatJet0_deltaR_YP[1]","Reco_M_jj","FatJet0_deltaR_YP[1]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],
        "FatJet1_deltaR_YP[1]_vs_Reco_Mjj" : ["FatJet1_deltaR_YP[1]","Reco_M_jj","FatJet1_deltaR_YP[1]","Reco_M_jj",100,50,-0.1,0,2,6500,"",[],[]],

        }

    for ikey,ival in pdraw2d.items():
        print(ikey)
        #tree, selcuts, vnamex, vnamey, xtitle, ytitle, nXbins, nYbins, xmin, ymin, xmax, ymax, color = 1
        h_ex0 = get2dhist(l_ex0,ival[10],ival[0],ival[1],ival[2],ival[3],ival[4],ival[5],ival[6],ival[7],ival[8],ival[9],ival[11],ival[12])
        print(ikey,h_ex0.GetEntries())

        ROOT.gROOT.SetBatch()
        ROOT.gStyle.SetOptStat(0)
        canv = ROOT.TCanvas('{here}'.format(here=ikey), 'bar', 600, 600)
        canv.cd()
        ymax=max(h_ex0.GetMaximum(),0)#h_ex1.GetMaximum())
        # h_ex0.GetYaxis().SetRangeUser(0.,ymax+ymax/2.0)
        #h_ex0.GetXaxis().SetTitle(ikey)
        #h_ex0.GetYaxis().SetTitle("a.u.")
        #h_ex0.GetXaxis().SetRangeUser(-2.5,2.5)
        #h_ex0.GetYaxis().SetNdivisions(505)
        #h_ex0.GetXaxis().SetNdivisions(510)
        h_ex0.GetYaxis().SetTitleOffset(1.0)
        h_ex0.GetXaxis().SetTitleOffset(1.0)
        h_ex0.GetYaxis().SetTitleSize(0.045)
        h_ex0.GetXaxis().SetTitleSize(0.045)

        # if 'log' in ikey:
        #     print("log plot found")
        #     ROOT.gPad.SetLogx(1)
        # else:
        #     ROOT.gPad.SetLogx(0)

        h_ex0.Draw('colz')

        leg = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)

        leg.SetTextSize(0.045); leg.SetShadowColor(0);        leg.SetFillStyle(0)
        leg.SetTextFont(42);   leg.SetBorderSize(0);
        #leg.SetNColumns(3)
        leg.AddEntry(h_ex0, 'XToYYPrime'  , 'l')
        if "ratio" in ikey:
            dR1 = ROOT.TF1("dR1","x",logBinning_YX[0],logBinning_YX[-1])
            dR1.SetLineColor(46)
            dR1.Draw("SAME")
            leg.AddEntry(dR1,'y=x','l')
        leg.AddEntry('NULL', ival[10],'')

        leg.SetLineColor(ROOT.kWhite)
        leg.SetFillColor(ROOT.kWhite)
        leg.Draw('same')
        
        lat = ROOT.TLatex()
        lat.SetNDC()
        #fout.WriteTObject(canv)
        canv.SaveAs('{od}/{mww}.pdf'.format(od=outdir,mww=ikey))
        canv.SaveAs('{od}/{mww}.png'.format(od=outdir,mww=ikey))



def overlayshapesSpecial(outdir,l_ex0):
    # outdir='genplots'

    # basedir='./'
    # outdir='{Here}{when}{here}'.format(Here=basedir,when=date,here='gen')
    # if outdir not in os.listdir(basedir):
    #     os.system('mkdir -p {od}'.format(od=outdir))

    # f_ex0 = ROOT.TFile(f'{fname}.root', 'READ'); l_ex0 = f_ex0.Get('tree')




    # special plot

    # List of histograms to draw
    #         [[hist,title,subtitle/sel,drawArg], ...]
    s_hists = []

    s_hist_sel_M_Y = get2dhist(l_ex0, sel_str_Gen_matching_paper, "M_X", "M_Y", "M_X", "M_Y", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    s_hist_M_Y = get2dhist(l_ex0, "","M_X", "M_Y", "M_X", "M_Y", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    s_hist_sel_M_Y.Divide(s_hist_M_Y)
    s_hists.append([s_hist_sel_M_Y,'M_sel_div_M__Y','','colz'])

    s_hist_sel_M_YP = get2dhist(l_ex0, sel_str_Gen_matching_paper, "M_X", "M_YP", "M_X", "M_YP", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    s_hist_M_YP = get2dhist(l_ex0, "","M_X", "M_YP", "M_X", "M_YP", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    s_hist_sel_M_YP.Divide(s_hist_M_YP)
    s_hists.append([s_hist_sel_M_YP,'M_sel_div_M__YP','','colz'])  


    s_hist_Gen_M_jj = gethist(l_ex0,sel_str_Gen_matching_paper,"Gen_M_jj","",326,0,7000,ROOT.kGreen+2)
    s_hist_Gen_M_jj_sel = gethist(l_ex0,sel_str_Gen_matching_paper,"Gen_M_jj","Gen_M_jj_match_paper",326,0,7000,ROOT.kGreen+2)

    s_hist_Gen_M_j0Y = gethist(l_ex0,"((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1))","GenJetAK8_mass[0]","Gen_M_jY",100,0,1000,ROOT.kGreen+2)
    s_hist_Gen_M_j1Y = gethist(l_ex0,"((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1))","GenJetAK8_mass[1]","Gen_M_j1Y",100,0,1000,ROOT.kGreen+2)
    s_hist_Gen_M_j0Y.Add(s_hist_Gen_M_j1Y)
    s_hists.append([s_hist_Gen_M_j0Y,'Gen_M_jY','','hist'])
    
    s_hist_Gen_M_j0YP = gethist(l_ex0,"((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1))","GenJetAK8_mass[0]","Gen_M_jYP",100,0,1000,ROOT.kGreen+2)
    s_hist_Gen_M_j1YP = gethist(l_ex0,"((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1))","GenJetAK8_mass[1]","Gen_M_j1YP",100,0,1000,ROOT.kGreen+2)
    s_hist_Gen_M_j0YP.Add(s_hist_Gen_M_j1YP)
    s_hists.append([s_hist_Gen_M_j0YP,'Gen_M_jYP','','hist'])


    s_hist_Gen_sel_M_Y = get2dhist(l_ex0, f"{sel_str_Gen_matching_paper} && ((Gen_matched_Y[Iteration$]==1) && (Gen_matched_YP[abs(Iteration$-1)]==1))", "Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "GenJetAK8_mass", 326,100,0,0,7000,1000, [], [])
    s_hist_Gen_M_Y = get2dhist(l_ex0, "((Gen_matched_Y[Iteration$]==1) && (Gen_matched_YP[abs(Iteration$-1)]==1))","Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "GenJetAK8_mass", 326,100,0,0,7000,1000, [], [])
    s_hist_Gen_sel_M_YP = get2dhist(l_ex0, f"{sel_str_Gen_matching_paper} && ((Gen_matched_YP[Iteration$]==1) && (Gen_matched_Y[abs(Iteration$-1)]==1))", "Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "GenJetAK8_mass", 326,100,0,0,7000,1000, [], [])
    s_hist_Gen_M_YP = get2dhist(l_ex0, "((Gen_matched_YP[Iteration$]==1) && (Gen_matched_Y[abs(Iteration$-1)]==1))","Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "GenJetAK8_mass", 326,100,0,0,7000,1000, [], [])

    s_hist_Gen_Mjj_vs_MjY = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Gen_M_jj.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Gen_M_j0Y.GetNbinsX() + 1):
            content_x = s_hist_Gen_M_jj.GetBinContent(bin_x)
            content_y = s_hist_Gen_M_j0Y.GetBinContent(bin_y)
            bin_center_x = s_hist_Gen_M_jj.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Gen_M_j0Y.GetXaxis().GetBinCenter(bin_y)
            s_hist_Gen_Mjj_vs_MjY.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Gen_Mjj_vs_MjY.GetXaxis().SetTitle("Gen_M_jj")
    s_hist_Gen_Mjj_vs_MjY.GetYaxis().SetTitle("Gen_M_jY")
    s_hists.append([s_hist_Gen_Mjj_vs_MjY, 'Gen_Mjj_vs_MjY', 'no selection', 'colz'])

    s_hist_Gen_Mjj_vs_MjY_sel = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Gen_M_jj_sel.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Gen_M_j0Y.GetNbinsX() + 1):
            content_x = s_hist_Gen_M_jj_sel.GetBinContent(bin_x)
            content_y = s_hist_Gen_M_j0Y.GetBinContent(bin_y)
            bin_center_x = s_hist_Gen_M_jj_sel.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Gen_M_j0Y.GetXaxis().GetBinCenter(bin_y)
            s_hist_Gen_Mjj_vs_MjY_sel.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Gen_Mjj_vs_MjY_sel.GetXaxis().SetTitle("Gen_M_jj")
    s_hist_Gen_Mjj_vs_MjY_sel.GetYaxis().SetTitle("Gen_M_jY")
    s_hists.append([s_hist_Gen_Mjj_vs_MjY_sel, 'Gen_Mjj_vs_MjY_sel', 'matching & VV/VH paper', 'colz'])

    s_hist_Gen_Mjj_vs_MjYP = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Gen_M_jj.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Gen_M_j0YP.GetNbinsX() + 1):
            content_x = s_hist_Gen_M_jj.GetBinContent(bin_x)
            content_y = s_hist_Gen_M_j0YP.GetBinContent(bin_y)
            bin_center_x = s_hist_Gen_M_jj.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Gen_M_j0YP.GetXaxis().GetBinCenter(bin_y)
            s_hist_Gen_Mjj_vs_MjYP.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Gen_Mjj_vs_MjYP.GetXaxis().SetTitle("Gen_M_jj")
    s_hist_Gen_Mjj_vs_MjYP.GetYaxis().SetTitle("Gen_M_jYP")
    s_hists.append([s_hist_Gen_Mjj_vs_MjYP, 'Gen_Mjj_vs_MjYP', 'no selection', 'colz'])

    s_hist_Gen_Mjj_vs_MjYP_sel = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Gen_M_jj_sel.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Gen_M_j0YP.GetNbinsX() + 1):
            content_x = s_hist_Gen_M_jj_sel.GetBinContent(bin_x)
            content_y = s_hist_Gen_M_j0YP.GetBinContent(bin_y)
            bin_center_x = s_hist_Gen_M_jj_sel.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Gen_M_j0YP.GetXaxis().GetBinCenter(bin_y)
            s_hist_Gen_Mjj_vs_MjYP_sel.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Gen_Mjj_vs_MjYP_sel.GetXaxis().SetTitle("Gen_M_jj")
    s_hist_Gen_Mjj_vs_MjYP_sel.GetYaxis().SetTitle("Gen_M_jYP")
    s_hists.append([s_hist_Gen_Mjj_vs_MjYP_sel, 'Gen_Mjj_vs_MjYP_sel', 'matching & VV/VH paper', 'colz'])

    s_hist_Gen_sel_M_Y.Divide(s_hist_Gen_M_Y)
    s_hists.append([s_hist_Gen_sel_M_Y,'Gen_M_sel_div_M__Y','','colz'])

    s_hist_Gen_sel_M_YP.Divide(s_hist_Gen_M_YP)
    s_hists.append([s_hist_Gen_sel_M_YP,'Gen_M_sel_div_M__YP','','colz'])


    # s_hist_Reco_M_j0Y = gethist(l_ex0,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1))","FatJet_mass[0]","Reco_M_jY",100,0,1000,ROOT.kGreen+2)
    # s_hist_Reco_M_j1Y = gethist(l_ex0,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1))","FatJet_mass[1]","Reco_M_j1Y",100,0,1000,ROOT.kGreen+2)
    # s_hist_Reco_M_j0Y.Add(s_hist_Reco_M_j1Y)
    # s_hists.append([s_hist_Reco_M_j0Y,'Reco_M_jY','','hist'])
    
    # s_hist_Reco_M_j0YP = gethist(l_ex0,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1))","FatJet_mass[0]","Reco_M_jYP",100,0,1000,ROOT.kGreen+2)
    # s_hist_Reco_M_j1YP = gethist(l_ex0,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1))","FatJet_mass[1]","Reco_M_j1YP",100,0,1000,ROOT.kGreen+2)
    # s_hist_Reco_M_j0YP.Add(s_hist_Reco_M_j1YP)
    # s_hists.append([s_hist_Reco_M_j0YP,'Reco_M_jYP','','hist'])


    # s_hist3 = get2dhist(l_ex0, sel_str_Reco_matching_paper, "M_X", "M_Y", "M_X", "M_Y", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    # s_hist4 = get2dhist(l_ex0, "","M_X", "M_Y", "M_X", "M_Y", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    # s_hist3.Divide(s_hist4)
    # s_hists.append([s_hist3,'Reco_M_sel_div_M__Y','','colz'])

    # s_hist4 = get2dhist(l_ex0, sel_str_Reco_matching_paper, "M_X", "M_YP", "M_X", "M_YP", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    # s_hist5 = get2dhist(l_ex0, "","M_X", "M_YP", "M_X", "M_YP", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    # s_hist4.Divide(s_hist5)
    # s_hists.append([s_hist4,'Reco_M_sel_div_M__YP','','colz'])  
    s_hist_Reco_M_jj = gethist(l_ex0,sel_str_Reco_matching_paper,"Reco_M_jj","",326,0,7000,ROOT.kGreen+2)
    s_hist_Reco_M_jj_sel = gethist(l_ex0,sel_str_Reco_matching_paper,"Reco_M_jj","Reco_M_jj_match_paper",326,0,7000,ROOT.kGreen+2)

    s_hist_Reco_M_j0Y = gethist(l_ex0,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1))","FatJet_mass[0]","Reco_M_jY",100,0,1000,ROOT.kGreen+2)
    s_hist_Reco_M_j1Y = gethist(l_ex0,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1))","FatJet_mass[1]","Reco_M_j1Y",100,0,1000,ROOT.kGreen+2)
    s_hist_Reco_M_j0Y.Add(s_hist_Reco_M_j1Y)
    s_hists.append([s_hist_Reco_M_j0Y,'Reco_M_jY','','hist'])
    
    s_hist_Reco_M_j0YP = gethist(l_ex0,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1))","FatJet_mass[0]","Reco_M_jYP",100,0,1000,ROOT.kGreen+2)
    s_hist_Reco_M_j1YP = gethist(l_ex0,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1))","FatJet_mass[1]","Reco_M_j1YP",100,0,1000,ROOT.kGreen+2)
    s_hist_Reco_M_j0YP.Add(s_hist_Reco_M_j1YP)
    s_hists.append([s_hist_Reco_M_j0YP,'Reco_M_jYP','','hist'])


    s_hist_Reco_sel_M_Y = get2dhist(l_ex0, f"{sel_str_Reco_matching_paper} && ((Reco_matched_Y[Iteration$]==1) && (Reco_matched_YP[abs(Iteration$-1)]==1))", "Reco_M_jj", "FatJet_mass", "Reco_M_jj", "FatJet_mass", 326,100,0,0,7000,1000, [], [])
    s_hist_Reco_M_Y = get2dhist(l_ex0, "((Reco_matched_Y[Iteration$]==1) && (Reco_matched_YP[abs(Iteration$-1)]==1))","Reco_M_jj", "FatJet_mass", "Reco_M_jj", "FatJet_mass", 326,100,0,0,7000,1000, [], [])
    s_hist_Reco_sel_M_YP = get2dhist(l_ex0, f"{sel_str_Reco_matching_paper} && ((Reco_matched_YP[Iteration$]==1) && (Reco_matched_Y[abs(Iteration$-1)]==1))", "Reco_M_jj", "FatJet_mass", "Reco_M_jj", "FatJet_mass", 326,100,0,0,7000,1000, [], [])
    s_hist_Reco_M_YP = get2dhist(l_ex0, "((Reco_matched_YP[Iteration$]==1) && (Reco_matched_Y[abs(Iteration$-1)]==1))","Reco_M_jj", "FatJet_mass", "Reco_M_jj", "FatJet_mass", 326,100,0,0,7000,1000, [], [])

    s_hist_Reco_Mjj_vs_MjY = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Reco_M_jj.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Reco_M_j0Y.GetNbinsX() + 1):
            content_x = s_hist_Reco_M_jj.GetBinContent(bin_x)
            content_y = s_hist_Reco_M_j0Y.GetBinContent(bin_y)
            bin_center_x = s_hist_Reco_M_jj.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Reco_M_j0Y.GetXaxis().GetBinCenter(bin_y)
            s_hist_Reco_Mjj_vs_MjY.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Reco_Mjj_vs_MjY.GetXaxis().SetTitle("Reco_M_jj")
    s_hist_Reco_Mjj_vs_MjY.GetYaxis().SetTitle("Reco_M_jY")
    s_hists.append([s_hist_Reco_Mjj_vs_MjY, 'Reco_Mjj_vs_MjY', 'no selection', 'colz'])

    s_hist_Reco_Mjj_vs_MjY_sel = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Reco_M_jj_sel.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Reco_M_j0Y.GetNbinsX() + 1):
            content_x = s_hist_Reco_M_jj_sel.GetBinContent(bin_x)
            content_y = s_hist_Reco_M_j0Y.GetBinContent(bin_y)
            bin_center_x = s_hist_Reco_M_jj_sel.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Reco_M_j0Y.GetXaxis().GetBinCenter(bin_y)
            s_hist_Reco_Mjj_vs_MjY_sel.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Reco_Mjj_vs_MjY_sel.GetXaxis().SetTitle("Reco_M_jj")
    s_hist_Reco_Mjj_vs_MjY_sel.GetYaxis().SetTitle("Reco_M_jY")
    s_hists.append([s_hist_Reco_Mjj_vs_MjY_sel, 'Reco_Mjj_vs_MjY_sel', 'matching & VV/VH paper', 'colz'])

    s_hist_Reco_Mjj_vs_MjYP = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Reco_M_jj.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Reco_M_j0YP.GetNbinsX() + 1):
            content_x = s_hist_Reco_M_jj.GetBinContent(bin_x)
            content_y = s_hist_Reco_M_j0YP.GetBinContent(bin_y)
            bin_center_x = s_hist_Reco_M_jj.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Reco_M_j0YP.GetXaxis().GetBinCenter(bin_y)
            s_hist_Reco_Mjj_vs_MjYP.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Reco_Mjj_vs_MjYP.GetXaxis().SetTitle("Reco_M_jj")
    s_hist_Reco_Mjj_vs_MjYP.GetYaxis().SetTitle("Reco_M_jYP")
    s_hists.append([s_hist_Reco_Mjj_vs_MjYP, 'Reco_Mjj_vs_MjYP', 'no selection', 'colz'])

    s_hist_Reco_Mjj_vs_MjYP_sel = ROOT.TH2F("hist2d" ,'',326,0,7000,100,0,1000)
    for bin_x in range(1, s_hist_Reco_M_jj_sel.GetNbinsX() + 1):
        for bin_y in range(1, s_hist_Reco_M_j0YP.GetNbinsX() + 1):
            content_x = s_hist_Reco_M_jj_sel.GetBinContent(bin_x)
            content_y = s_hist_Reco_M_j0YP.GetBinContent(bin_y)
            bin_center_x = s_hist_Reco_M_jj_sel.GetXaxis().GetBinCenter(bin_x)
            bin_center_y = s_hist_Reco_M_j0YP.GetXaxis().GetBinCenter(bin_y)
            s_hist_Reco_Mjj_vs_MjYP_sel.Fill(bin_center_x, bin_center_y, content_x * content_y)
    s_hist_Reco_Mjj_vs_MjYP_sel.GetXaxis().SetTitle("Reco_M_jj")
    s_hist_Reco_Mjj_vs_MjYP_sel.GetYaxis().SetTitle("Reco_M_jYP")
    s_hists.append([s_hist_Reco_Mjj_vs_MjYP_sel, 'Reco_Mjj_vs_MjYP_sel', 'matching & VV/VH paper', 'colz'])

    s_hist_Reco_sel_M_Y.Divide(s_hist_Reco_M_Y)
    s_hists.append([s_hist_Reco_sel_M_Y,'Reco_M_sel_div_M__Y','','colz'])

    s_hist_Reco_sel_M_YP.Divide(s_hist_Reco_M_YP)
    s_hists.append([s_hist_Reco_sel_M_YP,'Reco_M_sel_div_M__YP','','colz'])





    for s_hist in s_hists:
        print(s_hist[1])
        ROOT.gROOT.SetBatch()
        ROOT.gStyle.SetOptStat(0)
        canv = ROOT.TCanvas(s_hist[1], 'bar', 600, 600)
        canv.cd()
        s_hist[0].Draw(s_hist[3])
        leg = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
        leg.SetTextSize(0.045); leg.SetShadowColor(0);        leg.SetFillStyle(0)
        leg.SetTextFont(42);   leg.SetBorderSize(0);
        leg.AddEntry(s_hist[0], 'XToYYPrime'  , 'l')
        leg.AddEntry('NULL', s_hist[2],'')
        leg.SetLineColor(ROOT.kWhite)
        leg.SetFillColor(ROOT.kWhite)
        s_hist[0].GetYaxis().SetTitleOffset(1.0)
        s_hist[0].GetXaxis().SetTitleOffset(1.0)
        s_hist[0].GetYaxis().SetTitleSize(0.045)
        s_hist[0].GetXaxis().SetTitleSize(0.045)
        leg.Draw('same')

        lat = ROOT.TLatex()
        lat.SetNDC()
        #fout.WriteTObject(canv)
        canv.SaveAs('{od}/{title}.pdf'.format(od=outdir,title=s_hist[1]))
        canv.SaveAs('{od}/{title}.png'.format(od=outdir,title=s_hist[1]))


def overlayshapesQCD(outdir,l_ex0,l_ex1):
    pdraw={
        #### xtitle : [xname, nXbins, xmin, xmax, sel0, sel1, _info]
              "GenJetAK8_pt[0]" : [  "GenJetAK8_pt[0]",100 ,0   ,4000,"","","no selection"],
             "GenJetAK8_eta[0]" : [ "GenJetAK8_eta[0]",100 ,-3.5,3.5 ,"","","no selection"],
             "GenJetAK8_phi[0]" : [ "GenJetAK8_phi[0]",100 ,-3.5,3.5 ,"","","no selection"],
            "GenJetAK8_mass[0]" : ["GenJetAK8_mass[0]",100 ,0   ,700 ,"","","no selection"],
          "GenJetAK8_pt_sel[0]" : [  "GenJetAK8_pt[0]",100 ,0   ,4000,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],
         "GenJetAK8_eta_sel[0]" : [ "GenJetAK8_eta[0]",100 ,-3.5,3.5 ,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],
         "GenJetAK8_phi_sel[0]" : [ "GenJetAK8_phi[0]",100 ,-3.5,3.5 ,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],
        "GenJetAK8_mass_sel[0]" : ["GenJetAK8_mass[0]",100 ,0   ,700 ,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],

              "GenJetAK8_pt[1]" : [  "GenJetAK8_pt[1]",100 ,0   ,4000,"","","no selection"],
             "GenJetAK8_eta[1]" : [ "GenJetAK8_eta[1]",100 ,-3.5,3.5 ,"","","no selection"],
             "GenJetAK8_phi[1]" : [ "GenJetAK8_phi[1]",100 ,-3.5,3.5 ,"","","no selection"],
            "GenJetAK8_mass[1]" : ["GenJetAK8_mass[1]",100 ,0   ,700 ,"","","no selection"],
          "GenJetAK8_pt_sel[1]" : [  "GenJetAK8_pt[1]",100 ,0   ,4000,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],
         "GenJetAK8_eta_sel[1]" : [ "GenJetAK8_eta[1]",100 ,-3.5,3.5 ,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],
         "GenJetAK8_phi_sel[1]" : [ "GenJetAK8_phi[1]",100 ,-3.5,3.5 ,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],
        "GenJetAK8_mass_sel[1]" : ["GenJetAK8_mass[1]",100 ,0   ,700 ,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],

            "Gen_M_jj" : ["Gen_M_jj",326 ,0   ,7000,"","","no selection"],
        "Gen_M_jj_sel" : ["Gen_M_jj",326 ,0   ,7000,sel_str_Gen_paper,sel_str_Gen_paper,"paper sel"],

 
                   "FatJet_pt[0]" : [       "FatJet_pt[0]",100 ,0   ,4000,"","","no selection"],
                  "FatJet_eta[0]" : [      "FatJet_eta[0]",100 ,-3.5,3.5 ,"","","no selection"],
                  "FatJet_phi[0]" : [      "FatJet_phi[0]",100 ,-3.5,3.5 ,"","","no selection"],
                 "FatJet_mass[0]" : [     "FatJet_mass[0]",100 ,0   ,700 ,"","","no selection"],
            "FatJet_msoftdrop[0]" : ["FatJet_msoftdrop[0]",100 ,0   ,700 ,"","","no selection"],
               "FatJet_pt_sel[0]" : [       "FatJet_pt[0]",100 ,0   ,4000,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
              "FatJet_eta_sel[0]" : [      "FatJet_eta[0]",100 ,-3.5,3.5 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
              "FatJet_phi_sel[0]" : [      "FatJet_phi[0]",100 ,-3.5,3.5 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
             "FatJet_mass_sel[0]" : [     "FatJet_mass[0]",100 ,0   ,700 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_msoftdrop_sel[0]" : ["FatJet_msoftdrop[0]",100 ,0   ,700 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],

                   "FatJet_pt[1]" : [       "FatJet_pt[1]",100 ,0   ,4000,"","","no selection"],
                  "FatJet_eta[1]" : [      "FatJet_eta[1]",100 ,-3.5,3.5 ,"","","no selection"],
                  "FatJet_phi[1]" : [      "FatJet_phi[1]",100 ,-3.5,3.5 ,"","","no selection"],
                 "FatJet_mass[1]" : [     "FatJet_mass[1]",100 ,0   ,700 ,"","","no selection"],
            "FatJet_msoftdrop[1]" : ["FatJet_msoftdrop[1]",100 ,0   ,700 ,"","","no selection"],
               "FatJet_pt_sel[1]" : [       "FatJet_pt[1]",100 ,0   ,4000,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
              "FatJet_eta_sel[1]" : [      "FatJet_eta[1]",100 ,-3.5,3.5 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
              "FatJet_phi_sel[1]" : [      "FatJet_phi[1]",100 ,-3.5,3.5 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
             "FatJet_mass_sel[1]" : [     "FatJet_mass[1]",100 ,0   ,700 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_msoftdrop_sel[1]" : ["FatJet_msoftdrop[1]",100 ,0   ,700 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],

            "Reco_M_jj" : ["Reco_M_jj",326 ,0   ,7000,"","","no selection"],
        "Reco_M_jj_sel" : ["Reco_M_jj",326 ,0   ,7000,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],


               "nFatJet" : ["nFatJet"       ,6   ,0   ,5   ,"","","no selection"],
        "FatJet_n2b1[0]" : ["FatJet_n2b1[0]",100 ,0   ,0.6 ,"","","no selection"],
        "FatJet_n3b1[0]" : ["FatJet_n3b1[0]",100 ,0   ,5   ,"","","no selection"],
        "FatJet_tau1[0]" : ["FatJet_tau1[0]",100 ,1e-3,0.7 ,"","","no selection"],
        "FatJet_tau2[0]" : ["FatJet_tau2[0]",100 ,1e-3,0.3 ,"","","no selection"],
        "FatJet_tau3[0]" : ["FatJet_tau3[0]",100 ,1e-3,0.2 ,"","","no selection"],
        "FatJet_tau4[0]" : ["FatJet_tau4[0]",100 ,1e-4,0.2 ,"","","no selection"],

        "FatJet_n2b1[1]" : ["FatJet_n2b1[1]",100 ,0   ,0.6 ,"","","no selection"],
        "FatJet_n3b1[1]" : ["FatJet_n3b1[1]",100 ,0   ,5   ,"","","no selection"],
        "FatJet_tau1[1]" : ["FatJet_tau1[1]",100 ,1e-3,0.7 ,"","","no selection"],
        "FatJet_tau2[1]" : ["FatJet_tau2[1]",100 ,1e-3,0.3 ,"","","no selection"],
        "FatJet_tau3[1]" : ["FatJet_tau3[1]",100 ,1e-3,0.2 ,"","","no selection"],
        "FatJet_tau4[1]" : ["FatJet_tau4[1]",100 ,1e-4,0.2 ,"","","no selection"],

               "nSubJet" : ["nSubJet"       ,6   ,0   ,5   ,"","","no selection"],
        "SubJet_mass[0]" : ["SubJet_mass[0]",100 ,1e-7,140 ,"","","no selection"],
        "SubJet_n2b1[0]" : ["SubJet_n2b1[0]",100 ,0   ,0.5 ,"","","no selection"],
        "SubJet_n3b1[0]" : ["SubJet_n3b1[0]",100 ,0   ,8   ,"","","no selection"],
        "SubJet_tau1[0]" : ["SubJet_tau1[0]",100 ,1e-5,0.3 ,"","","no selection"],
        "SubJet_tau2[0]" : ["SubJet_tau2[0]",100 ,1e-6,0.3 ,"","","no selection"],
        "SubJet_tau3[0]" : ["SubJet_tau3[0]",100 ,1e-6,0.3 ,"","","no selection"],
        "SubJet_tau4[0]" : ["SubJet_tau4[0]",100 ,1e-6,0.3 ,"","","no selection"],
        
        "SubJet_mass[1]" : ["SubJet_mass[1]",100 ,1e-7,140 ,"","","no selection"],
        "SubJet_n2b1[1]" : ["SubJet_n2b1[1]",100 ,0   ,0.5 ,"","","no selection"],
        "SubJet_n3b1[1]" : ["SubJet_n3b1[1]",100 ,0   ,8   ,"","","no selection"],
        "SubJet_tau1[1]" : ["SubJet_tau1[1]",100 ,1e-5,0.3 ,"","","no selection"],
        "SubJet_tau2[1]" : ["SubJet_tau2[1]",100 ,1e-6,0.3 ,"","","no selection"],
        "SubJet_tau3[1]" : ["SubJet_tau3[1]",100 ,1e-6,0.3 ,"","","no selection"],
        "SubJet_tau4[1]" : ["SubJet_tau4[1]",100 ,1e-6,0.3 ,"","","no selection"],


               "nFatJet_sel" : ["nFatJet"       ,6   ,0   ,5   ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_n2b1_sel[0]" : ["FatJet_n2b1[0]",100 ,0   ,0.6 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_n3b1_sel[0]" : ["FatJet_n3b1[0]",100 ,0   ,5   ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau1_sel[0]" : ["FatJet_tau1[0]",100 ,1e-3,0.7 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau2_sel[0]" : ["FatJet_tau2[0]",100 ,1e-3,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau3_sel[0]" : ["FatJet_tau3[0]",100 ,1e-3,0.2 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau4_sel[0]" : ["FatJet_tau4[0]",100 ,1e-4,0.2 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],

        "FatJet_n2b1_sel[1]" : ["FatJet_n2b1[1]",100 ,0   ,0.6 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_n3b1_sel[1]" : ["FatJet_n3b1[1]",100 ,0   ,5   ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau1_sel[1]" : ["FatJet_tau1[1]",100 ,1e-3,0.7 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau2_sel[1]" : ["FatJet_tau2[1]",100 ,1e-3,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau3_sel[1]" : ["FatJet_tau3[1]",100 ,1e-3,0.2 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "FatJet_tau4_sel[1]" : ["FatJet_tau4[1]",100 ,1e-4,0.2 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],

               "nSubJet_sel" : ["nSubJet"       ,6   ,0   ,5   ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_mass_sel[0]" : ["SubJet_mass[0]",100 ,1e-7,140 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_n2b1_sel[0]" : ["SubJet_n2b1[0]",100 ,0   ,0.5 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_n3b1_sel[0]" : ["SubJet_n3b1[0]",100 ,0   ,8   ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau1_sel[0]" : ["SubJet_tau1[0]",100 ,1e-5,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau2_sel[0]" : ["SubJet_tau2[0]",100 ,1e-6,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau3_sel[0]" : ["SubJet_tau3[0]",100 ,1e-6,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau4_sel[0]" : ["SubJet_tau4[0]",100 ,1e-6,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],

        "SubJet_mass_sel[1]" : ["SubJet_mass[1]",100 ,1e-7,140 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_n2b1_sel[1]" : ["SubJet_n2b1[1]",100 ,0   ,0.5 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_n3b1_sel[1]" : ["SubJet_n3b1[1]",100 ,0   ,8   ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau1_sel[1]" : ["SubJet_tau1[1]",100 ,1e-5,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau2_sel[1]" : ["SubJet_tau2[1]",100 ,1e-6,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau3_sel[1]" : ["SubJet_tau3[1]",100 ,1e-6,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        "SubJet_tau4_sel[1]" : ["SubJet_tau4[1]",100 ,1e-6,0.3 ,sel_str_Reco_paper,sel_str_Reco_paper,"paper sel"],
        }
    for ikey,ival in pdraw.items():
        print(f'{ikey} with bkg')
        ROOT.gROOT.SetBatch()
        ROOT.gStyle.SetOptStat(0)
        canv = ROOT.TCanvas(f"{ikey}_with_bkg", 'bar', 600, 600)
        canv.cd()
        if ('SubJet_mass' in ikey) or ('tau' in ikey):
            print("log plot found")
            ROOT.gPad.SetLogx(1)
            logBins = np.logspace(np.log10(ival[2]) if ival[2]>0 else np.log10(1e-7),np.log10(ival[3])+1,ival[1])
            #               tree,seli,xname,xtitle,nXbins,xmin,xmax, color=1, Xbins=[]
            h_ex0 = gethist(l_ex0,ival[4],ival[0],f"{ikey}_with_bkg",ival[1],ival[2],ival[3],ROOT.kAzure+1, logBins)
            h_ex1 = gethist(l_ex1,ival[5],ival[0],f"{ikey}_with_bkg",ival[1],ival[2],ival[3],ROOT.kOrange+5, logBins)
        else:
            ROOT.gPad.SetLogx(0)
            h_ex0 = gethist(l_ex0,ival[4],ival[0],f"{ikey}_with_bkg",ival[1],ival[2],ival[3],ROOT.kAzure+1)
            h_ex1 = gethist(l_ex1,ival[5],ival[0],f"{ikey}_with_bkg",ival[1],ival[2],ival[3],ROOT.kOrange+5)

        h_ex1.Draw('hist')
        h_integral_1 = h_ex1.Integral()
        h_ex1.Scale(1/h_integral_1 if h_integral_1 > 0 else 0)
        h_ex1.SetFillColorAlpha(ROOT.kOrange+5, 0.2)
        h_ex0.Draw('histsame')
        h_integral_0 = h_ex0.Integral()
        h_ex0.Scale(1/h_integral_0 if h_integral_0 > 0 else 0)
        h_ex0.SetFillColorAlpha(ROOT.kAzure+1, 0.5)
        ymax=max(h_ex0.GetMaximum(),h_ex1.GetMaximum(),0)
        h_ex0.GetYaxis().SetRangeUser(0.,ymax+ymax*0.1)
        h_ex1.GetYaxis().SetRangeUser(0.,ymax+ymax*0.1)
        h_ex0.GetYaxis().SetTitleOffset(1)
        h_ex1.GetYaxis().SetTitleOffset(1)

        leg = ROOT.TLegend(0.65, 0.78, 0.9, 0.9)
        leg.SetBorderSize(0); leg.SetShadowColor(0); leg.SetFillStyle(0)
        leg.AddEntry(h_ex0, 'Signal', 'l')
        leg.AddEntry(h_ex1, 'Background', 'l')
        leg.AddEntry('NULL', ival[6],'')
        leg.AddEntry('NULL', str(h_ex0.GetEntries())+" entries",'')
        leg.Draw('same')

        canv.SaveAs('{od}/{mww}_with_bkg.pdf'.format(od=outdir,mww=ikey))
        canv.SaveAs('{od}/{mww}_with_bkg.png'.format(od=outdir,mww=ikey))


    
##############################################


if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('--overlayshapes', action='store_true', dest='overlayshapes'        , default=True, help='blah blah ')
    global opts
    (opts, args) = parser.parse_args()


    if opts.overlayshapes:
        outdir='genplots'
        basedir='./'
        outdir='{Here}{when}{here}'.format(Here=basedir,when=date,here='gen')
        if outdir not in os.listdir(basedir):
            os.system('mkdir -p {od}'.format(od=outdir))

        # f_ex0 = ROOT.TFile('RunIIAutumn18_signal.root', 'READ'); l_ex0 = f_ex0.Get('tree')
        # f_ex1 = ROOT.TFile('RunIIAutumn18_background.root', 'READ'); l_ex1 = f_ex1.Get('tree')
        f_ex0 = ROOT.TFile('RunIIAutumn18_debug.root', 'READ'); l_ex0 = f_ex0.Get('tree')
        f_ex1 = ROOT.TFile('RunIIAutumn18_QCD_debug.root', 'READ'); l_ex1 = f_ex1.Get('tree')

        overlayshapes(outdir,l_ex0)
        overlay2dshapes(outdir,l_ex0)
        overlayshapesSpecial(outdir,l_ex0)
        overlayshapesQCD(outdir,l_ex0,l_ex1)
