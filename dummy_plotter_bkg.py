import ROOT, os, subprocess, sys, optparse, datetime
from array import array
import math
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


def gethist(tree,selcuts,vname,xtitle,nXbins,xmin,xmax, color = 1):
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
        hist2 = ROOT.TH2F("hist2",'',len(Xbins)-1,Xbins,len(Ybins)-1,Ybins)
    elif Xbins!=[]: # Xbins array instead of xmin & xmax
        hist2 = ROOT.TH2F("hist2",'',len(Xbins)-1,Xbins,nYbins,ymin,ymax)
    elif Ybins!=[]: # Ybins array instead of ymin & ymax
        hist2 = ROOT.TH2F("hist2",'',nXbins,xmin,xmax,len(Ybins)-1,Ybins)
    else: # usual
        hist2 = ROOT.TH2F("hist2" ,'',nXbins,xmin,xmax,nYbins,ymin,ymax)
    tree.Draw("{herey}:{herex}>>hist2".format(herex=vnamex,herey=vnamey), selcuts, "goff")
    hist2.GetXaxis().SetTitle(xtitle)
    hist2.GetYaxis().SetTitle(ytitle)
    print(hist2.Integral())
    hist2.SetDirectory(0)
    return hist2


def overlayshapes():
    outdir='genplots'

    basedir='./'
    outdir='{Here}{when}{here}'.format(Here=basedir,when=date,here='gen')
    if outdir not in os.listdir(basedir):
        os.system('mkdir -p {od}'.format(od=outdir))

#
#   change filename!
#
    f_ex0 = ROOT.TFile('RunIIAutumn18_debug.root', 'READ'); l_ex0 = f_ex0.Get('tree')
    # f_ex1 = ROOT.TFile('/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/lhetools/LatinoTreesLHE/WpmWmp_LL_genlevel.root', 'READ'); l_ex1 = f_ex1.Get('tree')
    # f_ex2 = ROOT.TFile('/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/lhetools/LatinoTreesLHE/WpmWmp_TLnLT_genlevel.root', 'READ'); l_ex2 = f_ex2.Get('tree')
    
    sel_str_Gen_matching = "((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1)) || ((Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1) && (Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0))"
    sel_str_Gen_paper = "GenJetAK8_pt[0]>200 && GenJetAK8_pt[1]>200 && abs(GenJetAK8_eta[0])<2.5 && abs(GenJetAK8_eta[1])<2.5 && Gen_M_jj>1250 && GenJetAK8_mass[0]>30 && GenJetAK8_mass[1]>30"
    sel_str_Gen_matching_paper = "(" + sel_str_Gen_paper + ") && (" + sel_str_Gen_matching + ")"
    sel_str_Reco_matching = "((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1)) || ((Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1) && (Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0))"
    sel_str_Reco_paper = "FatJet_pt[0]>200 && FatJet_pt[1]>200 && abs(FatJet_eta[0])<2.5 && abs(FatJet_eta[1])<2.5 && Reco_M_jj>1250 && FatJet_mass[0]>30 && FatJet_mass[1]>30"
    sel_str_Reco_matching_paper = "(" + sel_str_Reco_paper + ") && (" + sel_str_Reco_matching + ")"
    pdraw={
        #### vname : [vname, nXbins, xmin, xmax, sel, _info]
        "Gen_M_jj" : ["Gen_M_jj",326,0,6500,"","no selection"],
        "Gen_M_jj_match" : ["Gen_M_jj",326,0,6500,sel_str_Gen_matching,'matching'],
        "Gen_M_jj_paper" : ["Gen_M_jj",326,0,6500,sel_str_Gen_paper, "VV/VH paper"],
        "Gen_M_jj_match_paper" : ["Gen_M_jj",326,0,6500,sel_str_Gen_matching_paper,"matching & VV/VH paper"],
        # "Gen_M_jY" : ["GenJetAK8_mass",326,0,6500,"((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1)) || ((Gen_matched_Y[1]==1 && Gen_matched_Y[0]==0) && (Gen_matched_YP[1]==0 && Gen_matched_YP[0]==1))",''],
        # "Gen_M_jYP" : ["GenJetAK8_mass",326,0,6500,"((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1)) || ((Gen_matched_YP[1]==1 && Gen_matched_YP[0]==0) && (Gen_matched_Y[1]==0 && Gen_matched_Y[0]==1))",''],
        "Reco_M_jj" : ["Reco_M_jj",326,0,6500,"","no selection"],
        "Reco_M_jj_match" : ["Reco_M_jj",326,0,6500,sel_str_Reco_matching,'matching'],
        "Reco_M_jj_paper" : ["Reco_M_jj",326,0,6500,sel_str_Reco_paper, "VV/VH paper"],
        "Reco_M_jj_match_paper" : ["Reco_M_jj",326,0,6500,sel_str_Reco_matching_paper,"matching & VV/VH paper"],
        # "Reco_M_jY" : ["FatJet_mass",326,0,6500,"((Reco_matched_Y[0]==1 && Reco_matched_Y[1]==0) && (Reco_matched_YP[0]==0 && Reco_matched_YP[1]==1)) || ((Reco_matched_Y[1]==1 && Reco_matched_Y[0]==0) && (Reco_matched_YP[1]==0 && Reco_matched_YP[0]==1))",''],
        # "Reco_M_jYP" : ["FatJet_mass",326,0,6500,"((Reco_matched_YP[0]==1 && Reco_matched_YP[1]==0) && (Reco_matched_Y[0]==0 && Reco_matched_Y[1]==1)) || ((Reco_matched_YP[1]==1 && Reco_matched_YP[0]==0) && (Reco_matched_Y[1]==0 && Reco_matched_Y[0]==1))",''],
        "nFatJet" : ["nFatJet",6,0,5,"",""],
        "FatJet_pt" : ["FatJet_pt",100,0,4000,"",""],
        "FatJet_eta" : ["FatJet_eta",100,-3.5,3.5,"",""],
        "FatJet_phi" : ["FatJet_phi",100,-3.5,3.5,"",""],
        "FatJet_mass" : ["FatJet_mass",100,0,200,"",""],
        "FatJet_msoftdrop" : ["FatJet_msoftdrop",100,0,200,"",""],
        "FatJet_n2b1" : ["FatJet_n2b1",100,-3.5,3.5,"",""],
        "FatJet_n3b1" : ["FatJet_n3b1",100,-3.5,3.5,"",""],
        "FatJet_tau1" : ["FatJet_tau1",100,-3.5,3.5,"",""],
        "FatJet_tau2" : ["FatJet_tau2",100,-3.5,3.5,"",""],
        "FatJet_tau3" : ["FatJet_tau3",100,-3.5,3.5,"",""],
        "FatJet_tau4" : ["FatJet_tau4",100,-3.5,3.5,"",""],
        "nSubJet" : ["nSubJet",6,0,5,"",""],
        "SubJet_mass" : ["SubJet_mass",100,0,200,"",""],
        "SubJet_n2b1" : ["SubJet_n2b1",100,-3.5,3.5,"",""],
        "SubJet_n3b1" : ["SubJet_n3b1",100,-3.5,3.5,"",""],
        "SubJet_tau1" : ["SubJet_tau1",100,-3.5,3.5,"",""],
        "SubJet_tau2" : ["SubJet_tau2",100,-3.5,3.5,"",""],
        "SubJet_tau3" : ["SubJet_tau3",100,-3.5,3.5,"",""],
        "SubJet_tau4" : ["SubJet_tau4",100,-3.5,3.5,"",""],
        "Reco_matched_Y" : ["Reco_matched_Y",100,-3.5,3.5,"",""],
        "Reco_matched_Y" : ["Reco_matched_Y",100,-3.5,3.5,"",""],
        }  
    sel=''
    for ikey,ival in pdraw.items():
        print(ikey)
        h_ex0 = gethist(l_ex0,ival[4],ival[0],ikey,ival[1],ival[2],ival[3],ROOT.kGreen+2)
        # h_ex1 = gethist(l_ex1,sel,ikey,ival[0],ival[1],ival[2],ival[3],'W^{#pm}W^{#mp}_LL',ROOT.kAzure+4)
        # h_ex2 = gethist(l_ex2,sel,ikey,ival[0],ival[1],ival[2],ival[3],'W^{#pm}W^{#mp}_TLnLT',ROOT.kOrange+10)
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
        # h_ex1.Draw('histsame')
        # h_ex2.Draw('histsame')

        leg = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
        leg.SetTextSize(0.045); leg.SetShadowColor(0);        leg.SetFillStyle(0)
        leg.SetTextFont(42);   leg.SetBorderSize(0);
        #leg.SetNColumns(3)
        # leg.AddEntry('NULL', 'W^{#pm}W^{#mp}'  , '')
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
def overlay2dshapes():
    outdir='genplots'

    basedir='./'
    outdir='{Here}{when}{here}'.format(Here=basedir,when=date,here='gen')
    if outdir not in os.listdir(basedir):
        os.system('mkdir -p {od}'.format(od=outdir))

    f_ex0 = ROOT.TFile('RunIIAutumn18_debug.root', 'READ'); l_ex0 = f_ex0.Get('tree')
    # f_ex1 = ROOT.TFile('/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/lhetools/LatinoTreesLHE/WpmWmp_LL_genlevel.root', 'READ'); l_ex1 = f_ex1.Get('tree')
    # f_ex2 = ROOT.TFile('/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/lhetools/LatinoTreesLHE/WpmWmp_TLnLT_genlevel.root', 'READ'); l_ex2 = f_ex2.Get('tree')
    binsForX = array('f',[1000, 1500, 2000, 3000, 4000, 5000, 6000])
    binsForY = array('f',[ 30,  60,  90, 120, 150, 200, 250, 300, 400, 500])
    binsForYP =array('f',[ 30,  60,  90, 120, 150, 200, 250, 300, 400, 500])
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
        "delta_R_qq_Y23_vs_delta_M" : ["delta_R_qq_Y23", "M_X-(M_Y+M_YP)", "delta_R_qq_Y23", "delta M", 200, 200, 0, 0, 2, 6500, "",[],[]],
        "delta_R_qq_YP25_vs_delta_M" : ["delta_R_qq_YP25", "M_X-(M_Y+M_YP)", "delta_R_qq_YP25", "delta M", 200, 200, 0, 0, 2, 6500, "",[],[]],
        "delta_R_qq_Y23_vs_M_Y_ratio" : ["delta_R_qq_Y23", "M_Y/M_X", "delta_R_qq_Y23", "M_Y_ratio", 200, 200, 0, 0, 2, 0.3, "V_pdgId==23",[],[]],
        "delta_R_qq_YP25_vs_M_YP_ratio" : ["delta_R_qq_YP25", "M_YP/M_X", "delta_R_qq_YP25", "M_YP_ratio", 200, 200, 0, 0, 2, 0.3, "V_pdgId==25",[],[]],

        # # "Gen_M_jj_vs_M_jY" : ["Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "M_jY", 50, 50, 0, 0, 6500, 1000, "((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1)) || ((Gen_matched_Y[1]==1 && Gen_matched_Y[0]==0) && (Gen_matched_YP[1]==0 && Gen_matched_YP[0]==1))"],
        # # "Gen_M_jj_vs_M_jYP" : ["Gen_M_jj", "GenJetAK8_mass", "Gen_M_jj", "M_jYP", 50, 50, 0, 0, 6500, 1000, "((Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0) && (Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1)) || ((Gen_matched_YP[1]==1 && Gen_matched_YP[0]==0) && (Gen_matched_Y[1]==0 && Gen_matched_Y[0]==1))"]
        }

    for ikey,ival in pdraw2d.items():
        print(ikey)
        #tree, selcuts, vnamex, vnamey, xtitle, ytitle, nXbins, nYbins, xmin, ymin, xmax, ymax, color = 1
        h_ex0 = get2dhist(l_ex0,ival[10],ival[0],ival[1],ival[2],ival[3],ival[4],ival[5],ival[6],ival[7],ival[8],ival[9],ival[11],ival[12])
        # h_ex1 = gethist(l_ex1,sel,ikey,ival[0],ival[1],ival[2],ival[3],'W^{#pm}W^{#mp}_LL',ROOT.kAzure+4)
        # h_ex2 = gethist(l_ex2,sel,ikey,ival[0],ival[1],ival[2],ival[3],'W^{#pm}W^{#mp}_TLnLT',ROOT.kOrange+10)
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
        h_ex0.GetYaxis().SetTitleOffset(1.2)
        h_ex0.GetXaxis().SetTitleOffset(1.0)
        h_ex0.GetYaxis().SetTitleSize(0.045)
        h_ex0.GetXaxis().SetTitleSize(0.045)


        h_ex0.Draw('colz')
        # h_ex1.Draw('histsame')
        # h_ex2.Draw('histsame')

        leg = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
        leg.SetTextSize(0.045); leg.SetShadowColor(0);        leg.SetFillStyle(0)
        leg.SetTextFont(42);   leg.SetBorderSize(0);
        #leg.SetNColumns(3)
        # leg.AddEntry('NULL', 'W^{#pm}W^{#mp}'  , '')
        leg.AddEntry(h_ex0, 'XToYYPrime'  , 'l')
        leg.AddEntry('NULL', ival[10],'')
        # leg.AddEntry(h_ex1, 'LL'  , 'l')
        # leg.AddEntry(h_ex2, 'TLnLT'  , 'l')
        leg.SetLineColor(ROOT.kWhite)
        leg.SetFillColor(ROOT.kWhite)
        leg.Draw('same')
        
        lat = ROOT.TLatex()
        lat.SetNDC()
        #fout.WriteTObject(canv)
        canv.SaveAs('{od}/{mww}.pdf'.format(od=outdir,mww=ikey))
        canv.SaveAs('{od}/{mww}.png'.format(od=outdir,mww=ikey))


    # special plot
    # "M_X_vs_M_Y" : ["M_X", "M_Y", "M_X", "M_Y", 50, 50, 0, 0, 6500, 500, "",binsForX,binsForY],
    s_hist1 = get2dhist(l_ex0, "(GenJetAK8_pt[0]>200 && abs(GenJetAK8_eta[0])<2.5 && abs(GenJetAK8_eta[1])<2.5 && Gen_M_jj>1250 \
                                       && GenJetAK8_mass[0]>30 && GenJetAK8_mass[1]>30) \
                                      && (((Gen_matched_Y[0]==1 && Gen_matched_Y[1]==0) && (Gen_matched_YP[0]==0 && Gen_matched_YP[1]==1)) \
                                           || ((Gen_matched_Y[0]==0 && Gen_matched_Y[1]==1) && (Gen_matched_YP[0]==1 && Gen_matched_YP[1]==0)))", "M_X", "M_Y", "M_X", "M_Y", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    s_hist2 = get2dhist(l_ex0, "","M_X", "M_Y", "M_X", "M_Y", 0, 0, 0, 0, 0, 0, binsForX, binsForY)
    s_hist1.Divide(s_hist2)

    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)
    canv = ROOT.TCanvas('M_sel_div_M', 'bar', 600, 600)
    canv.cd()
    s_hist1.GetYaxis().SetTitleOffset(1.2)
    s_hist1.GetXaxis().SetTitleOffset(1.0)
    s_hist1.GetYaxis().SetTitleSize(0.045)
    s_hist1.GetXaxis().SetTitleSize(0.045)
    s_hist1.Draw('colz')
    leg = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
    leg.SetTextSize(0.045); leg.SetShadowColor(0);        leg.SetFillStyle(0)
    leg.SetTextFont(42);   leg.SetBorderSize(0);
    leg.AddEntry(s_hist1, 'XToYYPrime'  , 'l')
    leg.SetLineColor(ROOT.kWhite)
    leg.SetFillColor(ROOT.kWhite)
    leg.Draw('same')
    
    lat = ROOT.TLatex()
    lat.SetNDC()
    #fout.WriteTObject(canv)
    canv.SaveAs('{od}/M_sel_div_M.pdf'.format(od=outdir))
    canv.SaveAs('{od}/M_sel_div_M.png'.format(od=outdir))
##############################################


if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('--overlayshapes', action='store_true', dest='overlayshapes'        , default=True, help='blah blah ')
    global opts
    (opts, args) = parser.parse_args()


    if opts.overlayshapes:
        overlayshapes()
        overlay2dshapes()
