import os,sys, re 
from hvt_mass_hypotheses_util import mass_points
import json
import ROOT
from ROOT import TLorentzVector
from array import array
ROOT.gROOT.SetBatch()


try:
    os.mkdir('mass_parser_jsons')
    # os.mkdir('/afs/desy.de/user/s/schmelch/GenLevelStudies/mass_parser_jsons')
except OSError as error:
    pass

max_events=-1
if len(sys.argv)>2:
    input_file=sys.argv[1]
    job_id=sys.argv[2]
    fIn=ROOT.TFile.Open(input_file)
    fName="RunIIAutumn18_{}".format(job_id)
else:
    fIn=ROOT.TFile.Open("/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/XToYYprime_YToQQ_YprimeToQQ_narrow_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIIAutumn18/B2BEAC83-CEA3-A241-B0D3-2E1461C38E09.root") 
    job_id = "debug"
    fName="RunIIAutumn18_debug"
    max_events=8000


mass_dict_X = {'total': 0}
mass_dict_Y = {'total': 0}
mass_dict_YP = {'total': 0}
mass_dict_branch = {}

ttree=fIn.Get('Events')
fOut = ROOT.TFile("%s.root"%fName,"recreate")
h_events = ROOT.TH1F('h_events','h_events',2,0,2)
outTree = ROOT.TTree('tree',"Events")
maxn=4
nV      = array("i",[-99])
V_pt    = array("f",[-9999., -9999., -9999., -9999.])
V_eta   = array("f",maxn*[-9999.])
V_phi   = array("f",maxn*[-9999.])
V_mass   = array("f",maxn*[-9999.])
V_pdgId = array("i",maxn*[-99])
M_X = array("f", [-9999.])
M_Y = array("f", [-9999.])
M_YP = array("f", [-9999.])
delta_R_qq_Y23 = array("f", [-9999.])
delta_R_qq_YP25 = array("f", [-9999.])
nGenJetAK8 = array("i",[-99])
GenJetAK8_pt = array("f", 2*[-9999.])
GenJetAK8_eta = array("f", 2*[-9999.])
GenJetAK8_phi = array("f", 2*[-9999.])
GenJetAK8_mass = array("f", 2*[-9999.])
Gen_matched_Y = array("i", [-1,-1])
Gen_matched_YP = array("i", [-1,-1])
Gen_M_jj = array("f", [-9999.])
nFatJet = array("i",[-99])
FatJet_pt = array("f", 2*[-9999.])
FatJet_eta = array("f", 2*[-9999.])
FatJet_phi = array("f", 2*[-9999.])
FatJet_mass = array("f", 2*[-9999.])
FatJet_msoftdrop = array("f", 2*[-9999.])
FatJet_n2b1 = array("f", 2*[-9999.])
FatJet_n3b1 = array("f", 2*[-9999.])
FatJet_tau1 = array("f", 2*[-9999.])
FatJet_tau2 = array("f", 2*[-9999.])
FatJet_tau3 = array("f", 2*[-9999.])
FatJet_tau4 = array("f", 2*[-9999.])
FatJet_subJetIdx1 = array("f", 2*[-9999.])
FatJet_subJetIdx2 = array("f", 2*[-9999.])
nSubJet = array("i",[-99])
SubJet_mass = array("f", 4*[-9999.])
SubJet_n2b1 = array("f", 4*[-9999.])
SubJet_n3b1 = array("f", 4*[-9999.])
SubJet_tau1 = array("f", 4*[-9999.])
SubJet_tau2 = array("f", 4*[-9999.])
SubJet_tau3 = array("f", 4*[-9999.])
SubJet_tau4 = array("f", 4*[-9999.])
Reco_matched_Y = array("i", [-1,-1])
Reco_matched_YP = array("i", [-1,-1])
Reco_M_jj = array("f", [-9999.])
GenJet0_deltaR_Y = array("f", [-9999.])
GenJet1_deltaR_Y = array("f", [-9999.])
GenJet0_deltaR_YP = array("f", [-9999.])
GenJet1_deltaR_YP = array("f", [-9999.])
FatJet0_deltaR_Y = array("f", [-9999.])
FatJet1_deltaR_Y = array("f", [-9999.])
FatJet0_deltaR_YP = array("f", [-9999.])
FatJet1_deltaR_YP = array("f", [-9999.])


outTree.Branch("nV",nV,"nV/I")
outTree.Branch("V_pt",V_pt,"V_pt[nV]/F")
outTree.Branch("V_eta",V_eta,"V_eta[nV]/F")
outTree.Branch("V_phi",V_phi,"V_phi[nV]/F")
outTree.Branch("V_mass",V_mass,"V_mass[nV]/F")
outTree.Branch("V_pdgId",V_pdgId,"LepGood_pdgId[nV]/I")
outTree.Branch("M_X",M_X,"M_X/F")
outTree.Branch("M_Y",M_Y,"M_Y/F")
outTree.Branch("M_YP",M_YP,"M_YP/F")
outTree.Branch("delta_R_qq_Y23",delta_R_qq_Y23,"delta_R_qq_Y23/F")
outTree.Branch("delta_R_qq_YP25",delta_R_qq_YP25,"delta_R_qq_YP25/F")
outTree.Branch("nGenJetAK8",nGenJetAK8,"nGenJetAK8/I")
outTree.Branch("GenJetAK8_pt",GenJetAK8_pt,"GenJetAK8_pt[nGenJetAK8]/F")
outTree.Branch("GenJetAK8_eta",GenJetAK8_eta,"GenJetAK8_eta[nGenJetAK8]/F")
outTree.Branch("GenJetAK8_phi",GenJetAK8_phi,"GenJetAK8_phi[nGenJetAK8]/F")
outTree.Branch("GenJetAK8_mass",GenJetAK8_mass,"GenJetAK8_mass[nGenJetAK8]/F")
outTree.Branch("Gen_matched_Y",Gen_matched_Y,"Gen_matched_Y[nGenJetAK8]/I")
outTree.Branch("Gen_matched_YP",Gen_matched_YP,"Gen_matched_YP[nGenJetAK8]/I")
outTree.Branch("Gen_M_jj",Gen_M_jj,"Gen_M_jj/F")
outTree.Branch("nFatJet",nFatJet,"nFatJet/I")
outTree.Branch("FatJet_pt",FatJet_pt,"FatJet_pt[nFatJet]/F")
outTree.Branch("FatJet_eta",FatJet_eta,"FatJet_eta[nFatJet]/F")
outTree.Branch("FatJet_phi",FatJet_phi,"FatJet_phi[nFatJet]/F")
outTree.Branch("FatJet_mass",FatJet_mass,"FatJet_mass[nFatJet]/F")
outTree.Branch("FatJet_msoftdrop",FatJet_msoftdrop,"FatJet_msoftdrop[nFatJet]/F")
outTree.Branch("FatJet_n2b1",FatJet_n2b1,"FatJet_n2b1[nFatJet]/F")
outTree.Branch("FatJet_n3b1",FatJet_n3b1,"FatJet_n3b1[nFatJet]/F")
outTree.Branch("FatJet_tau1",FatJet_tau1,"FatJet_tau1[nFatJet]/F")
outTree.Branch("FatJet_tau2",FatJet_tau2,"FatJet_tau2[nFatJet]/F")
outTree.Branch("FatJet_tau3",FatJet_tau3,"FatJet_tau3[nFatJet]/F")
outTree.Branch("FatJet_tau4",FatJet_tau4,"FatJet_tau4[nFatJet]/F")
outTree.Branch("FatJet_subJetIdx1",FatJet_subJetIdx1,"FatJet_subJetIdx1[nFatJet]/F")
outTree.Branch("FatJet_subJetIdx2",FatJet_subJetIdx2,"FatJet_subJetIdx2[nFatJet]/F")
outTree.Branch("nSubJet",nSubJet,"nSubJet/I")
outTree.Branch("SubJet_mass",SubJet_mass,"SubJet_mass[nSubJet]/F")
outTree.Branch("SubJet_n2b1",SubJet_n2b1,"SubJet_n2b1[nSubJet]/F")
outTree.Branch("SubJet_n3b1",SubJet_n3b1,"SubJet_n3b1[nSubJet]/F")
outTree.Branch("SubJet_tau1",SubJet_tau1,"SubJet_tau1[nSubJet]/F")
outTree.Branch("SubJet_tau2",SubJet_tau2,"SubJet_tau2[nSubJet]/F")
outTree.Branch("SubJet_tau3",SubJet_tau3,"SubJet_tau3[nSubJet]/F")
outTree.Branch("SubJet_tau4",SubJet_tau4,"SubJet_tau4[nSubJet]/F")
outTree.Branch("Reco_matched_Y",Reco_matched_Y,"Reco_matched_Y[2]/I")
outTree.Branch("Reco_matched_YP",Reco_matched_YP,"Reco_matched_YP[2]/I")
outTree.Branch("Reco_M_jj",Reco_M_jj,"Reco_M_jj/F")
outTree.Branch("GenJet0_deltaR_Y",GenJet0_deltaR_Y,"GenJet0_deltaR_Y/F")
outTree.Branch("GenJet1_deltaR_Y",GenJet1_deltaR_Y,"GenJet1_deltaR_Y/F")
outTree.Branch("GenJet0_deltaR_YP",GenJet0_deltaR_YP,"GenJet0_deltaR_YP/F")
outTree.Branch("GenJet1_deltaR_YP",GenJet1_deltaR_YP,"GenJet1_deltaR_YP/F")
outTree.Branch("FatJet0_deltaR_Y",FatJet0_deltaR_Y,"FatJet0_deltaR_Y/F")
outTree.Branch("FatJet1_deltaR_Y",FatJet1_deltaR_Y,"FatJet1_deltaR_Y/F")
outTree.Branch("FatJet0_deltaR_YP",FatJet0_deltaR_YP,"FatJet0_deltaR_YP/F")
outTree.Branch("FatJet1_deltaR_YP",FatJet1_deltaR_YP,"FatJet1_deltaR_YP/F")


nevts=0
for ev in ttree:
    if nevts > max_events and max_events>0: break
    leptons=[];quarks=[]; bosons=[];
    if ev.nFatJet < 2: continue

    for i in range(ev.nGenPart):
        if abs(ev.GenPart_pdgId[i]) in [23,25] and ev.GenPart_status[i] == 22:
            mpid=ev.GenPart_pdgId[ev.GenPart_genPartIdxMother[i]] if ev.GenPart_genPartIdxMother[i] != -1 else 0
            info={'pdgId':ev.GenPart_pdgId[i],'pt':ev.GenPart_pt[i],'eta':ev.GenPart_eta[i],'phi':ev.GenPart_phi[i],'mass':ev.GenPart_mass[i],'mpdgId':mpid}
            bosons.append(info)
        if abs(ev.GenPart_pdgId[i]) < 6:
            mpid=ev.GenPart_pdgId[ev.GenPart_genPartIdxMother[i]] if ev.GenPart_genPartIdxMother[i] != -1 else 0
            info={'pdgId':ev.GenPart_pdgId[i],'pt':ev.GenPart_pt[i],'eta':ev.GenPart_eta[i],'phi':ev.GenPart_phi[i],'mass':ev.GenPart_mass[i],'mpdgId':mpid}
            quarks.append(info)
        else : continue

    if len(bosons) < 2: continue
    h_events.Fill(1)
    nevts+=1
    nV[0]=len(bosons)
    for j in range(nV[0]):
        V_pt[j]             =bosons[j]['pt'];              
        V_eta[j]            =bosons[j]['eta'];   
        V_pdgId[j]          =bosons[j]['pdgId'];             
        V_phi[j]            =bosons[j]['phi'];
        V_mass[j]            =bosons[j]['mass'];
    for mX,mY,mYP in mass_points():
        mass_hyp_branch = "GenModel_XToYYprime_MX%i_MY%i_MYprime%i"%(mX, mY, mYP)
        if not hasattr(ev, mass_hyp_branch):continue
        if getattr(ev,mass_hyp_branch) == 1:
            M_X[0] = mX
            M_Y[0] = mY
            M_YP[0] = mYP
            if mX in mass_dict_X:
                mass_dict_X[mX] += 1
            else:
                mass_dict_X[mX] = 1
            mass_dict_X['total'] += 1

            if mY in mass_dict_Y:
                mass_dict_Y[mY] += 1
            else:
                mass_dict_Y[mY] = 1
            mass_dict_Y['total'] += 1

            if mYP in mass_dict_YP:
                mass_dict_YP[mYP] += 1
            else:
                mass_dict_YP[mYP] = 1
            mass_dict_YP['total'] += 1
            if mass_hyp_branch in mass_dict_branch:
                mass_dict_branch[mass_hyp_branch] += 1
            else:
                mass_dict_branch[mass_hyp_branch] = 1
            break

    nq = len(quarks)
    Y_quarks_vectors = []
    for iq in range(nq):
        if quarks[iq]["mpdgId"] == 23:
            fourvector = TLorentzVector()
            fourvector.SetPtEtaPhiM(quarks[iq]["pt"], quarks[iq]["eta"], quarks[iq]["phi"], quarks[iq]["mass"])
            Y_quarks_vectors.append(fourvector)
    delta_R_qq_Y23[0] = Y_quarks_vectors[0].DeltaR(Y_quarks_vectors[1])
    # q = __import__("functools").partial(__import__("os")._exit, 0)  # FIXME
    # __import__("IPython").embed()  # FIXME
    YP_quarks_vectors = []
    for iq in range(nq):
        if quarks[iq]["mpdgId"] == 25:
            fourvector = TLorentzVector()
            fourvector.SetPtEtaPhiM(quarks[iq]["pt"], quarks[iq]["eta"], quarks[iq]["phi"], quarks[iq]["mass"])
            YP_quarks_vectors.append(fourvector)
    delta_R_qq_YP25[0] = YP_quarks_vectors[0].DeltaR(YP_quarks_vectors[1])

    jets_vectors = []
    for j in range(2):
        jet_fourvector = TLorentzVector()
        jet_fourvector.SetPtEtaPhiM(ev.GenJetAK8_pt[j], ev.GenJetAK8_eta[j], ev.GenJetAK8_phi[j], ev.GenJetAK8_mass[j])
        jets_vectors.append(jet_fourvector)

    nGenJetAK8[0] = len(jets_vectors)
    for j in range(2):
        GenJetAK8_pt[j] = ev.GenJetAK8_pt[j]
        GenJetAK8_eta[j] = ev.GenJetAK8_eta[j]
        GenJetAK8_phi[j] = ev.GenJetAK8_phi[j]
        GenJetAK8_mass[j] = ev.GenJetAK8_mass[j]

    if len(jets_vectors) < 2: continue
    Gen_matched_Y[0] = (jets_vectors[0].DeltaR(Y_quarks_vectors[0]) < 0.6) and (jets_vectors[0].DeltaR(Y_quarks_vectors[1]) < 0.6)
    Gen_matched_Y[1] = (jets_vectors[1].DeltaR(Y_quarks_vectors[0]) < 0.6) and (jets_vectors[1].DeltaR(Y_quarks_vectors[1]) < 0.6)
    Gen_matched_YP[0] = (jets_vectors[0].DeltaR(YP_quarks_vectors[0]) < 0.6) and (jets_vectors[0].DeltaR(YP_quarks_vectors[1]) < 0.6)
    Gen_matched_YP[1] = (jets_vectors[1].DeltaR(YP_quarks_vectors[0]) < 0.6) and (jets_vectors[1].DeltaR(YP_quarks_vectors[1]) < 0.6)
    GenJet0_deltaR_Y = [jets_vectors[0].DeltaR(Y_quarks_vectors[0]), jets_vectors[0].DeltaR(Y_quarks_vectors[1])]
    GenJet1_deltaR_Y = [jets_vectors[1].DeltaR(Y_quarks_vectors[0]), jets_vectors[1].DeltaR(Y_quarks_vectors[1])]
    GenJet0_deltaR_YP = [jets_vectors[0].DeltaR(YP_quarks_vectors[0]), jets_vectors[0].DeltaR(YP_quarks_vectors[1])]
    GenJet1_deltaR_YP = [jets_vectors[1].DeltaR(YP_quarks_vectors[0]), jets_vectors[1].DeltaR(YP_quarks_vectors[1])]


    Gen_M_jj[0] = (jets_vectors[0] + jets_vectors[1]).M()


    fatJet_vectors = []
    for j in range(2):
        jet_fourvector = TLorentzVector()
        jet_fourvector.SetPtEtaPhiM(ev.FatJet_pt[j], ev.FatJet_eta[j], ev.FatJet_phi[j], ev.FatJet_mass[j])
        fatJet_vectors.append(jet_fourvector)

    nFatJet[0] = len(fatJet_vectors)
    subjet_idx = []
    for j in range(2):
        subjet_idx.append(ev.FatJet_subJetIdx1[j])
        subjet_idx.append(ev.FatJet_subJetIdx2[j])
        FatJet_subJetIdx1[j] = ev.FatJet_subJetIdx1[j]
        FatJet_subJetIdx2[j] = ev.FatJet_subJetIdx2[j]
        FatJet_msoftdrop[j] = ev.FatJet_msoftdrop[j]
        FatJet_n2b1[j] = ev.FatJet_n2b1[j]
        FatJet_n3b1[j] = ev.FatJet_n3b1[j]
        FatJet_tau1[j] = ev.FatJet_tau1[j]
        FatJet_tau2[j] = ev.FatJet_tau2[j]
        FatJet_tau3[j] = ev.FatJet_tau3[j]
        FatJet_tau4[j] = ev.FatJet_tau4[j]
        FatJet_pt[j] = ev.FatJet_pt[j]
        FatJet_eta[j] = ev.FatJet_eta[j]
        FatJet_phi[j] = ev.FatJet_phi[j]
        FatJet_mass[j] = ev.FatJet_mass[j]

    if len(fatJet_vectors) < 2: continue
    Reco_matched_Y[0] = (fatJet_vectors[0].DeltaR(Y_quarks_vectors[0]) < 0.6) and (fatJet_vectors[0].DeltaR(Y_quarks_vectors[1]) < 0.6)
    Reco_matched_Y[1] = (fatJet_vectors[1].DeltaR(Y_quarks_vectors[0]) < 0.6) and (fatJet_vectors[1].DeltaR(Y_quarks_vectors[1]) < 0.6)
    Reco_matched_YP[0] = (fatJet_vectors[0].DeltaR(YP_quarks_vectors[0]) < 0.6) and (fatJet_vectors[0].DeltaR(YP_quarks_vectors[1]) < 0.6)
    Reco_matched_YP[1] = (fatJet_vectors[1].DeltaR(YP_quarks_vectors[0]) < 0.6) and (fatJet_vectors[1].DeltaR(YP_quarks_vectors[1]) < 0.6)
    FatJet0_deltaR_Y = [fatJet_vectors[0].DeltaR(Y_quarks_vectors[0]), fatJet_vectors[0].DeltaR(Y_quarks_vectors[1])]
    FatJet1_deltaR_Y = [fatJet_vectors[1].DeltaR(Y_quarks_vectors[0]), fatJet_vectors[1].DeltaR(Y_quarks_vectors[1])]
    FatJet0_deltaR_YP = [fatJet_vectors[0].DeltaR(YP_quarks_vectors[0]), fatJet_vectors[0].DeltaR(YP_quarks_vectors[1])]
    FatJet1_deltaR_YP = [fatJet_vectors[1].DeltaR(YP_quarks_vectors[0]), fatJet_vectors[1].DeltaR(YP_quarks_vectors[1])]

    Reco_M_jj[0] = (fatJet_vectors[0] + fatJet_vectors[1]).M()

    
    if ev.nSubJet < 4: continue
    nSubJet[0] = ev.nSubJet
    for i,j in enumerate(subjet_idx):
        SubJet_mass[i] = ev.SubJet_mass[j]
        SubJet_n2b1[i] = ev.SubJet_n2b1[j]
        SubJet_n3b1[i] = ev.SubJet_n3b1[j]
        SubJet_tau1[i] = ev.SubJet_tau1[j]
        SubJet_tau2[i] = ev.SubJet_tau2[j]
        SubJet_tau3[i] = ev.SubJet_tau3[j]
        SubJet_tau4[i] = ev.SubJet_tau4[j]
        

    outTree.Fill()
#h_events.Write();
#h_mTT.Write();
fOut.Write();    fOut.Close()

with open("mass_parser_jsons/mass_dict_X_{}.json".format(job_id), "w") as outfile:
    json.dump(mass_dict_X, outfile)
with open("mass_parser_jsons/mass_dict_Y_{}.json".format(job_id), "w") as outfile:
    json.dump(mass_dict_Y, outfile)
with open("mass_parser_jsons/mass_dict_YP_{}.json".format(job_id), "w") as outfile:
    json.dump(mass_dict_YP, outfile)
with open("mass_parser_jsons/mass_dict_branch_{}.json".format(job_id), "w") as outfile:
    json.dump(mass_dict_branch, outfile)
#h_mTT.Draw("ehist");

#c.SaveAs("/eos/user/a/anmehta/www/VVsemilep/mTT.pdf");
#c.SaveAs("/eos/user/a/anmehta/www/VVsemilep/mTT.png");
