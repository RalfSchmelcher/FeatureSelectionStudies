import os,sys, re
import ROOT
from ROOT import TLorentzVector
from array import array
ROOT.gROOT.SetBatch()


max_events=-1
first_event=0
last_event=-1
if len(sys.argv)>2:
    input_file=sys.argv[1]
    job_id=sys.argv[2]
    fIn=ROOT.TFile.Open(input_file)
    fName="RunIIAutumn18_{}".format(job_id)
    if len(sys.argv)>4:
        first_event=int(sys.argv[3])
        last_event=int(sys.argv[4])
else:
    fIn=ROOT.TFile.Open("/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18/F6B34C92-581B-034F-8283-65DC1F68BD7F.root") 
    job_id = "debug"
    fName="RunIIAutumn18_QCD_debug"
    max_events=500

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
outTree.Branch("Reco_matched_Y",Reco_matched_Y,"Reco_matched_Y[nSubJet]/I")
outTree.Branch("Reco_matched_YP",Reco_matched_YP,"Reco_matched_YP[nSubJet]/I")
outTree.Branch("Reco_M_jj",Reco_M_jj,"Reco_M_jj/F")


nevts=0
if last_event<0:
    last_event = ttree.GetEntries()

for ievent in range(first_event, last_event):
    ttree.GetEntry(ievent)
    ev = ttree
    if nevts > max_events and max_events>0: break
    nevts+=1
    if ev.nFatJet < 2: continue
    h_events.Fill(1)
    if max_events>0  and nevts%100==0: print("processing event",nevts)
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
        FatJet_pt[j] = ev.FatJet_pt[j]
        FatJet_eta[j] = ev.FatJet_eta[j]
        FatJet_phi[j] = ev.FatJet_phi[j]
        FatJet_mass[j] = ev.FatJet_mass[j]

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


