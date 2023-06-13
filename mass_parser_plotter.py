import ROOT
import glob
from hvt_mass_hypotheses_util import mass_points
import os, sys
import json
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)


mass_dict_X = {'total': 0}
mass_dict_Y = {'total': 0}
mass_dict_YP = {'total': 0}
mass_dict_branch = {}


jsonFiles_X = glob.glob("/afs/desy.de/user/s/schmelch/GenLevelStudies/customTrees/workdir/mass_parser_jsons/mass_dict_X_*.json")
for file_index, file_name in enumerate(jsonFiles_X):
    temp_X = json.load(open(file_name))
    for i, (key, value) in enumerate(temp_X.items()):
        # print(file_index, mass_dict_X['total'])
        if key in mass_dict_X:
            mass_dict_X[key] += temp_X[key]
        else:
            mass_dict_X[key] = temp_X[key]

hist_X = ROOT.TH1F("hist", "M_X", len(mass_dict_X)-1, 0, len(mass_dict_X))

for i, (key, value) in enumerate(mass_dict_X.items()):
    if key=='total': continue
    hist_X.Fill(i, value)
    hist_X.GetXaxis().SetBinLabel(i, str(key).zfill(4))
hist_X.LabelsOption('a')


canv_X = ROOT.TCanvas('M_X', 'bar', 600, 600)
canv_X.cd()
hist_X.Draw('hist')
leg_X = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
leg_X.SetTextSize(0.045); leg_X.SetShadowColor(0);        leg_X.SetFillStyle(0)
leg_X.SetTextFont(42);   leg_X.SetBorderSize(0);
leg_X.AddEntry('NULL', str(mass_dict_X['total'])+" total",'')
leg_X.SetLineColor(ROOT.kWhite)
leg_X.SetFillColor(ROOT.kWhite)
leg_X.Draw('same')
# canv_X.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_X_count.pdf')
canv_X.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_X_count.png')




jsonFiles_Y = glob.glob("/afs/desy.de/user/s/schmelch/GenLevelStudies/customTrees/workdir/mass_parser_jsons/mass_dict_Y_*.json")
for file_index, file_name in enumerate(jsonFiles_Y):
    temp_Y = json.load(open(file_name))
    for i, (key, value) in enumerate(temp_Y.items()):
        # print(file_index, mass_dict_Y['total'])
        if key in mass_dict_Y:
            mass_dict_Y[key] += temp_Y[key]
        else:
            mass_dict_Y[key] = temp_Y[key]

hist_Y = ROOT.TH1F("hist", "M_Y", len(mass_dict_Y)-1, 0, len(mass_dict_Y))

for i, (key, value) in enumerate(mass_dict_Y.items()):
    if key=='total': continue
    hist_Y.Fill(i, value)
    hist_Y.GetXaxis().SetBinLabel(i, str(key).zfill(3))
hist_Y.LabelsOption('a')


canv_Y = ROOT.TCanvas('M_Y', 'bar', 600, 600)
canv_Y.cd()
hist_Y.Draw('hist')
leg_Y = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
leg_Y.SetTextSize(0.045); leg_Y.SetShadowColor(0);        leg_Y.SetFillStyle(0)
leg_Y.SetTextFont(42);   leg_Y.SetBorderSize(0);
leg_Y.AddEntry('NULL', str(mass_dict_Y['total'])+" total",'')
leg_Y.SetLineColor(ROOT.kWhite)
leg_Y.SetFillColor(ROOT.kWhite)
leg_Y.Draw('same')
# canv_Y.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_Y_count.pdf')
canv_Y.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_Y_count.png')



jsonFiles_YP = glob.glob("/afs/desy.de/user/s/schmelch/GenLevelStudies/customTrees/workdir/mass_parser_jsons/mass_dict_YP_*.json")
for file_index, file_name in enumerate(jsonFiles_YP):
    temp_YP = json.load(open(file_name))
    for i, (key, value) in enumerate(temp_YP.items()):
        # print(file_index, mass_dict_YP['total'])
        if key in mass_dict_YP:
            mass_dict_YP[key] += temp_YP[key]
        else:
            mass_dict_YP[key] = temp_YP[key]

hist_YP = ROOT.TH1F("hist", "M_YP", len(mass_dict_YP)-1, 0, len(mass_dict_YP))

for i, (key, value) in enumerate(mass_dict_YP.items()):
    if key=='total': continue
    hist_YP.Fill(i, value)
    hist_YP.GetXaxis().SetBinLabel(i, str(key).zfill(3))
hist_YP.LabelsOption('a')


canv_YP = ROOT.TCanvas('M_YP', 'bar', 600, 600)
canv_YP.cd()
hist_YP.Draw('hist')
leg_YP = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
leg_YP.SetTextSize(0.045); leg_YP.SetShadowColor(0);        leg_YP.SetFillStyle(0)
leg_YP.SetTextFont(42);   leg_YP.SetBorderSize(0);
leg_YP.AddEntry('NULL', str(mass_dict_YP['total'])+" total",'')
leg_YP.SetLineColor(ROOT.kWhite)
leg_YP.SetFillColor(ROOT.kWhite)
leg_YP.Draw('same')
# canv_YP.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_YP_count.pdf')
canv_YP.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_YP_count.png')



jsonFiles_branch = glob.glob("/afs/desy.de/user/s/schmelch/GenLevelStudies/customTrees/workdir/mass_parser_jsons/mass_dict_branch_*.json")
for file_index, file_name in enumerate(jsonFiles_branch):
    temp_branch = json.load(open(file_name))
    for i, (key, value) in enumerate(temp_branch.items()):
        # print(file_index, mass_dict_branch['total'])
        if key in mass_dict_branch:
            mass_dict_branch[key] += temp_branch[key]
        else:
            mass_dict_branch[key] = temp_branch[key]
with open("/afs/desy.de/user/s/schmelch/GenLevelStudies/mass_parser_jsons/mass_dict_branch_total.json", "w") as outfile:
    json.dump(mass_dict_branch, outfile)

# hist_branch = ROOT.TH1F("hist", "M_branch", len(mass_dict_branch)-1, 0, len(mass_dict_branch))

# for i, (key, value) in enumerate(mass_dict_branch.items()):
#     hist_branch.Fill(i, value)
#     hist_branch.GetXaxis().SetBinLabel(i, key)
# # hist_branch.LabelsOption('a')


# canv_branch = ROOT.TCanvas('M_branch', 'bar', 600, 600)
# canv_branch.cd()
# hist_branch.Draw('hist')
# leg_branch = ROOT.TLegend(0.5, 0.65, 0.8, 0.875)
# leg_branch.SetTextSize(0.045); leg_branch.SetShadowColor(0);        leg_branch.SetFillStyle(0)
# leg_branch.SetTextFont(42);   leg_branch.SetBorderSize(0);
# leg_branch.SetLineColor(ROOT.kWhite)
# leg_branch.SetFillColor(ROOT.kWhite)
# leg_branch.Draw('same')
# # canv_branch.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_branch_count.pdf')
# canv_branch.SaveAs('/afs/desy.de/user/s/schmelch/GenLevelStudies/M_branch_count.png')
