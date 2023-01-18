from TauPOG.TauIDSFs.ScoreSFTool import ScoreSFTool
from collections import OrderedDict
from array import array
import argparse
import ROOT
import copy
import json
import os

### Example commands ###

# python test/runTauIDScoreCorrections.py --mc_eff="data/MCEffTauIDScoreCorrections_DM.root" --sf="data/TauID_SF_dm_DeepTau2017v2p1VSjet_2017ReReco.root" --var_name="Decay Mode" --output_folder="score_reweighting_example_DM" --name="TauIDScoreSF_DM.root"
# python test/runTauIDScoreCorrections.py --mc_eff="data/MCEffTauIDScoreCorrections_pT.root" --sf="data/TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco.root" --var_name="p_{T}" --output_folder="score_reweighting_example_pT" --name="TauIDScoreSF_pT.root" --logx

### Parsers ###

parser = argparse.ArgumentParser()
parser.add_argument('--sf',help= 'The name and location of the tau ID scale factor file', default='data/TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco.root')
parser.add_argument('--mc_eff',help= 'The name and location of the tau ID MC efficiency file', default='data/MCEffTauIDScoreCorrections_pT.root')
parser.add_argument('--wp_scores',help= 'The name and location of the json file with the WP to score conversion in', default='data/TauID_WPtoScore_DeepTau2017v2p1VSjet.json')
parser.add_argument('--var_name',help= 'Variable name for plotting', default='p_{T}')
parser.add_argument('--output_folder',help= 'Name of output folder for plots and root file', default='.')
parser.add_argument('--name',help= 'Name of output root file', default='TauIDScore_SF.root')
parser.add_argument("--logx",help="Use logx for plots with variable on x axis",action='store_true')
args = parser.parse_args()

### functions ###

def RemoveSameBins(hist):
  non_zero_diff_bins = []
  for i in range(1, hist.GetNbinsX()+1):
    if not (hist.GetBinContent(i) == 0 or hist.GetBinContent(i) == hist.GetBinContent(i-1)):
      non_zero_diff_bins.append(hist.GetBinLowEdge(i))
  non_zero_diff_bins.append(1000)
  bins = array('f', map(float,non_zero_diff_bins))
  hout = ROOT.TH1D('hout','',len(bins)-1, bins)
  for i in range(0,hout.GetNbinsX()+1):
    hout.SetBinContent(i,hist.GetBinContent(hist.FindBin(hout.GetBinLowEdge(i))))
  return hout

### Make output directory ###

if not os.path.isdir(args.output_folder): os.system("mkdir {}".format(args.output_folder))

### Load in scale factors ###

sf_h = OrderedDict()
sf_file = ROOT.TFile(args.sf)
for ind, key in enumerate(sf_file.GetListOfKeys()):
  name = key.GetName()
  rt = sf_file.Get(name)

  if type(rt) == type(ROOT.TH1F()):
    sf_h[name] = rt.Clone()

  elif type(rt) == type(ROOT.TF1()):
    if "_cent" not in name: continue
    name = name.replace("_cent","")
    # convert to histogram
    h = rt.DoCreateHistogram(0,500) # need to find a better way to do this but for now it is fine
    h = RemoveSameBins(h)
    h.SetName(name)
    sf_h[name] = copy.deepcopy(h)

    up = sf_file.Get(name+"_up")
    up_h = up.DoCreateHistogram(0,500)
    up_h = RemoveSameBins(up_h)

    for i in range(0,sf_h[name].GetNbinsX()+2): 
      sf_h[name].SetBinError(i,up_h.GetBinContent(i)-sf_h[name].GetBinContent(i))
    del up_h

### Load in mc efficiencies ###

mc_eff_h = OrderedDict()
mc_eff_file = ROOT.TFile(args.mc_eff)
for key in mc_eff_file.GetListOfKeys(): mc_eff_h[key.GetName()] = copy.deepcopy(mc_eff_file.Get(key.GetName()))

### Load into class ###

ssf = ScoreSFTool()

ssf.InputHistograms(sf_h,"sf")
ssf.PlotSFs(logx=args.logx,x_label=args.var_name,with_spline=False,folder=args.output_folder)
ssf.InputHistograms(mc_eff_h,"mc")

### Calculate SFs with respect to the loosest WP ##
loosest_sf = copy.deepcopy(ssf.GetHistogram("sf",ssf.sorted_wp[0]))
ssf.ScaleHistogramsByHistogram("sf",copy.deepcopy(loosest_sf),divide=True)

### Get data efficiencies ###
ssf.CalculateHistograms("data")
ssf.PlotEfficienciesAndSFs(logx=args.logx,x_label=args.var_name,folder=args.output_folder)

### Convert to WP binned histograms ###
ssf.ConvertToWPBinnedHistograms(rebin_threshold=0.1)

# Change bin names and plot
for k, v in ssf.rebinned_bins.iteritems():
  replace = range(0,len(v))
  replace = [x+0.5 for x in replace]
  replace_labels = []
  for ind, i in enumerate(v):
    if ind+1 != len(v):
      replace_labels.append(i+"&!"+v[ind+1])
    else:
      replace_labels.append(i)

  ssf.PlotEfficienciesAndSFs(logx=False,ratio_range=[0,2],extra_ratio_line=[0.5,1.5],replace=replace,replace_labels=replace_labels,individual=k,folder=args.output_folder)

### Scale back to full SF and move to score bins ###

ssf.ScaleWPHistogramsByHistogram("sf",loosest_sf)
with open(args.wp_scores) as json_file: score_dict = json.load(json_file, object_pairs_hook=OrderedDict)
ssf.ConvertToScoreBinnedHistograms(score_dict)
#ssf.FitSpline("sf")
ssf.PlotSFs(x_label="Score",with_spline=False,folder=args.output_folder)

### Get TF2s from histograms ###

cent,up,down = ssf.DumpSFTF2(variation="cent"),ssf.DumpSFTF2(variation="up"),ssf.DumpSFTF2(variation="down")
print "SF formula:"
cent.Print("all")
print "Up variation:"
up.Print("all")
print "Down variation:"
down.Print("all")

### Write to file ###
fout = ROOT.TFile(args.output_folder+"/"+args.name, 'RECREATE')
cent.Write()
up.Write()
down.Write()
fout.Close()
print "Created {}/{}".format(args.output_folder,args.name)
