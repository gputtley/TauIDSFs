from TauPOG.TauIDSFs.ScoreSFTool import ScoreSFTool
from collections import OrderedDict
from array import array
import argparse
import ROOT
import copy
import json

parser = argparse.ArgumentParser()
parser.add_argument('--sf',help= 'The name and location of the tau ID scale factor file', default='data/TauID_SF_pt_DeepTau2017v2p1VSjet_2018ReReco.root')
parser.add_argument('--mc_eff',help= 'The name and location of the tau ID MC efficiency file', default='/vols/cms/gu18/CrabCMSSW/CMSSW_10_2_19/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/mc_efficiencies.root')
parser.add_argument('--wp_scores',help= 'The name and location of the json file with the WP to score conversion in', default='data/TauID_WPtoScore_DeepTau2017v2p1VSjet.json')
parser.add_argument('--var_name',help= 'Variable name for plotting', default='p_{T}')
parser.add_argument('--shift',help= 'Shift to use can be up or down or unset for nominal', default='cent')
parser.add_argument('--output_folder',help= 'Name of output folder for plots and root file', default='.')
parser.add_argument('--name',help= 'Name of output root file', default='TauIDScore_SF.root')
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

### Load in scale factors ###

sf_h = OrderedDict()
sf_file = ROOT.TFile(args.sf)
for ind, key in enumerate(sf_file.GetListOfKeys()):
  name = key.GetName()
  rt = sf_file.Get(name)
  if type(rt) == type(ROOT.TH1F()):
    sf_h[name] = rt.Clone()
    for b in range(0,sf_h[name].GetNbinsX()+2):
      if args.shift == "up":
        sf_h[name].SetBinContent(i,sf_h[name].GetBinContent(i)+sf_h[name].GetBinError(i))
      elif args.shift == "down":
        sf_h[name].SetBinContent(i,sf_h[name].GetBinContent(i)-sf_h[name].GetBinError(i))
  elif type(rt) == type(ROOT.TF1()):
    if "_"+args.shift not in name: continue
    name = name.replace("_"+args.shift,"")
    # convert to histogram
    h = rt.DoCreateHistogram(0,500) # need to find a better way to do this but for now it is fine
    h = RemoveSameBins(h)
    h.SetName(name)
    sf_h[name] = copy.deepcopy(h)

# set uncertainties to zero
for k, v in sf_h.iteritems():
  for i in range(0,sf_h[k].GetNbinsX()+2):
    sf_h[k].SetBinError(i,0)

### Load in mc efficiencies ###

mc_eff_h = OrderedDict()
mc_eff_file = ROOT.TFile(args.mc_eff)
for key in mc_eff_file.GetListOfKeys(): mc_eff_h[key.GetName()] = copy.deepcopy(mc_eff_file.Get(key.GetName()))

### Load into class ###

ssf = ScoreSFTool()

ssf.InputHistograms(sf_h,"sf")
ssf.InputHistograms(mc_eff_h,"mc")

#ssf.PrintHistograms()

loosest_sf = copy.deepcopy(ssf.GetHistogram("sf",ssf.sorted_wp[0]))
ssf.ScaleHistogramsByHistogram("sf",copy.deepcopy(loosest_sf),divide=True)

ssf.CalculateHistograms("data")

ssf.PlotEfficienciesAndSFs(x_label=args.var_name,folder=args.output_folder)

ssf.ConvertToWPBinnedHistograms()

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

ssf.ScaleWPHistogramsByHistogram("sf",loosest_sf)

with open(args.wp_scores) as json_file: score_dict = json.load(json_file, object_pairs_hook=OrderedDict)

ssf.ConvertToScoreBinnedHistograms(score_dict)

#ssf.FitSpline("sf")

ssf.PlotSFs(x_label="Score",with_spline=False)

tf2 = ssf.DumpSFTF2()
print "SF formula:"
tf2.Print("all")

### Write to file ###
fout = ROOT.TFile(args.output_folder+"/"+args.name, 'RECREATE')
tf2.Write()
fout.Close()
print "Created {}/{}".format(args.output_folder,args.name)
