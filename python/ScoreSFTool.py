import ROOT
import json
from collections import OrderedDict
from array import array
import copy
ROOT.gROOT.SetBatch(1)

class ScoreSFTool:

  def __init__(self):

    self.histograms_od = OrderedDict()
    self.spline_od = OrderedDict()
    self.sorted_wp = list()
    self.rebinned_bins = OrderedDict()

  def CheckType(self,input_type):
    if input_type not in ["mc","data","sf"]:
      print "ERROR: Input type must be either mc, data or sf"
      return False
    else:
      return True
      
  def InputHistograms(self,hists,input_type):
  
    if not self.CheckType(input_type): return None

    self.histograms_od[input_type] = copy.deepcopy(hists)
    self.GetSortedWPs(hists)
    for k,v in self.histograms_od.iteritems():
      for wp in v.keys():
        if wp not in self.sorted_wp:
          print "WARNING: {} not in {}, removing this bin".format(wp,k)
          del self.histograms_od[k][wp]
      for wp in self.sorted_wp:
        if wp not in v.keys():
          print "WARNING: {} not in {}, removing this bin".format(wp,input_type)
          del self.histograms_od[input_type][wp]
       
  def GetSortedWPs(self,hists):

    loose_keys = list()
    medium_keys = list()
    tight_keys = list()
    for k in hists.keys():
      if "Loose" in k:
        loose_keys.append(k)
      elif "Medium" in k:
        medium_keys.append(k)
      elif "Tight" in k:
        tight_keys.append(k)

    loose_keys = sorted(loose_keys)
    loose_keys.reverse()
    tight_keys = sorted(tight_keys)

    self.sorted_wp = loose_keys + medium_keys + tight_keys

  def CalculateHistograms(self,output_type):
  
    if not self.CheckType(output_type): return None

    search_keys = ["mc","data","sf"]
    search_keys.remove(output_type) 
    ret = False
    for sk in search_keys:
      if sk not in self.histograms_od.keys():
        print "ERROR: {} not inputted".format(sk)
        ret = True
    if ret: return None

    if output_type == "mc":
      self.histograms_od["mc"] = OrderedDict()
      for k,v in self.histograms_od["data"].iteritems():
        self.histograms_od["mc"][k] = copy.deepcopy(v)
        self.histograms_od["mc"][k].Divide(self.histograms_od["sf"][k])
    elif output_type == "data":
      self.histograms_od["data"] = OrderedDict()
      for k,v in self.histograms_od["mc"].iteritems():
        self.histograms_od["data"][k] = copy.deepcopy(v)
        self.histograms_od["data"][k].Multiply(self.histograms_od["sf"][k])
    elif output_type == "sf":
      self.histograms_od["sf"] = OrderedDict()
      for k,v in self.histograms_od["data"].iteritems():
        self.histograms_od["sf"][k] = copy.deepcopy(v)
        self.histograms_od["sf"][k].Divide(self.histograms_od["mc"][k])

  def PrintHistograms(self,print_type=["mc","data","sf"]):
    if type(print_type) == str: print_type = [print_type]
    
    for tk, tv in self.histograms_od.iteritems():    
      if tk in print_type:
        print "<< {} >>".format(tk)
        for hk, hv in tv.iteritems():
          hv.Print("all")

  def ScaleHistogramsByHistogram(self,input_type,hist, divide=False):

    for k,v in self.histograms_od[input_type].iteritems():
      if not divide:
        self.histograms_od[input_type][k].Multiply(hist)
      else:
        self.histograms_od[input_type][k].Divide(hist)
    
  def GetHistograms(self,output_type):

    return copy.deepcopy(self.histograms_od[output_type])

  def GetHistogram(self,output_type,key):

    return copy.deepcopy(self.histograms_od[output_type][key])

  def ConvertToWPBinnedHistograms(self,rebin_threshold=0.1):
    
    temp_histograms_od = OrderedDict()
    temp_histograms_od["data"] = OrderedDict()
    temp_histograms_od["mc"] = OrderedDict()

    h = list(self.histograms_od["sf"].items())[0][1]
    bins = [[h.GetBinLowEdge(i),h.GetBinLowEdge(i+1)] for i in range(1,h.GetNbinsX()+1)] 
    
    # get data efficiencies first
    eff_bin_data_val = OrderedDict()
    eff_bin_wp = OrderedDict()
    for b in bins:
      bin_name = str(b[0]) + "to" + str(b[1])
      eff_bin_data_val[bin_name] = list()
      eff_bin_wp[bin_name] = list()
      for ind, wp in enumerate(self.sorted_wp):
        bin_number = self.histograms_od["data"][wp].FindBin((b[0]+b[1])/2)
        if ind + 1 != len(self.sorted_wp):
          eff_bin_data_val[bin_name].append(self.histograms_od["data"][wp].GetBinContent(bin_number) - self.histograms_od["data"][self.sorted_wp[ind+1]].GetBinContent(bin_number))
          eff_bin_wp[bin_name].append(self.sorted_wp[ind])
        else:
          eff_bin_data_val[bin_name].append(self.histograms_od["data"][wp].GetBinContent(bin_number))
          eff_bin_wp[bin_name].append(self.sorted_wp[ind])

    # rebin
    for k, v in eff_bin_data_val.iteritems():
      while not all(i >= rebin_threshold  for i in eff_bin_data_val[k]):
        for b, val in enumerate(eff_bin_data_val[k][:-1]):
          if val < rebin_threshold:
            eff_bin_data_val[k][b] += eff_bin_data_val[k][b+1]
            del eff_bin_data_val[k][b+1], eff_bin_wp[k][b+1]
            break

    # draw into histograms    
    temp_histograms_od = OrderedDict()
    for tk in ["data","mc"]:
      temp_histograms_od[tk] = OrderedDict()
      for b in bins:
        bin_name = str(b[0]) + "to" + str(b[1])
        wp_bins = array('f', map(float,range(0,len(eff_bin_wp[bin_name])+1)))
        wp_hist = ROOT.TH1D(bin_name,"",len(wp_bins)-1, wp_bins)
        temp_histograms_od[tk][bin_name] = copy.deepcopy(wp_hist)
        for ind, wp in enumerate(eff_bin_wp[bin_name]):
          bin_number = self.histograms_od[tk][wp].FindBin((b[0]+b[1])/2)
          if ind + 1 != len(eff_bin_wp[bin_name]):
            temp_histograms_od[tk][bin_name].SetBinContent(ind+1, self.histograms_od[tk][wp].GetBinContent(bin_number) - self.histograms_od[tk][eff_bin_wp[bin_name][ind+1]].GetBinContent(bin_number))
            temp_histograms_od[tk][bin_name].SetBinError(ind+1, (self.histograms_od[tk][wp].GetBinError(bin_number)**2 + self.histograms_od[tk][eff_bin_wp[bin_name][ind+1]].GetBinError(bin_number)**2)**0.5)
          else:
            temp_histograms_od[tk][bin_name].SetBinContent(ind+1, self.histograms_od[tk][wp].GetBinContent(bin_number))
            temp_histograms_od[tk][bin_name].SetBinError(ind+1, self.histograms_od[tk][wp].GetBinError(bin_number))

    self.histograms_od = copy.deepcopy(temp_histograms_od)
    self.CalculateHistograms("sf")
    self.rebinned_bins = copy.deepcopy(eff_bin_wp)

  def ScaleWPHistogramsByHistogram(self,input_type,hist,divide=False):

    for k,v in self.histograms_od[input_type].iteritems():
      scale_val = hist.GetBinContent(hist.FindBin((float(k.split("to")[0]) + float(k.split("to")[1]))/2))
      if not divide:
        self.histograms_od[input_type][k].Scale(scale_val)
      else:
        self.histograms_od[input_type][k].Scale(1.0/scale_val)
  
  def ConvertToScoreBinnedHistograms(self,score_dict):
    temp_histograms_od = OrderedDict()
    for tk, tv in self.histograms_od.iteritems():
      temp_histograms_od[tk] = OrderedDict()
      for bk, bv in tv.iteritems():
        bins = []
        for k in self.rebinned_bins[bk]: bins.append(score_dict[k])
        bins.append(1.0)
        score_bins = array('f', map(float,bins))
        temp_histograms_od[tk][bk] = copy.deepcopy(ROOT.TH1D(bv.GetName(),"",len(score_bins)-1, score_bins))
        for i in range(0,bv.GetNbinsX()+1):
          temp_histograms_od[tk][bk].SetBinContent(i,bv.GetBinContent(i))
          temp_histograms_od[tk][bk].SetBinError(i,bv.GetBinError(i))

    self.histograms_od = copy.deepcopy(temp_histograms_od)
        
  def FitSpline(self,input_type):
    self.spline_od[input_type] = OrderedDict()
    for tk, tv in self.histograms_od["sf"].iteritems():
      self.spline_od[input_type][tk] = ROOT.TSpline3(tv)

  def PlotEfficienciesAndSFs(self,logx=True,title_left="",title_right="",x_label="",ratio_range=[0.8,1.2],extra_ratio_line=[0.9,1.1],replace=[],replace_labels=[],individual=None,folder="."):

    for k,v in self.histograms_od["sf"].iteritems():
      if individual == None or individual == k:
        self.DrawHistogramsWithRatio(
                     [self.histograms_od["mc"][k],self.histograms_od["data"][k]],
                     ["MC","Data"],
                     y_label="Efficiency",
                     y_ratio_label="SF",
                     x_label=x_label,
                     title_left="{} {}".format(k,title_left),
                     title_right=title_right,
                     save_name="{}/effandsf_{}".format(folder,k),
                     logx=logx,
                     ratio_range=ratio_range,
                     extra_ratio_line=extra_ratio_line,
                     replace=replace,
                     replace_labels=replace_labels
                     )    

  def PlotSFs(self,logx=False,title_left="",title_right="",x_label="",with_spline=False,folder="."):
    for k,v in self.histograms_od["sf"].iteritems():
      if with_spline: spline = self.spline_od["sf"][k]
      else: spline = None
      self.DrawHistogramsWithRatio(
                   [self.histograms_od["sf"][k]],
                   [],
                   y_label="Scale Factor",
                   x_label=x_label,
                   title_left="{} {}".format(k,title_left),
                   title_right=title_right,
                   save_name="{}/sf_{}".format(folder,k),
                   logx=logx,
                   ratio=False,
                   fit=spline
                   )

  def DumpSFTF2(self):
    func_str = "("
    for k,v in self.histograms_od["sf"].iteritems():

      if k != self.histograms_od["sf"].keys()[-1]:
        func_str += "((x>{} && x<={})*(".format(k.split("to")[0],k.split("to")[1])
      else:
        func_str += "((x>{})*(".format(k.split("to")[0])

      for ind_wp, wp in enumerate(self.rebinned_bins[k]):
        if ind_wp + 1 != len(self.rebinned_bins[k]):
          func_str += "({}*(y>{} && y<={})) + ".format(v.GetBinContent(ind_wp+1),v.GetBinLowEdge(ind_wp+1),v.GetBinLowEdge(ind_wp+2))
        else:
          func_str += "({}*(y>{})))) + ".format(v.GetBinContent(ind_wp+1),v.GetBinLowEdge(ind_wp+1))

    func_str = func_str[:-3]+")"
    
    tf2 = ROOT.TF2("sf",func_str,0.0,1000.0,0.0,1.0)
    return tf2

  def DrawHistogramsWithRatio(self, hists, titles,x_label="",y_label="",y_ratio_label="",colours=[2,6,42,46,39,49], title_left="", title_right="",ratio_range=[0.8,1.2], extra_ratio_line=[0.9,1.1], save_name="plot", anchor_to_zero=True, logx=True, logy=False,replace=[],replace_labels=[],do_split_ratio_uncert=False,ratio=True,fit=None):

    c = ROOT.TCanvas('c','c',600,600)

    if ratio:
      c.SetBottomMargin(0.)
      pad1 = ROOT.TPad("pad1","pad1",0,0.47,1,1)
      pad1.SetBottomMargin(0.15)
      pad1.SetLeftMargin(0.15)
      if logx: pad1.SetLogx()
      if logy: pad1.SetLogy()
      pad1.Draw()
      pad1.cd()
    else:
      c.SetBottomMargin(0.15)
      c.SetLeftMargin(0.15)
      if logx: c.SetLogx()
      if logy: c.SetLogy()

    hists[0].Draw("E")
    hists[0].SetLineColor(1)
    hists[0].SetFillColor(1)
    hists[0].SetStats(0)
    hists[0].GetXaxis().SetTitle("")
    hists[0].GetYaxis().SetTitle(y_label)
    if ratio:
      hists[0].GetYaxis().SetTitleOffset(0.8)
      hists[0].GetYaxis().SetTitleSize(0.08)
      hists[0].GetYaxis().SetLabelSize(0.07)
    else:
      hists[0].GetYaxis().SetTitleOffset(1.5)
      hists[0].GetYaxis().SetTitleSize(0.045)
      hists[0].GetYaxis().SetLabelSize(0.04)
      hists[0].GetXaxis().SetTitle(x_label)
      hists[0].GetXaxis().SetTitleSize(0.045)
      hists[0].GetXaxis().SetTitleOffset(1.3)
    hists[0].SetTitle("")

    if ratio: hists[0].GetXaxis().SetLabelSize(0)

    for ind,val in enumerate(hists[1:]):
      hists[ind+1].Draw("E SAME")
      hists[ind+1].SetMarkerColor(colours[ind])
      hists[ind+1].SetLineColor(colours[ind])
      hists[ind+1].SetLineWidth(2)
      hists[ind+1].SetMarkerStyle(2)

    if fit != None: 
      fit.SetLineColor(2)
      fit.Draw("same")

    minimum = 999999.0
    maximum = 0.0
    for ind, h in enumerate(hists):
      for i in range(0,h.GetNbinsX()+1):
        if h.GetBinError(i) != 0 and h.GetBinContent(i)-h.GetBinError(i) < minimum:
          minimum = h.GetBinContent(i)-h.GetBinError(i)
        if h.GetBinError(i) != 0 and h.GetBinContent(i)+h.GetBinError(i) > maximum:
          maximum = h.GetBinContent(i)+h.GetBinError(i)

    hists[0].SetMaximum(1.5*maximum)
    if anchor_to_zero and not logy:
      hists[0].SetMinimum(0.0)
    else:
      hists[0].SetMinimum(0.8*minimum)

    if titles != []:
      l = ROOT.TLegend(0.7,0.65,0.88,0.85);
      l.SetBorderSize(0)
      for ind,h in enumerate(hists):
        l.AddEntry(hists[ind],titles[ind],"lep")
      l.Draw()


    if ratio:
      c.cd()
      pad2 = ROOT.TPad("pad2","pad2",0,0.0,1,0.45)
      pad2.SetLeftMargin(0.15)
      pad2.SetTopMargin(0.03)
      pad2.SetBottomMargin(0.55)
      if logx: pad2.SetLogx()
      pad2.Draw()
      pad2.cd()
  
      hist_divide = hists[0].Clone()
      if do_split_ratio_uncert:
        for i in range(0,hist_divide.GetNbinsX()+2): hist_divide.SetBinError(i,0)
  
      hists_ratio = []
      for ind,h in enumerate(hists):
        hists_ratio.append(hists[ind].Clone())
        hists_ratio[ind].Divide(hist_divide)
  
      hists_ratio[0].SetTitle("")
      hists_ratio[0].SetMarkerSize(0)
      hists_ratio[0].SetFillColorAlpha(12,0.5)
      hists_ratio[0].SetLineWidth(0)
      hists_ratio[0].SetAxisRange(ratio_range[0],ratio_range[1],'Y')
      hists_ratio[0].GetYaxis().SetNdivisions(4)
      hists_ratio[0].SetStats(0)
      hists_ratio[0].GetXaxis().SetLabelSize(0.08)
      hists_ratio[0].GetYaxis().SetLabelSize(0.08)
      hists_ratio[0].GetXaxis().SetTitle(x_label)
      hists_ratio[0].GetYaxis().SetTitle(y_ratio_label)
      hists_ratio[0].GetYaxis().SetTitleColor(1)
      hists_ratio[0].GetYaxis().SetTitleSize(0.1)
      hists_ratio[0].GetYaxis().SetTitleOffset(0.6)
      hists_ratio[0].GetXaxis().SetTitleSize(0.1)
      hists_ratio[0].GetXaxis().SetTitleOffset(1.2)
      if logx:
        hists_ratio[0].GetXaxis().SetMoreLogLabels()
        hists_ratio[0].GetXaxis().SetNoExponent()
  
      if do_split_ratio_uncert: 
        hists_ratio[0].Draw("e2")
      else:
        hists_ratio[0].Draw("l")
  
      if len(replace) != 0:
        self.ChooseAxisLabels(hists_ratio[0],replace,size=0.06,logx=logx,pad=pad2,offset=0.04,replace=replace_labels)
  
      for ind,h in enumerate(hists_ratio[1:]):
        hists_ratio[ind+1].SetLineColor(colours[ind])
        hists_ratio[ind+1].SetMarkerColor(colours[ind])
        hists_ratio[ind+1].SetMarkerStyle(2)
        hists_ratio[ind+1].SetLineWidth(1)
        hists_ratio[ind+1].Draw("E same")
  
      ratio_line = ROOT.TLine(hists[0].GetBinLowEdge(1),1,hists[0].GetBinLowEdge(hists[0].GetNbinsX()+1),1)
      ratio_line.SetLineStyle(3)
      ratio_line.Draw("l same")
  
      extra_ratio_l = []
      for ind, i in enumerate(extra_ratio_line):
        extra_ratio_l.append(ROOT.TLine(hists[0].GetBinLowEdge(1),i,hists[0].GetBinLowEdge(hists[0].GetNbinsX()+1),i))
        extra_ratio_l[ind].SetLineStyle(3)
        extra_ratio_l[ind].Draw("l same")

    if ratio:
      self.DrawTitle(pad1, title_left, 1, scale=1)
      self.DrawTitle(pad1, title_right, 3, scale=1)
    else:
      self.DrawTitle(c, title_left, 1, scale=.7)
      self.DrawTitle(c, title_right, 3, scale=.7)

    c.Update()
    c.SaveAs(save_name+".pdf")
    c.Close()

  def DrawTitle(self, pad, text, align, scale=1):
    pad_backup = ROOT.gPad
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()

    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if pad_ratio < 1.:
        pad_ratio = 1.

    textSize = 0.6
    textOffset = 0.2

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(textSize * t * pad_ratio * scale)

    y_off = 1 - t + textOffset * t
    if align == 1:
        latex.SetTextAlign(11)
    if align == 1:
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
    if align == 2:
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
    if align == 3:
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()

  def ChooseAxisLabels(self,axis,labels,size=0.04,logx=False,pad=None,offset=0.04,replace=[],angle=315):
    if replace != [] and len(labels) != len(replace):
      print "Replace and labels are not the same length"
      return

    axis.GetXaxis().SetLabelSize(0)
    axis_min = axis.GetXaxis().GetXmin()
    axis_max = axis.GetXaxis().GetXmax()
    axis.SetNdivisions(len(labels))
    if pad==None:
      b = ROOT.gPad.GetBottomMargin()
      l = ROOT.gPad.GetLeftMargin()
      r = 1 -ROOT.gPad.GetRightMargin()
    else:
      b = pad.GetBottomMargin()
      l = pad.GetLeftMargin()
      r = 1 -pad.GetRightMargin()

    for ind, num in enumerate(labels):
      latex = ROOT.TLatex()
      latex.SetNDC()
      latex.SetTextAngle(angle)
      latex.SetTextColor(ROOT.kBlack)
      latex.SetTextFont(42)
      latex.SetTextSize(size)
      latex.SetTextAlign(11)
      if not logx:
        x_shift = ((num - axis_min)*(r - l))/(axis_max - axis_min)
      else:
        x_shift = (math.log(num,10) - math.log(axis_min,10))*(r - l)/(math.log(axis_max,10) - math.log(axis_min,10))
      if replace != []:
        latex.DrawLatex(l+x_shift, b-offset, str(replace[ind]))
      else:
        latex.DrawLatex(l+x_shift, b-offset, str(num))
