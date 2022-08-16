import numpy as np
import ROOT
import scipy
from scipy import interpolate

import sys;
sys.path.append("JPyPlotRatio");
import JPyPlotRatio

fload = ROOT.TFile("output/h2dCorrFit.root","read"); #Opens figs

dataTypePlotParams = 	[
							{'plotType':'data','color':'k','fmt':'o','markersize':2.0},
							{'plotType':'theory','facecolor':'b','edgecolor':'b','alpha':0.3,'linestyle':'solid','linecolor':'b'},
							{'plotType':'data','color':'r','fmt':'s','markersize':1.0},
							{'plotType':'theory','facecolor':'y','edgecolor':'y','alpha':0.3,'linestyle':'solid','linecolor':'y'},
							{'plotType':'theory','facecolor':'g','edgecolor':'g','alpha':0.3,'linestyle':'solid','linecolor':'g'},
							{'plotType':'theory','facecolor':'r','edgecolor':'r','alpha':0.3,'linestyle':'solid','linecolor':'r'}
						];

xtitle = ["$\\Delta\\varphi (\\mathrm{rad})$"];
ytitle = ["$\\frac{1}{N_{\\mathrm{trig}}}\\frac{\\mathrm{d}N^{\\mathrm{pair}}}{\\mathrm{d}\\Delta\\varphi}$"];
nrow = 1;
ncol = 1;
TypeName =["Signal (0--0.1\%)", "Fit","$FY_{\\mathrm{LM}} + G$", "$G(1+2v_{2,2}cos(2\\Delta\\varphi))$ \n $+ FY_{\\mathrm{LM,min}}$","$G(1+2v_{3,3}cos(3\\Delta\\varphi))$ \n $+ FY_{\\mathrm{LM,min}}$"];
'''
# Limits the panel 
xlimits = {0:(-1.4,4.5)};
ylimits = {0:(0.903,0.99)};
'''
rlimits = [(0.8,1.2)];

plot = JPyPlotRatio.JPyPlotRatio(panels=(1,1),
	panelsize=(13,13), # change the size
	#panelLabel="",
	#rowBounds=ylimits,  # for nrow
	#colBounds=xlimits,  # for ncol
	ratioBounds=rlimits,# for nrow
	panelLabelLoc=(0.2,1.1),panelLabelSize=13,panelLabelAlign="left",
	legendPanel=0,
	legendLoc=(0.2,1.12),legendSize=11,xlabel=xtitle[0],ylabel=ytitle[0]); # x- and y-coordinate labels

fig = fload.Get("hDphiHM;1"); # HM
fit = fload.Get("fFit_best"); # fit
fit_v2 = fload.Get("fit_s_v2;1");# v22
fit_v3 = fload.Get("fit_s_v3;1"); # v33
fit_G = fload.Get("hY_a_G;1"); # 

data = plot.Add(0, fig, **dataTypePlotParams[0], label=TypeName[0]);
data_fit = plot.Add(0, fit, **dataTypePlotParams[1], label=TypeName[1]);
data_fit_v2 = plot.Add(0, fit_v2, **dataTypePlotParams[3], label=TypeName[2]);
data_fit_v3 = plot.Add(0, fit_v3, **dataTypePlotParams[4], label=TypeName[3]);
data_fit_G = plot.Add(0, fit_G, **dataTypePlotParams[2], label=TypeName[4]);

fload.Close();

plot.GetPlot().text(0.35,0.98,"Peripheral 2-particle correlation",fontsize=12);
plot.GetPlot().text(0.1,0.90,"pp $\\sqrt{s}$ = ",fontsize=11);
plot.GetPlot().text(0.1, 0.92,"$ < p_\\mathrm{T,trig(assoc)} < \\,\\mathrm{GeV}/c$",fontsize=10);
plot.GetPlot().text(0.5, 0.925,"$-4.0 < \\eta < 4.0 $");
plot.GetRatioAxes(0).xaxis.set_ticks_position('both');
plot.GetRatioAxes(0).yaxis.set_ticks_position('both');
plot.Ratio(data, data_fit); # Plots theory vs data ratio

plot.Plot();


# plot.Save("figs/Fig14_FlowExt.pdf");
# plot.Save("figs/Fig14_FlowExt.png");
plot.Show();