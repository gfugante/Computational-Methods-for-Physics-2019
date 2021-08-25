


import ROOT
import numpy as np

PATH1 = "30min_b1E6"
PATH2 = "30min_2"
PATH3 = "30min_1E6"

data1 = np.loadtxt(PATH1+".txt", usecols=[0, 1, 2])
data2 = np.loadtxt(PATH2+".txt", usecols=[0, 1, 2])
data3 = np.loadtxt(PATH3+".txt", usecols=[0, 1, 2])

time1 = ROOT.TGraph()
time2 = ROOT.TGraph()
time3 = ROOT.TGraph()
error1 = ROOT.TGraph()
error2 = ROOT.TGraph()
error3 = ROOT.TGraph()

time1.SetTitle("Time")
error1.SetTitle("Error")

for i in range(len(data2.T[0])):
	time1.SetPoint(i, data1.T[0][i], data1.T[1][i])
	time2.SetPoint(i, data2.T[0][i], data2.T[1][i])
	time3.SetPoint(i, data3.T[0][i], data3.T[1][i])

	error1.SetPoint(i, data1.T[0][i], data1.T[2][i]/1E6)
	error2.SetPoint(i, data2.T[0][i], data2.T[2][i])
	error3.SetPoint(i, data3.T[0][i], data3.T[2][i])
	


c1 = ROOT.TCanvas("c1", "c1")
c2 = ROOT.TCanvas("c2", "c2")

legend = ROOT.TLegend(0.1,0.9,0.30,0.7)
legend.AddEntry(time1, PATH1)
legend.AddEntry(time2, PATH2)
legend.AddEntry(time3, PATH3)


c1.cd()

time2.SetMarkerColor(4)
time3.SetMarkerColor(2)
#time1.Fit("pol2")
time1.Draw("AP*")
time2.Draw("*SAME")
time3.Draw("*SAME")
legend.Draw("SAME")


c2.cd()

error2.SetMarkerColor(4)
error3.SetMarkerColor(2)
#error1.Fit("pol2")
error1.Draw("AP*")
error2.Draw("*SAME")
error3.Draw("*SAME")
legend.Draw("SAME")





raw_input("Press enter to close.")