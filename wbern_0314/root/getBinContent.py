import ROOT

# Open the ROOT file
file = ROOT.TFile("ut_coeff.root")

# Get the TH1D histogram by name (replace "your_histogram_name" with the actual name of your histogram)
hist = file.Get("sm_nominal")

# Loop over the bins of the histogram and print the bin number and its corresponding value
for i in range(1, hist.GetNbinsX()+1):
    print("Bin {}: {}".format(i, hist.GetBinContent(i)))

# Close the ROOT file
file.Close()
