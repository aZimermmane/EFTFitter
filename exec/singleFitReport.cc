#include "../src/PlotUtil.h"
#include "TFile.h"

// usage: root -l -b 'singleFitReport.cc++("/afs/desy.de/user/z/zimermma/work/EFTFitter/ut_1var_shape_sm.pdf", 0)'

// check if file exists - works only for full absolute path
bool is_nonexistent(const std::string& name) {
  return ifstream(name.c_str()).fail();
}

std::array<std::unique_ptr<TGraphAsymmErrors>, 3> makeConstraintGraph(const std::string &varName, const int nVar, const int nStep) {
  // read the graph and transform them to the single point
  const std::string fileName = "/afs/desy.de/user/z/zimermma/work/EFTFitter/output-2016-1D/sigleObsFit/Systematics_AllVars_1D_228x228_1000PE/TotalStatSystCovMatrix_AllVarNorm_rebinnedA/ut/ut_dchi2.root";

  // first check if file is there - MUST be direct path, method can't handle links
  const bool isNull = is_nonexistent(fileName);
  std::cout << "0" << std::endl;

  // make the graph for each point
  std::array<std::unique_ptr<TGraphAsymmErrors>, 3> a_graph = {nullptr, nullptr, nullptr};
  if (isNull)
    return a_graph;
    std::cout << "1" << std::endl;

  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));
  std::cout << "2" << std::endl;
  std::cout << "varName: " << varName << std::endl;
  // read the real best fit and sigmas
  //auto g_src1 = std::unique_ptr<TGraphAsymmErrors>(dynamic_cast<TGraphAsymmErrors *>(( file->Get(("linear_"+varName+"_sigma1").c_str() ))->Clone() ));
  //auto g_src2 = std::unique_ptr<TGraphAsymmErrors>(dynamic_cast<TGraphAsymmErrors *>(( file->Get(("linear_"+varName+"_sigma2").c_str() ))->Clone() ));
  if ( file->Get(("linear_"+varName+"_sigma1").c_str()) && file->Get(("linear_"+varName+"_sigma2").c_str()) ) {
  auto g_src1 = std::unique_ptr<TGraphAsymmErrors>(dynamic_cast<TGraphAsymmErrors *>(( file->Get(("linear_"+varName+"_sigma1").c_str() ))->Clone() ));
  auto g_src2 = std::unique_ptr<TGraphAsymmErrors>(dynamic_cast<TGraphAsymmErrors *>(( file->Get(("linear_"+varName+"_sigma2").c_str() ))->Clone() ));

  double best_fit = 0., dummy = 0.;
  g_src1->GetPoint(0, best_fit, dummy);
  const double sig1U = g_src1->GetErrorXhigh(0), sig1D = g_src1->GetErrorXlow(0);
  const double sig2U = g_src2->GetErrorXhigh(0), sig2D = g_src2->GetErrorXlow(0);
  std::cout << "3" << std::endl;

  a_graph.at(0) = std::make_unique<TGraphAsymmErrors>(1);
  a_graph.at(0)->SetName((varName + "_sig2_" + toStr(nVar) + "_" + toStr(nStep)).c_str());
  a_graph.at(0)->SetPoint(0, nVar, best_fit + nStep);
  a_graph.at(0)->SetPointError(0, 0., 0., sig2D, sig2U);
  stylePlot(a_graph.at(0).get(), kOrange, 1., 1001, 0, 1.5, 1, 3);
  std::cout << "4" << std::endl;

  a_graph.at(1) = std::make_unique<TGraphAsymmErrors>(1);
  a_graph.at(1)->SetName((varName + "_sig1_" + toStr(nVar) + "_" + toStr(nStep)).c_str());
  a_graph.at(1)->SetPoint(0, nVar, best_fit + nStep);
  a_graph.at(1)->SetPointError(0, 0., 0., sig1D, sig1U);
  stylePlot(a_graph.at(1).get(), kGreen + 1, 1., 1001, 0, 1.5, 1, 3);

  a_graph.at(2) = std::make_unique<TGraphAsymmErrors>(1);
  a_graph.at(2)->SetName((varName + "_best_" + toStr(nVar) + "_" + toStr(nStep)).c_str());
  a_graph.at(2)->SetPoint(0, nVar, best_fit + nStep);
  a_graph.at(2)->SetPointError(0, 0., 0., 0., 0.);
  stylePlot(a_graph.at(2).get(), kBlack, 1., 0, 0, 1.5, 1, 2);
  };
  return a_graph;

}

std::string variableString(const std::vector<std::string> &v_src, std::vector<int> &v_ind, const int index) {
  // strings up the variables in the correct permutation (without check of misuse)
  if (v_ind.empty()) return v_src.at(index);

  std::sort(std::begin(v_ind), std::end(v_ind));
  const auto iUp = std::upper_bound(std::begin(v_ind), std::end(v_ind), index);
  std::string varStr;

  // case where index is smaller than all
  if (iUp == std::begin(v_ind)) {
    varStr = v_src.at(index);
    for (long unsigned int iI = 0; iI < v_ind.size(); ++iI)
      varStr = varStr + "_" + v_src.at(v_ind.at(iI));

    return varStr;
  }

  // case where index is larger than all
  if (iUp == std::end(v_ind)) {
    varStr = v_src.at(v_ind.at(0));
    for (long unsigned int iI = 1; iI < v_ind.size(); ++iI)
      varStr = varStr + "_" + v_src.at(v_ind.at(iI));

    return varStr + "_" + v_src.at(index);
  }

  // ok fine it's somewhere in the middle
  const int nLo = std::distance(std::begin(v_ind), iUp);
  varStr = v_src.at(v_ind.at(0));
  for (int iI = 1; iI < nLo; ++iI)
    varStr = varStr + "_" + v_src.at(v_ind.at(iI));

  varStr = varStr + "_" + v_src.at(index);

  for (long unsigned int iI = nLo; iI < v_ind.size(); ++iI)
    varStr = varStr + "_" + v_src.at(v_ind.at(iI));

  return varStr;
}

std::array<double, 2> iterativeConstraint(const std::vector<std::string> &v_var, std::vector<int> &v_index, std::vector<int> &v_unsorted,
                                          std::vector<std::unique_ptr<TGraphAsymmErrors>> &v_graph, const std::string &best) {
  int iMin = -1;
  double minSig1 = 9999., minSig2 = 9999.;
  for (long unsigned int iV = 0; iV < v_var.size(); ++iV) {
    if (std::count(std::begin(v_index), std::end(v_index), iV)) continue;

    auto a_tmp = makeConstraintGraph(variableString(v_var, v_index, iV), (iV * 5) + 3, v_index.size());

    if (a_tmp.at(0) != nullptr) {
      for (auto &tmp : a_tmp)
        tmp->SetPointError(0, 2.999, 1.999, tmp->GetErrorYlow(0), tmp->GetErrorYhigh(0));

      if (a_tmp.at(0)->GetErrorYlow(0) + a_tmp.at(0)->GetErrorYhigh(0) < minSig2) {
        iMin = iV;
        minSig1 = a_tmp.at(1)->GetErrorYlow(0) + a_tmp.at(1)->GetErrorYhigh(0);
        minSig2 = a_tmp.at(0)->GetErrorYlow(0) + a_tmp.at(0)->GetErrorYhigh(0);
      }
    }

    std::move(std::begin(a_tmp), std::end(a_tmp), std::back_inserter(v_graph));
  }

  std::cout << "Best in iteration " << v_index.size() + 1 << ": " << v_var.at(iMin)
            << ", 1 sigma width " << minSig1 << ", 2 sigma width " << minSig2 << " at index " << iMin << std::endl;
  v_index.push_back(std::distance(std::begin(v_var), std::find(std::begin(v_var), std::end(v_var), best)));
  v_unsorted.push_back(v_index.back());

  return {minSig1, minSig2};
}

void singleFitReport(const std::string &pdfName, const uint showSnake = 1) {
  gROOT->Reset();
  TGaxis::SetMaxDigits(4);

  // start by making the snake graph at each step
  std::vector<std::unique_ptr<TGraphAsymmErrors>> v_graph;
  const int nStep = 1;

  double snakeSig1, snakeSig2;
  // make the snake graph and modify its x error to cover entire range
  std::string varAsSnake = (showSnake) ? "snake_spinHeli" : "cHel";
  auto a_tmp = makeConstraintGraph(varAsSnake, 55, 0);

  // and some manipulations to get these to look right
  for (auto &tmp : a_tmp)
    tmp->SetPointError(0, 110., 110., tmp->GetErrorYlow(0), tmp->GetErrorYhigh(0));
  stylePlot(a_tmp.at(0).get(), kGray + 1, 1., 1001, 0, 1.5, 1, 3);
  stylePlot(a_tmp.at(1).get(), kGray + 2, 1., 1001, 0, 1.5, 1, 3);
  stylePlot(a_tmp.at(2).get(), kBlack, 1., 0, 0, 1.5, 2, 1);

  snakeSig1 = a_tmp.at(1)->GetErrorYlow(0) + a_tmp.at(1)->GetErrorYhigh(0);
  snakeSig2 = a_tmp.at(0)->GetErrorYlow(0) + a_tmp.at(0)->GetErrorYhigh(0);

  std::move(std::begin(a_tmp), std::end(a_tmp), std::back_inserter(v_graph));

  // list of variables
  const std::vector<std::string> v_hStr = {"b1k", "b2k", "b1r", "b2r", "b1n", "b2n", "b1j", "b2j", "b1q", "b2q",
                                           "ckk", "crr", "cnn",
                                           "cP_rk", "cM_rk", "cP_nr", "cM_nr", "cP_nk", "cM_nk",
                                           "cHel","chan","ctra","csca",
                                           "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM"/*,
                                           "cLab", "LL_dPhi"*/};

  // make the graphs for each step - keep track of two vectors, one sorted one not...
  std::vector<int> v_hInd, v_nInd;
  std::vector<std::array<double, 2>> v_band;

  v_band.push_back( iterativeConstraint(v_hStr, v_hInd, v_nInd, v_graph, "cHel") );

  std::cout << "snake_spinCorr fit: 1 sigma width " <<  snakeSig1 << ", 2 sigma width " << snakeSig2 << std::endl;

  /*/ modify the axis of the first graph to string labels
  const std::vector<std::string> v_bStr = {"b_{1k}", "b_{2k}", "b_{1r}", "b_{2r}", "b_{1n}", "b_{2n}", "b_{1j}", "b_{2j}", "b_{1q}", "b_{2q}",
                                           "c_{kk}", "c_{rr}", "c_{nn}",
                                           "c_{rk} + c_{kr}", "c_{rk} - c_{kr}", "c_{nr} + c_{rn}", "c_{nr} - c_{rn}", "c_{nk} + c_{kn}", "c_{nk} - c_{kn}",
                                           "c_{hel}",
                                           "c_{lab}", "#Delta#phi"};

  const std::vector<std::string> v_bStr = {"b_{1k}", "b_{2k}", "b_{1r}", "b_{2r}", "b_{1n}", "b_{2n}", "b_{1k*}", "b_{2k*}", "b_{1r*}", "b_{2r*}",
                                           "c_{kk}", "c_{rr}", "c_{nn}",
                                           "c_{rk} + c_{kr}", "c_{rk} - c_{kr}", "c_{nr} + c_{rn}", "c_{nr} - c_{rn}", "c_{nk} + c_{kn}", "c_{nk} - c_{kn}",
                                           "cos #varphi",
                                           "cos #varphi_{lab}", "|#Delta#phi|"};

  const std::vector<std::string> v_bStr = {"B^{1}_{k}", "B^{2}_{k}", "B^{1}_{r}", "B^{2}_{r}", "B^{1}_{n}", "B^{2}_{n}", "B^{1}_{k*}", "B^{2}_{k*}", "B^{1}_{r*}", "B^{2}_{r*}",
                                           "C_{kk}", "C_{rr}", "C_{nn}",
                                           "C_{rk} + C_{kr}", "C_{rk} - C_{kr}", "C_{nr} + C_{rn}", "C_{nr} - C_{rn}", "C_{nk} + C_{kn}", "C_{nk} - C_{kn}",
                                           "D",
                                           "A_{cos #varphi}^{lab}", "A_{|#Delta#phi_{ll}|}"};
  */

  const std::vector<std::string> v_bStr = {"b^{1}_{k}", "b^{2}_{k}", "b^{1}_{r}", "b^{2}_{r}", "b^{1}_{n}", "b^{2}_{n}", "b^{1}_{k*}", "b^{2}_{k*}", "b^{1}_{r*}", "b^{2}_{r*}",
                                           "c_{kk}", "c_{rr}", "c_{nn}",
                                           "c_{rk} + c_{kr}", "c_{rk} - c_{kr}", "c_{nr} + c_{rn}", "c_{nr} - c_{rn}", "c_{nk} + c_{kn}", "c_{nk} - c_{kn}",
                                           "c_{hel}","c_{han}","c_{tra}","c_{sca}",
                                           "c_{rk+}", "c_{rk-}", "c_{nr+}", "c_{nr-}", "c_{nk+}", "c_{nk-}"/*,
                                           "A_{cos #varphi}^{lab}", "A_{|#Delta#phi_{ll}|}"*/};

  auto &axis = v_graph.at(2);
  axis->GetHistogram()->GetXaxis()->Set(5 * v_hStr.size(), 0., 5. * v_hStr.size());
  for (long unsigned int iV = 0; iV < 5 * v_hStr.size(); ++iV) {
    const std::string binStr = (iV % 5 != 2) ? "" : v_bStr.at(std::floor(iV / 5));
    axis->GetHistogram()->GetXaxis()->SetBinLabel(iV + 1, binStr.c_str());
  }

  // axisPlot(axis.get(),
  //          -1.999, 1.999, "#hat #mu _{t} / #Lambda^{2} [TeV^{-2}]", 0.057, 0.47, 0.067,
  //          0., 0., "", 0.053, 1.11, 0.053);
  // axis->GetXaxis()->SetNdivisions(-v_hStr.size());
  // axis->GetXaxis()->SetTickSize(0.01);

  axisPlot(axis.get(),
           -0.1999, 0.1999,510,
           "#hat#mu_{t}", 0.057, 0.47, 0.067,
           0., 0., 0,
           "", 0.053, 1.11, 0.053);
  //axis->GetXaxis()->SetNdivisions(-v_hStr.size());
  //axis->GetYaxis()->SetNdivisions(220);
  axis->GetXaxis()->SetTickSize(0.01);

  // need to find first non-null index
  int isFirstValid = -1;
  for (long unsigned int iV = 3 * nStep; iV < v_graph.size(); ++iV) {
    if (v_graph.at(iV) == nullptr) continue;
    isFirstValid = iV;
    break;
  }

  auto leg = std::make_unique<TLegend>();
  if (showSnake) leg->AddEntry(v_graph.at(2).get(), "Full best fit", "l");
  leg->AddEntry(v_graph.at(isFirstValid + 2).get(), "Best fit", "l");
  if (showSnake) leg->AddEntry(v_graph.at(1).get(), "68% CL", "f");
  leg->AddEntry(v_graph.at(isFirstValid + 1).get(), "68% CL", "f");
  if (showSnake) leg->AddEntry(v_graph.at(0).get(), "95% CL", "f");
  leg->AddEntry(v_graph.at(isFirstValid).get(), "95% CL", "f");

  // making the canvas
  setH1Style();

  auto can = std::make_unique<TCanvas>("can", "can", 200, 10, 1920, 1080);
  can->SetTopMargin(0.09);
  can->SetBottomMargin(0.13);
  can->SetLeftMargin(0.07);
  can->SetRightMargin(0.015);

  const std::string topRight = "35.9 fb^{-1} (13 TeV)";
  TLatex txt;
  txt.SetTextSize(0.063);
  txt.SetTextAlign(13);

  can->cd();

  if (showSnake) {
    styleLegend(leg.get(), 2, 0, 0, 42, 0.043, "");
    putLegend(leg.get(), 0.101, 0.520, 0.675, 0.985);
  }
  else {
    styleLegend(leg.get(), 3, 0, 0, 42, 0.053, "");
    putLegend(leg.get(), 0.101, 0.620, 0.775, 0.985);
  }

  axis->Draw("a l");

  leg->Draw();

  for (long unsigned int iV = 3 * nStep; iV < v_graph.size(); ++iV) {
    if (v_graph.at(iV) == nullptr) continue;

    if (iV % 3 != 2)
      v_graph.at(iV)->Draw("2");
    else
      v_graph.at(iV)->Draw("lz");
  }

  if (showSnake) {
    for (int iV = 0; iV < 3 * nStep; ++iV) {
      if (v_graph.at(iV) == nullptr) continue;

      if (iV % 3 != 2)
        v_graph.at(iV)->Draw("2");
      else
        v_graph.at(iV)->Draw("lz");
    }
  }

  txt.DrawLatexNDC(0.74, 0.99, topRight.c_str());

  can->RedrawAxis();
  can->SaveAs(pdfName.c_str());
  gROOT->ProcessLine(".q");
}
