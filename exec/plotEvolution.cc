#include "../src/PlotUtil.h"
#include "../src/TemplateUtil.h"
#include "TFile.h"


// usage: root -l -b 'plotEvolution.cc++("input.root", "output.pdf", "op", <isData>, <addCentral>)'
// evolution plot where the evolution is directly in the exec macro instead of manual fits
//for iOp in dt cmm cmp cvv cva cav c1 c3 c123 root -l -b 'plotEvolution.cc++("../fit_output/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/TotalStatCovMatrix_AllVarNorm_rebinnedA/'${iOp}'/'${iOp}'_evolution.root","../fit_output/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/TotalStatCovMatrix_AllVarNorm_rebinnedA/'${iOp}'_evolution.pdf","'${iOp}'",false,false)'; sleep 0.1; done
// for iOp in ut dt cmm cmp cvv cva cav c1 c3 c123; do root -l -b '../../../exec/plotEvolution.cc++("'${iOp}'/'${iOp}'_evolution.root", "'${iOp}'_evolution.pdf", "'${iOp}'", false, false)'; sleep 0.1; done
//for iOp in ut dt cmm cmp cvv cva cav c1 c3 c123; do root -l -b 'plotEvolution.cc++("../fit_output/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/Systematics_AllVars_1D_228x228_1000PE/TotalStatCovMatrix_AllVarNorm_rebinnedA/'${iOp}'/'${iOp}'_evolution.root","../evolution_output/'${iOp}'_stat_evolution.pdf","'${iOp}'",false,false)'; done


//for iOp in ut dt cmm cmp cvv cva cav c1 c3 c123; do root -l -b 'plotEvolution.cc++("../fit_output/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/Systematics_AllVars_1D_132x132/TotalStatCovMatrix_AllVarNorm_rebinnedA/'${iOp}'/'${iOp}'_evolution.root","../fit_output/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/Systematics_AllVars_1D_132x132/TotalStatCovMatrix_AllVarNorm_rebinnedA/'${iOp}'/'${iOp}'_evolution.pdf","'${iOp}'",false,false)'; done
void plotEvolution(const std::string &fileName = "", const std::string &pdfName = "./evo.pdf",
                   const std::string &opName = "", const bool isData = false, const bool addCentral = false) {
  gROOT->Reset();
  TGaxis::SetMaxDigits(2);

  if (fileName == "" or opName == "")
    gROOT->ProcessLine(".q");

  // load file, initial setup
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));
  const bool isAll = (std::string(file->GetListOfKeys()->At(0)->GetName()).find("_all_") != std::string::npos);

  // map for the plot settings - key is opName, value is y axis range
  std::map<std::string, std::array<double, 2>> m_plot_config;
  m_plot_config.insert({"ctG", {0.0001, 0.4999}});
  m_plot_config.insert({"ut", {0., 0.03}});
  m_plot_config.insert({"dt", {0.0, 0.03}});
  m_plot_config.insert({"cmm", {0.0, 0.08}});
  m_plot_config.insert({"cmp", {0., 0.02}});
  m_plot_config.insert({"cvv", {0., 0.06}});
  m_plot_config.insert({"cva", {0., 0.06}});
  m_plot_config.insert({"cav", {0., 0.120}});
  m_plot_config.insert({"c1", {0., 0.7}});
  m_plot_config.insert({"c3", {0., 1.}});
  m_plot_config.insert({"c123", {0., 0.4}});
  m_plot_config.insert({"cnn", {0., 0.1}});

  std::map<std::string, std::string> root_ope_config;
  root_ope_config.insert({"ctG", "c_{tG} / {#Lambda}^2"});
  root_ope_config.insert({"ut", "#hat#mu_{t}"});
  root_ope_config.insert({"dt", "#hatd_{t}"});
  root_ope_config.insert({"cmm", "#hat c_{--}"});
  root_ope_config.insert({"cmp", "#hat c_{-+}"});
  root_ope_config.insert({"cvv", "#hat c_{VV}"});
  root_ope_config.insert({"cva", "#hat c_{VA}"});
  root_ope_config.insert({"cav", "#hat c_{AV}"});
  root_ope_config.insert({"c1", "#hat c_{1}"});
  root_ope_config.insert({"c3", "#hat c_{3}"});
  root_ope_config.insert({"c123", "#hat{c}_{1}-#hat{c}_{2}+#hat{c}_{3}"});
  root_ope_config.insert({"cnn", "#hat{C}_{nn}"});





  const double &yMin = m_plot_config.at(opName).at(0), yMax = m_plot_config.at(opName).at(1);

  // list of variables
  const std::vector<std::string> v_hist_str = {"b1k", "b2k", "b1r", "b2r", "b1n", "b2n", "b1j", "b2j", "b1q", "b2q",
                                           "ckk", "crr", "cnn",
                                           "cP_rk", "cM_rk", "cP_nr", "cM_nr", "cP_nk", "cM_nk",
                                           "ckj","crq",
                                           "chan","ctra","csca",
                                           "cPrj","cM_rj",
                                           "c_kjL", "c_rqL",
                                           "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM",
                                           "cHel"};

  const std::vector<std::string> v_bin_str = {"b^{1}_{k}", "b^{2}_{k}", "b^{1}_{r}", "b^{2}_{r}", "b^{1}_{n}", "b^{2}_{n}",
                                              "b^{1}_{k*}", "b^{2}_{k*}", "b^{1}_{r*}", "b^{2}_{r*}",
                                              "c_{kk}", "c_{rr}", "c_{nn}",
                                              "c_{rk} + c_{kr}", "c_{rk} - c_{kr}",
                                              "c_{nr} + c_{rn}", "c_{nr} - c_{rn}",
                                              "c_{nk} + c_{kn}", "c_{nk} - c_{kn}",
                                              "c_{kk*}","c_{rr*}",
                                              "c_{han}","c_{tra}","c_{sca}",
                                              "c_{rk*} + c_{k*r}","c_{rk*} - c_{k*r}",
                                              "c_{kk*L}", "c_{rr*L}",
                                              "c_{rk+}", "c_{rk-}", "c_{nr+}", "c_{nr-}", "c_{nk+}", "c_{nk-}",
                                              "c_{hel}"};
                                              //"c_{lab}", "#Delta#phi"
  /*
  const std::vector<std::string> v_bin_str = {"B^{1}_{k}", "B^{2}_{k}", "B^{1}_{r}", "B^{2}_{r}", "B^{1}_{n}", "B^{2}_{n}",
                                              "B^{1}_{k*}", "B^{2}_{k*}", "B^{1}_{r*}", "B^{2}_{r*}",
                                              "C_{kk}", "C_{rr}", "C_{nn}",
                                              "C_{rk} + C_{kr}", "C_{rk} - C_{kr}",
                                              "C_{nr} + C_{rn}", "C_{nr} - C_{rn}",
                                              "C_{nk} + C_{kn}", "C_{nk} - C_{kn}",
                                              "D"};
                                              //"A_{cos #varphi}^{lab}", "A_{|#Delta#phi_{ll}|}"
  */
  // get the indices of the variables that are the best at each iter
  int iPos = 0;
  std::vector<int> v_index;
  std::vector<std::array<double, 3>> v_band;

  for (uint iIter = 0; iIter < v_hist_str.size(); ++iIter) {
    for (int iK = iPos; iK < file->GetListOfKeys()->GetSize(); ++iK) {
      std::string name = file->GetListOfKeys()->At(iK)->GetName(), exp = "iter_" + toStr(iIter) + "_best";
      if (name.find(exp) == std::string::npos) continue;
      if (name.find("_sigma2") != std::string::npos) continue;

      exp = (isAll) ? exp + "_all_" : exp + "_linear_";

      replace(name, exp, "");
      replace(name, "_sigma1", "");

      v_index.push_back( std::stoi(name) );
      iPos = iK;
      break;
    }
  }

  // having known this we just grab their band widths
  for (uint iI = 0; iI < v_index.size(); ++iI) {
    std::string key = "iter_" + toStr(iI) + "_best";
    key = (isAll) ? key + "_all_" + toStr(v_index.at(iI)) : key + "_linear_" + toStr(v_index.at(iI));

    auto g_src1 = std::unique_ptr<TGraphAsymmErrors>(dynamic_cast<TGraphAsymmErrors *>(( file->Get( (key + "_sigma1").c_str() ))->Clone() ));
    auto g_src2 = std::unique_ptr<TGraphAsymmErrors>(dynamic_cast<TGraphAsymmErrors *>(( file->Get( (key + "_sigma2").c_str() ))->Clone() ));

    double best_fit, dummy;
    g_src1->GetPoint(0, best_fit, dummy);
    const double sig1U = g_src1->GetErrorXhigh(0), sig1D = g_src1->GetErrorXlow(0);
    const double sig2U = g_src2->GetErrorXhigh(0), sig2D = g_src2->GetErrorXlow(0);

    std::cout << "Best in iteration " << iI + 1 << ": " << v_hist_str.at(v_index.at(iI)) << ", best fit value " << best_fit
              << ", 1 sigma width " << sig1U + sig1D << ", 2 sigma width " << sig2U + sig2D << " at index " << v_index.at(iI) << std::endl;

    v_band.push_back({best_fit, sig1U + sig1D, sig2U + sig2D});
  }

  const uint nPnt = v_band.size();

  std::vector<double> v_sig0, v_sig1, v_sig2, v_step = {0.5};
  for (uint iB = 0; iB < v_band.size(); ++iB) {
    v_sig0.push_back( v_band.at(iB).at(0) );
    v_sig1.push_back( v_band.at(iB).at(1) );
    v_sig2.push_back( v_band.at(iB).at(2) );

  }
  fillInterval(v_step, v_step.back(), -0.5 + nPnt, 1.);

  auto g_sig0 = std::make_unique<TGraph>(nPnt, v_step.data(), v_sig0.data());
  stylePlot(g_sig0.get(), kBlack, 1., 0, 0, 1.5, 1, 3);

  auto g_sig1 = std::make_unique<TGraph>(nPnt, v_step.data(), v_sig1.data());
  stylePlot(g_sig1.get(), kGreen + 1, 1., 0, 0, 1.5, 1, 3);

  auto g_sig2 = std::make_unique<TGraph>(nPnt, v_step.data(), v_sig2.data());
  stylePlot(g_sig2.get(), kOrange, 1., 0, 0, 1.5, 1, 3);

  g_sig1->GetHistogram()->GetXaxis()->Set(nPnt, 0., nPnt);

  for (uint iP = 0; iP < nPnt; ++iP) {
    g_sig1->GetHistogram()->GetXaxis()->SetBinLabel(iP + 1, v_bin_str.at( v_index.at(iP) ).c_str());
    g_sig1->GetHistogram()->GetXaxis()->LabelsOption("v");

  }

  axisPlot(g_sig1.get(),
           yMin, yMax,510, "Confidence interval width", 0.043, 1.13, 0.037,
           0., 0.,0, "", 0.049, 1.01, 0.049);

// void axisPlot(Plot *plot,
//               const double yMin, const double yMax, const int yDiv,
//               const std::string &yTxt, const double ySiz, const double yOff, const double yLab,
//               const double xMin, const double xMax, const int xDiv,
//               const std::string &xTxt, const double xSiz, const double xOff, const double xLab)

  g_sig1->GetXaxis()->SetTickSize(0.01);

  // make the legend
  auto leg = std::make_unique<TLegend>();
  if (addCentral)
    leg->AddEntry(g_sig0.get(), "Best fit", "l");
  leg->AddEntry(g_sig1.get(), "68\% CL", "l");
  leg->AddEntry(g_sig2.get(), "95\% CL", "l");

  // making the canvas
  setH1Style();

  auto can = std::make_unique<TCanvas>("can", "can", 200, 10, 1000, 1000);
  can->SetTopMargin(0.075);
  can->SetBottomMargin(0.14);
  can->SetLeftMargin(0.10);
  can->SetRightMargin(0.025);

  // for tacking on some text on the plot
  const std::string &legOperator =  root_ope_config.at(opName);
  const std::string topMid = "EFT coupling: "+ legOperator;
  //const std::string topMid = "EFT Coupling: "+ legOperator + ", #Lambda = 1 TeV";
  //const std::string topLeft = "#scale[1.2]{#bf{CMS}}";
  const std::string topLeft = "#bf{CMS} #it{Preliminary}";
  const std::string topRight = "35.9 fb^{-1} (13 TeV)";
  TLatex txt;
  //txt.SetTextSize(0.039); // standard
  txt.SetTextSize(0.043); // cms
  txt.SetTextAlign(13);

  can->cd();

  ;
  const std::string legHeader = (isData) ? ( "Observed Data") : (" SM expectation");
  const double legYOffset = (addCentral) ? 0.065 : 0.;
  styleLegend(leg.get(), 1, 0, 0, 42, 0.043, legHeader);
  putLegend(leg.get(), 0.635, 0.915, 0.655 - legYOffset, 0.815 + legYOffset); // top right

  g_sig1->Draw("a l");

  leg->Draw();
  txt.DrawLatexNDC(0.131, 0.895, topLeft.c_str());
  txt.DrawLatexNDC(0.5, 0.895, topMid.c_str());
  txt.DrawLatexNDC(0.676, 0.975, topRight.c_str());

  if (addCentral)
    g_sig0->Draw("l");
  g_sig1->Draw("l");
  g_sig2->Draw("l");

  can->RedrawAxis();
  can->SaveAs(pdfName.c_str());

  gROOT->ProcessLine(".q");
}
