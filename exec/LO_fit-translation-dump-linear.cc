 // -*- C++ -*-
// example code using the EFTFitter plugin
// plugins rely on ROOT6 and C++14 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh ./exec/LO_fit-translation-dump-linear.cc
// version for doing regular fits ie absolute or shape
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/stdc++.h>
#include <iostream>
#include "../src/EFTFitter.h"

int main() {
  using EFT = EFTFitter;

  // construction done with the key we want to treat as data, lambda, fit mode and stat mode
  EFT eft("data", 1., EFT::Fit::shape, EFT::Stat::xsec);

  // these are just to avoid writing them out repeatedly
  //const std::string obsName = "cHel";
  const std::string sumName = "TTbarSpinDensityMatrix/sumWgt_noCut";
  //const std::string input_dir = "/nfs/dust/cms/user/zimermma/histograms/makeHistograms/coffeaSpin/";
  //const std::string input_dir = "/nfs/dust/cms/user/zimermma/histograms/eftsamples_ttbar_eft_fcnc_Breuther_dim6top_reweighting_13TeV_extra/";
  const std::string input_dir = "/nfs/dust/cms/user/zimermma/histograms/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/";
  //const std::string mcSample = "ttbareft_4f_Breuther_dim6top_reweighting_13TeV";
  //const std::string mcSample = "eftsamples_ttbar_eft_fcnc_Breuther_dim6top_reweighting_13TeV";
  const std::string mcSample = "ttbareft_translationp_Breuther_dim6top_reweighting_13TeV";
  const int nRebin = 4;
  //const std::string fitMode = eft.dumpFitMode();
  const double k_nnlo_lo = 1.667296656, br_tt_2l = 0.041062412; // NNLO/LO k-factor
  //const double k_nnlo_lo = 1.0 , br_tt_2l = 0.06635231524673364;
  //const double k_nnlo_lo = 1.0 , br_tt_2l = 1.0;
  //const std::vector<std::string> opNameStr = {/*"ctu8","ctq8","ctd8","cQu8","cQq83","cQq81",*/
  //                                            "cQd8"
  //                                          /*,"ctG","ctGI"*/};
  //const std::vector<std::string> opNameStr = { "cmp" ,"cmmfcnc"}
  const std::vector<std::string> opNameStr = {"cAVp","cVAp","cVVp","cAAp","c1p","c3p",
                                             "c1mc2pc3p","ctG","ctGI"};
  //const std::vector<std::string> opNameStr ={"ctG"};
  //const std::vector<std::string> opLatexStr ={"c_{tG}"};
  //const std::vector<std::string> opLatexStr = { "c_{-+}" ,"c_{--}"}
  const std::vector<std::string> opLatexStr = {"c'_{AV}","c'_{VA}","c'_{VV}","c'_{AA}","c'_{1}","c'_{3}",
                                             "c'_{1}-c'_{2}+c'_{3}","c_{tG}","c_{tG}^{I}"};
  //const std::vector<std::string>  WCvalueList = {"-2"};
  const std::vector<std::string> WCvalueList = {"-5","-2","-1","1","2","5"};
  //const std::vector<std::string> WCvalueList = {"-0.1","-0.01","0.01","0.1"};
  for (int i = 0; i < opNameStr.size(); ++i) {
    const std::string &opName = opNameStr.at(i);
    const std::string &opLatex = opLatexStr.at(i);
    const std::vector<std::string> spinVarStr = {"b1k", "b2k", "b1j", "b2j","b1r", "b2r", "b1q", "b2q", "b1n", "b2n","ckk", "ckj", "crr", "crq", "cnn", "cHel", "cHan","cSca",
    "cTra","cPrk","cMrk", "cPrj", "cMrj", "cPnr", "cMnr", "cPnk", "cMnk"};
                   //"bPkk","bMkk","bPjj","bMjj","bPrr","bMrr","bPqq","bMqq","bPnn","bMnn",
                   //"ckjL","crqL",
                   //"cPqk", "cMqk",
                   //"crkP", "crkM", "crjP", "crjM", "cqkP", "cqkM", "cnrP", "cnrM","cnkP", "cnkM"

    // //const std::vector<std::string> spinVarStr = {"b1x","b2x","b1y", "b2y","b1z", "b2z","b1u", "b2u","b1v", "b2v","b1w", "b2w", "cxx","cyy","czz","cxu","cyv","czw"};
    //const std::vector<std::string> spinVarStr = {"b1x","b2x","b1u", "b2u"};

    //const std::vector<std::string> spinVarStr = {"ttbar_mass_no_cut","ttbar_pt_no_cut", "lepton_pt_no_cut","lepton_eta_no_cut"};
    //const std::vector<std::string> spinVarStr = {"b1k"};

     const std::vector<std::string> spinVarLatex = {"b_{k}^{1}", "b_{k}^{2}","b_{k*}^{1}", "b_{k*}^{2}", "b_{r}^{1}", "b_{r}^{2}","b_{r*}^{1}", "b_{r*}^{2}", "b_{n}^{1}", "b_{n}^{2}","c_{kk}", "c_{kj}","c_{rr}","c_{rq}","c_{nn}","c_{Hel}","c_{Han}","c_{Sca}","c_{Tra}","c_{rk}+c_{kr}"
     ,"c_{rk}-c_{kr}","c_{rk*}+c_{k*r}", "c_{rk*}-c_{k*r}", "c_{nr}+c_{rn}", "c_{nr}-c_{rn}", "c_{nk}+c_{kn}", "c_{nk}-c_{kn}" };

                                                   //"b_{k}^{1}+b_{k}^{2}","b_{k}^{1}-b_{k}^{2}","b_{k*}^{1}+b_{k*}^{2}","b_{k*}^{1}-b_{k*}^{2}","b_{r}^{1}+b_{r}^{2}","b_{r}^{1}-b_{r}^{2}","b_{r*}^{1}+b_{r*}^{2}","b_{r*}^{1}-b_{r*}^{2}","b_{r}^{n}+b_{n}^{2}","b_{n}^{1}-b_{n}^{2}",


                                                   //,"c_{kjL}","c_{rqL}",
                                                   //"c_{r*k}+c_{kr*}", "c_{r*k}-c_{kr*}",
                                                   //"c_{rk}^{+}", "c_{rk}^{-}", "c_{rj}^{+}", "c_{rj}^{-}", "c_{qk}^{+}", "c_{qk}^{-}", "c_{nr}^{+}", "c_{nr}^{-}","c_{nk}^{+}", "c_{nk}^{-}"

    //const std::vector<std::string> spinVarLatex = {"b_{x}^{1}", "b_{x}^{2}","b_{y}^{1}", "b_{y}^{2}","b_{z}^{1}", "b_{z}^{2}","b_{x*}^{1}", "b_{x*}^{2}","b_{y*}^{1}", "b_{y*}^{2}","b_{z*}^{1}", "b_{z*}^{2}","c_{xx}","c_{yy}","c_{zz}","c_{xx*}","c_{yy*}","c_{zz*}"};
    //const std::vector<std::string> spinVarLatex = {"b_{x}^{1}", "b_{x}^{2}","b_{x*}^{1}", "b_{x*}^{2}"};
    //const std::vector<std::string> spinVarLatex = {"m_{tt}","pt_{tt}", "pt_{l}","eta_{l}"};
    //const std::vector<std::string> spinVarLatex = {"b_{k}^{1}"};
  //const std::vector<std::array<int, 2>> v_snake = {{1, 6}, {7, 12}, {13, 18}, {19, 24}, {25, 30}, {31, 36}, {37, 42},
    //                                      {43, 48}, {49, 54}, {55, 60}, {61, 66}, {67, 72}, {73, 78}, {79, 84},
    //                                      {85, 90}, {91, 96}, {97, 102}, {103, 108}, {109, 114},
    //                                      {115, 120}/*, {121, 126}, {127, 132}*/};
  //std::vector< std::pair<std::vector<std::string>, std::vector<std::array<int, 2>>>>
    //v_hName = {{"snake_spinHeli", v_snake}};
    //v_hName = {{"cHel", v_snake}};


  //std::mkdir(output_dir)


  // Creating a directory
  mkdir(("./fit_output/"+mcSample+"/").c_str(), 0777);
  mkdir(("./fit_output/"+mcSample+"/dumplinear").c_str(), 0777);
  if (mkdir(("./fit_output/"+mcSample+"/dumplinear/").c_str(), 0777) != -1) std::cout << "Directory created";
  if (mkdir(("./fit_output/"+mcSample+"/dumplinear/"+opName+"/").c_str(), 0777) != -1) std::cout << "Directory created";


  for (int iVar = 0; iVar < spinVarStr.size(); ++iVar) {
    const std::string &hName = spinVarStr.at(iVar);
    const std::string &hNameLatex =spinVarLatex.at(iVar);
    //const std::vector<std::array<int, 2>> &v_endbin = v_hName.at(iVar).second;

    for (int iVar  = 0; iVar  < WCvalueList.size(); ++iVar) {
      const std::string &WCvalue = WCvalueList.at(iVar);

    const std::string output_dir = "./fit_output/"+mcSample+"/dumplinear/"+opName+"/"+hName +"/";
    { // Creating a directory
    if (mkdir((output_dir).c_str(), 0777) != -1) std::cout << "Directory created";}
    { // Creating a directory
    if (mkdir((output_dir+"c"+WCvalue).c_str(), 0777) != -1) std::cout << "Directory created";}
    const std::string histName = "gen_" + hName;
    //const std::string histName = hName;


    std::cout << "EFTFitter: starting fit on variable " << hName << std::endl;
    //eft.setShapeSum(v_endbin.size());

  //const std::string histName = "TTbarSpinDensityMatrix/some_histogram", sumName = "TTbarSpinDensityMatrix/sumWgt_noCut";

  //const double k_nnlo_lo = 1.0 , br_tt_2l =1.0;
  //const double r = 831.76/20.4865028262;

  // add the input file and hist names (including if necessary the sum of weight hist for normalization)
  // hist name is such that file->Get(name) works
  // please ensure there is exactly 1 input with dataName as in ctor (assign 0 xsec, error to deactivate normalization)
  //eft.addRawInput("data", EFT::Sample::all, input_dir + "/unfolded_data_190114.root.root",
  //                histName, "", 1, {0., 0.}, EFT::Stat::xsec);

  // and MC following the syntax op1_val1--op2_val2-- ... --opN_valN
  // all operators to be considered must be present in key - ie dont write c2_4 when doing c1 c2 c3 fits, write c1_0--c2_4--c3_0
  // only the Sample::all types are considered for interpolation
  // xsec given is some dummy values (with k-factor applied), last arg stands for the kind of histogram given: xsec vs xsec
  // more explanation in header
  // SM

  //eft.setShapeSum(v_endbin.size());//added to try linear fit

  //eft.addRawInput("data", EFT::Sample::all, input_dir + "covariance_matrix/unfolded_data_190114.root",
  //               "TTbarSpinDensityMatrix/" + hName + "_absolute_data_a", "", 1, {0., 0.}, EFT::Stat::xsec);

    eft.addRawInput(opName+"_0", EFT::Sample::all, input_dir + mcSample+"_SM_hist.root",
                        histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);

  //eft.addRawInput("c1_0", EFT::Sample::all, input_dir + "SM_hist.root",
  //                      histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);

  //eft.addRawInput("cVA_0", EFT::Sample::all, input_dir + "CtgMCData/c1_0-NLO.root",
                                              //histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);

  //eft.addRawInput("c1_0", EFT::Sample::all, input_dir + "c1_0.root",
  //                                histName, sumName, nRebin, {(k_nnlo_lo / br_tt_2l) * 20.4847, 0.}, EFT::Stat::xsec);

  // 1D inputs - here c1 = 10 and c1 = 272 is chosen as raw inputs
  //eft.addRawInput("u_-0.15", EFT::Sample::all, input_dir + "ttbareft_translation2_dim6top_reweighting_13TeV_u_m0p15_hist.root",
  //                histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);

    eft.addRawInput(opName+"_-4", EFT::Sample::all, input_dir + mcSample+"_"+opName+"_m4_hist.root",
                                  histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);
    eft.addRawInput(opName+"_8", EFT::Sample::all, input_dir + mcSample+"_"+opName+"_8_hist.root",
                  histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);


  // prepare the base for interpolation ie compute individual contribution at 1
    eft.prepareInterpolationBase();

  // in case of fitting to MC
  // assign as data a particular key of choice
   eft.assignAsData(opName+"_0", EFT::Sample::all);

  /*/ data needs to be in the list to be drawn - in the braces are key, type and legend text
  std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
  vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});
  vt_keySampleLegend.push_back({"c1_1", EFT::Sample::all, "c1 = 1"});
  vt_keySampleLegend.push_back({"c1_0", EFT::Sample::all, "SM"});

  // args are just filename and customizing the axes and flags for dividing by bin width etc
  // also can control what to put in ratio plot (more detail in header)
  eft.drawHistogram(vt_keySampleLegend,
                    output_dir + "cQq11_snake_spinCorr", "#sigma [pb]", "Index", 0.0001, 599.9999, 0.5001, 1.9999, false, "", false, "simple");
  */

  // grab the total stat error matrix - as usual matrix name is such that file->Get("some_matrix_name") works
  // can also partially extract along the diagonal, pass a vector of bin index range eg {{1, 6}, {115, 120}} as last arg
  //ANDRE  eft.readCovMatRoot("finalcov", input_dir + "covariance_matrix/covmat_190114.root", "TotalStatSyst_"+fitMode+"_a",{{iVar*6+1, iVar*6+6}});

  // draw and store all available covmats as TH2 (can also provide a vector of keys to draw only specific ones)
  // also produces the matrix in text format
  //eft.drawCovMat(output_dir);

    // std::cout << "Measured rate in data is " << eft.getDataRate() << std::endl;
    //
        std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
        vt_keySampleLegend.push_back({opName+"_"+WCvalue, EFT::Sample::all, opLatex+"="+WCvalue+" All"});
        vt_keySampleLegend.push_back({opName+"_"+WCvalue, EFT::Sample::linear, opLatex+"="+WCvalue+" SM+Linear"});
        vt_keySampleLegend.push_back({opName+"_"+WCvalue, EFT::Sample::quadratic, opLatex+"="+WCvalue+" SM+EFT^{2}"});
        //vt_keySampleLegend.push_back({opName+"_"+WCvalue2, EFT::Sample::all, "c'_{3} ="+WCvalue2+"All"});
        //vt_keySampleLegend.push_back({opName+"_"+WCvalue2, EFT::Sample::linear, "c'_{3} ="+WCvalue2+"Linear"});
  //   vt_keySampleLegend.push_back({opName+"_0.15", EFT::Sample::all, opName+" 0.15 A"});
  // //vt_keySampleLegend.push_back({opName+"_0.15", EFT::Sample::linear, opName+" 0.15 L"});
  //   vt_keySampleLegend.push_back({opName+"_0", EFT::Sample::all, "SM"});
  //
      vt_keySampleLegend.push_back({"data", EFT::Sample::all, "SM"});
  // //
      eft.drawHistogram(vt_keySampleLegend,
                      output_dir+"c"+WCvalue +"/"+ opName+"_" + hName + "", "Fraction", hNameLatex , 0.0001, 1.0, 0.96, 1.03, true,  "", false, "simple");
  //
  // void drawHistogram(const std::vector< std::tuple<std::string, Sample, std::string> > &vt_keySampleLegend,
  //                    const std::string &plotName, const std::string &yLabel, const std::string &xLabel,
  //                    const double &histMin, const double &histMax, const double &ratioMin, const double &ratioMax,
  //                    const bool drawLogY = false, const std::string &legHeader = "",
  //                    const bool divideBinWidth = false, const std::string &ratioMode = "simple");

  // make the range to compute the chi2; in this case [min, max: step]
  //ANDRE  const std::vector<double> v_opPoint = makeInterval(-10., 10., 0.0001);
  //fillInterval(v_opPoint, -3., 3., 0.00001); // alternatively for non-uniform intervals - second arg has to be the same as last vector element

  // make the list of keys to fit
  // in this case op1 and op2 over the specified list of points
  // the points are made gridwise considering all operators
  //ANDRE  eft.listKeyToFit({ {opName, v_opPoint} });

  //ANDRE  const std::vector<EFT::Sample> v_sample = {EFT::Sample::all, EFT::Sample::linear};
  //ANDRE eft.computeFitChi2(v_sample);

  // now we provide the op for which we wanna draw the dChi2 plot on
  // insert into map: op key, op string in plot, op range (if none, select by dChi2), y axis range, x axis range
  //ANDRE std::map<std::string, std::tuple<std::string, std::vector<double>, std::array<double, 2>, std::array<double, 2>>> m_1D_all;
  //m_1D.insert({"c1", { "c1", {/* op range in min, max */}, {0., 9.999}, {-1.499, 1.999} }});
  //const std::string m_1D;
  //m_1D.insert({"c1", { "c1", {/* op range in min, max */}, {0., 9.999}, {-1.499, 1.999} }});


  // in this case we just draw for cQq11 - also include a filename for the resulting plot
  //eft.draw1DChi2(m_1D, output_dir + "c1_constraint", v_sample);
  //ANDRE eft.draw1DChi2({ {m_1D, {opName, {/* op range in min, max */}, {0., 0.999}, {-0.5,0.5} }} }, output_dir + opName+"_constraint", v_sample); //added by Andre

    eft.clearContent();
    std::cout << "EFTFitter: finished histogramming on WC: " << WCvalue << std::endl << std::endl;
  }
    std::cout << "EFTFitter: finished fit on variable " << hName << std::endl << std::endl;
  }

   std::cout << "EFTFitter: finished fit on operator " << opName << std::endl << std::endl;
  }
  // same thing for 2D - but here only the contour is drawn
  // to brace-init the value, provide an extra pair of braces to obey the brace-elision rule
  // insert into map: op pair key (y-x), op1 string in plot, its axis range, op2 string in plot, its axis range
  //std::map<std::array<std::string, 2>, std::array<std::pair<std::string, std::array<double, 2>>, 2>> m_2D;
  //m_2D.insert({ {"c1", "c2"}, {{{"c1", {-0.799, 0.799}}, {"c2", {-0.799, 1.199}}}} });

  //eft.draw2DChi2(m_2D_all, output_dir + "c1_c2_constraint", v_sample);

  return 0;
}
