// example code using the EFTFitter plugin
// plugins rely on ROOT6 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh exec/simple_fit.cc
// exec config expanding the coeffs into shapes and fitting
// version to perform the evolution plot


#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"

#include "../src/EFTFitter.h"

//g++ `root-config --cflags --evelibs` -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I ./src -o ./exec/spin_mtt_to_1D.cc ./src/EFTFitter.cc ./exec/spin_mtt_to_1D.cc

// g++ `root-config --cflags --evelibs` -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -o ./exec/wbern_evo ./exec/wbern_evo.cc ./src/EFTFitter.cc
// for iOp in ut cvv c1 dt cmm cva c3 cav c123 cmp; do ${condorDir}/condorSubmit.sh -s ${condorDir}/condorParam.txt -w ${condorDir}/condorRun.sh -n wbern_evo_${iOp} -l ./wbern_0519/evolution_sm_tmp/${iOp} -e ./exec/wbern_evo -a "${iOp}"; sleep 0.1; done

//const int nRebin = 4;
//const double k_nnlo_lo = 1.0 , br_tt_2l = 0.06635231524673364;

std::array<std::unique_ptr<TGraphAsymmErrors>, 2> fitResult(const std::string &op = "", const std::string &fileName = "",
                                                            const int &iBin = 0, const int &iVar = 0,
                                                            const bool useAll = true) {
  using TG = TGraphAsymmErrors;
  std::array<std::unique_ptr<TG>, 2> a_graph = {nullptr, nullptr};
  if (op == "" or fileName == "")
    return a_graph;

  const std::string samp = (useAll) ? "all" : "linear", tag = toStr(iBin) + "_" + toStr(iVar);
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));

  if (!file->GetListOfKeys()->Contains((op + "_sigma2_" + samp).c_str()))
    return a_graph;

  // read the real best fit and sigmas
  a_graph.at(0) = std::unique_ptr<TG>(dynamic_cast<TG *>(( file->Get( (op + "_sigma1_" + samp).c_str() ))->Clone( (tag + "_s1").c_str() ) ));
  a_graph.at(1) = std::unique_ptr<TG>(dynamic_cast<TG *>(( file->Get( (op + "_sigma2_" + samp).c_str() ))->Clone( (tag + "_s2").c_str() ) ));

  return a_graph;
}

int main() {
//int main(int argc, char** argv) {
  // if (argc < 3)
  //   return 0;

  using EFT = EFTFitter;

  // common flags for an evolution test
  // fit with data or SM, check evolution of all or linear part, 2 sigma or 1 sigma interval width as figure of merit
  const bool useAll = false;
  //const bool useData = false;

  // construction done with the key we want to treat as data, lambda, fit mode and stat mode (optionally sum of shape templates)
  //EFT eft("data", 1., EFT::Fit::shape, EFT::Stat::xsec); // Modified by andre 01.02.23
  EFT eft("data", 1., EFT::Fit::shape, EFT::Stat::xsec,1.);
  //const std::string inDir = "/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/wbern_0519/root/", opName = argv[1];
  const std::string dataset = "2016UL";
  const std::string input_dir = "../analyseRoot/hists/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/";
  const std::string mcSample = "hist_Nevents10000000";
  //const std::string opName = argv[1], n_2DBins =argv[2];
  const std::string covMatrix_inDir = "/afs/desy.de/user/z/zimermma/work/EFTFitter/inputs/covariance_matrix/";
  const std::string covMatrix ="Systematics_AllVars_2D_100PE_MCStatFixed";// Modified by andre 01.02.23
  const std::string opName = "ctG";
  const int nCovMatBins = 6*4;
  const double eps = 0.0001;
  const int nRebin = 1;
  const int binToIgnore = 0;
  const int n_2DBins = 4;
  //const std::string covMatrix ="TotalStatCovMatrix_AllVarNorm_rebinnedA";
  //const std::string covMatrix ="TotalSystCovMatrix_AllVarNorm_rebinnedA";
  //const std::string covMatrix ="TotalSystStatCovMatrix_AllVarNorm_rebinnedA";
  //const std::string covMatrix = "TotalStat_shape_a_drop_bin_1";
  const std::string outDir = "/afs/desy.de/user/z/zimermma/work/EFTFitter/spin_mtt_to1D/" + mcSample + "_" + covMatrix + "/" + opName + "/";

  gSystem->Exec(("mkdir -p "+ outDir).c_str() );

  const std::map<std::string, double> m_range = {{"ctG",1.0},{"ut", 0.5}, {"cvv", 1.0}, {"c1", 15.0},
                                                 {"dt", 0.2}, {"cmm", 0.5},
                                                 {"cva", 1.0}, {"c3", 3.0},
                                                 {"cav", 1.0}, {"c123", 3.0},
                                                 {"cmp", 1.0},
                                                 {"b1k",2.0}, {"b2k",2.0}, {"b1r",2.0}, {"b2r",2.0},
                                                 {"b1n",2.0}, {"b2n",2.0}, {"b1j",2.0}, {"b2j",2.0},
                                                 {"b1q",2.0}, {"b2q",2.0}, {"ckk",2.0}, {"crr",2.0},
                                                 {"cnn",2.0}, {"cPrk",2.0}, {"cMrk",2.0}, {"cPnr",2.0},
                                                 {"cMnr",2.0}, {"cPnk",2.0}, {"cMnk",2.0}, {"ckj",2.0},
                                                 {"crq",2.0}, {"cPrj",2.0}, {"cMrj",2.0}, {"chel",1.0}};

  const double opRange = m_range.at(opName);




  // make the range to interpolate over; in this case [min, max: step]

  //const double eps = 0.00001;
  const std::vector<double> v_opPoint = makeInterval(-opRange, opRange, opRange / 10000.);
  const std::vector<EFT::Sample> v_sample = {EFT::Sample::linear};
  const std::string sample = (useAll) ? "all" : "linear";

  // const std::vector<std::string> v_hStr = {"b1k", "b2k", "b1r", "b2r", "b1n", "b2n", "b1kStar", "b2kStar", "b1rStar", "b2rStar",
  //                                          "ckk", "crr", "cnn",
  //                                          "cP_rk", "cM_rk", "cP_nr", "cM_nr", "cP_nk", "cM_nk",
  //                                          "ckkStar","crrStar",
  //                                          "chan","ctra","csca",
  //                                          "cP_rkStar","cM_rkStar",
  //                                          "c_kkStarL", "c_rrStarL",
  //                                          "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM",
  //                                          "cHel"};

  const std::vector<std::string> v_hStr = {"cHel"};
                                           //need to add an input file that contains b1j
  //const std::vector<int> v_all_var = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
  const std::vector<int> v_all_var = {34};

  // just in case of typos
  if (v_hStr.size() != v_all_var.size())
    throw std::range_error( "EFTFitter: List of variable names doesn't match list of variable index. Something is likely wrong." );

  // dummy rate as usual
  const std::array<double, 2> rate_zero = {0., 0.};

  // variable index
  auto evoFile = std::make_unique<TFile>((outDir + opName + "_dchi2.root").c_str(), "recreate");
  std::vector<int> v_fit_var;



  std::string histName;
  //std::string iBin;
  for (int iFit = 0; iFit < v_hStr.size(); ++iFit) 
  {
  int iVar = v_all_var.at(iFit);
  if (dataset == "2016preUL" && iVar > 18 && iVar < 34)
  {std::cout << "EFTFitter preUL: Excluding variable " << v_hStr.at(iFit) << " from the fit" << std::endl;
    continue;}
  else if (dataset == "2016preUL" && iVar == 34){
    std::cout << "EFTFitter preUL: we reached: "  << v_hStr.at(iFit) << std::endl;
    iVar = iVar - 15;
  }
  std::vector<int> v_iEB(v_fit_var);
  v_iEB.push_back( iVar );
//     std::cout << "EFTFitter: v_iEB debugging" << v_iEB << std::endl;
  std::sort(std::begin(v_iEB), std::end(v_iEB));

  std::cout << "EFTFitter: starting fit on variable " << v_hStr.at(iFit) << std::endl;
  std::cout << "DEBUG: n_2DBins: " << n_2DBins << std::endl;
  
  for (int iBin = 1; iBin < n_2DBins+1; ++iBin) {
  //for (int iBin = 0; iBin < 1; ++iBin) {
    std::cout << "EFTFitter: bin " << iBin << " starting..." << std::endl;
    std::vector<std::array<std::unique_ptr<TGraphAsymmErrors>, 2>> v_current_result;


    histName = "gen_" + v_hStr.at(iFit)+"_mtt"+"_bin" + iBin;

    std::cout << "DEBUG:eft.addRawInput: " << "input sample:" << input_dir + mcSample+"_baseline.root"  << std::endl;
    std::cout << "DEBUG:                 " << "Hist:" <<  histName+"_baseline"  << std::endl;
      
    eft.addRawInput(opName+"_0", EFT::Sample::all, input_dir + mcSample+"_baseline.root",
                         histName+"_baseline", "", nRebin, rate_zero, EFT::Stat::xsec);
    eft.addRawInput(opName+"_-4", EFT::Sample::all, input_dir + mcSample+"_"+opName+"_m4.root",
                                  histName+"_"+opName+"_m4", "", nRebin,rate_zero, EFT::Stat::xsec);
    eft.addRawInput(opName+"_8", EFT::Sample::all, input_dir + mcSample+"_"+opName+"_8.root",
                  histName+"_"+opName+"_8", "", nRebin, rate_zero, EFT::Stat::xsec);


    // prepare the base for interpolation
    eft.prepareInterpolationBase();

    // data
    // if (useData) {
    //   auto f_shape_data = [&v_iEB, &rate_zero, &binToIgnore] (const auto &v_binC) {
    //     // list of variables to include - ensure consistent with main fitter loop below
    //     const double shapeSum = 22.;

    //     const int nBin = v_binC.size();
    //     const int nBinEach = (nBin - 2) / int(shapeSum);

    //     // actual rate to be used in addHybridData()
    //     std::vector<std::array<double, 2>> v_tmpC = {{v_binC.at(0).at(0), v_binC.at(0).at(1)}, rate_zero};

    //     for (int iB = 2; iB < v_binC.size(); ++iB) {
    //       const int iBin = (iB - 2) / 6;
    //       if (!std::count(std::begin(v_iEB), std::end(v_iEB), iBin)) continue;
    //       if ((iB - 2) % nBinEach == binToIgnore) continue;

    //       // no shapeSum scaling unlike in MC since addHybridData doesn't normalize the binContent
    //       v_tmpC.push_back({v_binC.at(iB).at(0), v_binC.at(iB).at(1)});
    //     }

    //     return v_tmpC;
    //   };
    //   eft.addHybridData("/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/covariance_matrix/unfolded_data_190114.root", "TTbarSpinDensityMatrix/snake_spinCorr_shape_data_a",
    //                     rate_zero, EFT::Stat::xsec, f_shape_data);
    // }
    // else
    eft.assignAsData(opName + "_0", EFT::Sample::linear, false);

    // grab the cov matrix
    std::vector<std::array<int, 2>> covMat_binRange;
    for (auto &var : v_iEB)

    { covMat_binRange.push_back({(1+ (var) *nCovMatBins + (iBin-1)*6), ((var) *nCovMatBins + (iBin)*6)}); //modified by Andre 01.02. (always ignore the first bin of each dib)
      //covMat_binRange.push_back({((var + 1) * 6) - 5, ((var + 1) * 6 -1)}); //modified by Andre 01.02. (always ignore the last bin of each dib)

      std::cout << "DEBUG: var: " << var << std::endl;
      std::cout << "DEBUG: covMat_binRange: " << (1+ (var) *nCovMatBins + (iBin-1)*6)<< "," << ((var) *nCovMatBins + (iBin)*6) << std::endl;
      //covMat_binRange.push_back({((var + 1) * 5) - 4, ((var + 1) * 5)});
      };
    eft.readCovMatRoot("finalcov", covMatrix_inDir + covMatrix +".root", "TotalStatCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);

   
    eft.drawCovMat(outDir,{},true);
    gSystem->Exec( ("mv " + outDir + "cov_finalcov.txt " + outDir + "covMat_" + v_hStr.at(iFit) + "_iter_" + toStr(iBin) + ".txt").c_str() );
    gSystem->Exec( ("mv " + outDir + "cov_finalcov.pdf " + outDir + "covMat_" + v_hStr.at(iFit) + "_iter_" + toStr(iBin) + ".pdf").c_str() );
    gSystem->Exec( ("rm " + outDir + "cov_all.root").c_str() );
    std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
    vt_keySampleLegend.push_back({opName + "_-1", EFT::Sample::linear, opName + " -1"});
    vt_keySampleLegend.push_back({opName + "_1", EFT::Sample::linear, opName + " 1"});

    vt_keySampleLegend.push_back({opName + "_0", EFT::Sample::linear, "SM"});
    vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});

    eft.drawHistogram(vt_keySampleLegend,
                        outDir + opName + "_var_" + v_hStr.at(iFit) + "_bin_" + toStr(iBin) + "_shape",
                        "Fraction", "Index", -0.4999, 0.4999, 0.0001, 1.9999,
                        false, "", false, "none");

    eft.listKeyToFit({ {opName, v_opPoint}});
    eft.computeFitChi2(v_sample,binToIgnore);

    eft.draw1DChi2({ {opName, {opName, {/* op range in min, max */}, {0., 9.999}, {-opRange + eps, opRange - eps} }} }, outDir, v_sample);

    eft.clearContent();

    v_current_result.emplace_back(fitResult(opName, outDir + opName + "_dChi2.root", iBin, iFit, useAll));
    if (v_current_result.back().at(0) == nullptr)
      std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " in bin" << iBin << " doesn't result in a constraint!" << std::endl;
    else
      std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " in bin" << iBin << " saved." << std::endl;

    gSystem->Exec( ("rm " + outDir + opName + "_dChi2.root").c_str() );


    // // first find the best variable in a given iter by minimum sig2 or sig1 width
    int iMinF = -999, iMinV = -999;
    double width;
    int iRes = 0;
    // for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
        if (v_current_result.at(iRes).at(0) == nullptr) continue;
        bool useSig2 = false;
        auto &graph = (useSig2) ? v_current_result.at(iRes).at(1) : v_current_result.at(iRes).at(0);

    //    if (graph->GetErrorXlow(0) + graph->GetErrorXhigh(0) < width) {
         width = graph->GetErrorXlow(0) + graph->GetErrorXhigh(0);
         iMinF = iRes;
         iMinV = v_all_var.at(iRes);
    //    }
    //  }

    // if (iMinF == -999) {
    //   std::cout << "EFTFitter: bin" << iBin << " aborted; no further variables can make up an independent set." << std::endl;
    //   break;
    // }
    //if (v_current_result.at(iRes).at(0) == nullptr) continue;
    v_fit_var.push_back(iMinV);
    std::cout << "EFTFitter: variable " << v_hStr.at(iMinF) << " with 1sigma width " << width << " is saved in bin" << iBin << std::endl;

    evoFile->cd();
    for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
      if (v_current_result.at(iRes).at(0) == nullptr) continue;
      v_current_result.at(iRes).at(0)->SetName( ("bin_" + toStr(iBin) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma1").c_str() );
      v_current_result.at(iRes).at(0)->Write( ("bin_" + toStr(iBin) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma1").c_str() );
      v_current_result.at(iRes).at(1)->SetName( ("bin_" + toStr(iBin) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma2").c_str() );
      v_current_result.at(iRes).at(1)->Write( ("bin_" + toStr(iBin) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma2").c_str() );  
    }
    std::cout << "EFTFitter: bin" << iBin << " completed!" << std::endl;
    }
  }
  return 0;
}
