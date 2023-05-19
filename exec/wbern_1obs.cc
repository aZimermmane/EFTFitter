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
//usage 
//g++ `root-config --cflags --evelibs` -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I ./src -o ./exec/wbern_1obs ./src/EFTFitter.cc ./exec/wbern_1obs.cc
//./exec/wbern_1obs "operator" "dataset" "covariance"
//dataset -> <2016preUL|2016UL>
//covaraince -> <stat|total>

// g++ `root-config --cflags --evelibs` -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -o ./exec/wbern_evo ./exec/wbern_evo.cc ./src/EFTFitter.cc
// for iOp in ut cvv c1 dt cmm cva c3 cav c123 cmp; do ${condorDir}/condorSubmit.sh -s ${condorDir}/condorParam.txt -w ${condorDir}/condorRun.sh -n wbern_evo_${iOp} -l ./wbern_0519/evolution_sm_tmp/${iOp} -e ./exec/wbern_evo -a "${iOp}"; sleep 0.1; done

//const int nRebin = 4;
//const double k_nnlo_lo = 1.0 , br_tt_2l = 0.06635231524673364;

std::array<std::unique_ptr<TGraphAsymmErrors>, 2> fitResult(const std::string &op = "", const std::string &fileName = "",
                                                            const int &iIt = 0, const int &iVar = 0,
                                                            const bool useAll = true) {
  using TG = TGraphAsymmErrors;
  std::array<std::unique_ptr<TG>, 2> a_graph = {nullptr, nullptr};
  if (op == "" or fileName == "")
    return a_graph;

  const std::string samp = (useAll) ? "all" : "linear", tag = toStr(iIt) + "_" + toStr(iVar);
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));

  if (!file->GetListOfKeys()->Contains((op + "_sigma2_" + samp).c_str()))
    return a_graph;

  // read the real best fit and sigmas
  a_graph.at(0) = std::unique_ptr<TG>(dynamic_cast<TG *>(( file->Get( (op + "_sigma1_" + samp).c_str() ))->Clone( (tag + "_s1").c_str() ) ));
  a_graph.at(1) = std::unique_ptr<TG>(dynamic_cast<TG *>(( file->Get( (op + "_sigma2_" + samp).c_str() ))->Clone( (tag + "_s2").c_str() ) ));

  return a_graph;
}

int main(int argc, char** argv) {
  if (argc < 4){
      std::cout << "EFTFitter ERROR: insuficient number of positional arguments " << std::endl;
    return 0;}
  
  const std::string opName = argv[1], dataset=argv[2], covariance=argv[3], covMatrixTag=argv[4];
   
  if (dataset != "2016preUL" && argc < 5){
      std::cout << "EFTFitter ERROR: one positional arguments missing. If a UL dataset, remember to provide covMatrixTag " << std::endl;
    return 0;}

  using EFT = EFTFitter;

  // common flags for an evolution test
  // fit with data or SM, check evolution of all or linear part, 2 sigma or 1 sigma interval width as figure of merit
  const bool useData = false, useAll = false;
  // construction done with the key we want to treat as data, lambda, fit mode and stat mode (optionally sum of shape templates)
  //EFT eft("data", 1., EFT::Fit::shape, EFT::Stat::xsec); // Modified by andre 01.02.23
  EFT eft("data", 1., EFT::Fit::hybrid, EFT::Stat::xsec);
  //const std::string inDir = "/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/wbern_0519/root/", opName = argv[1];
    std::string inDir;
  //const std::string inDir = "../histograms/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/", opName = argv[1];
//  const std::string covMatrix ="Systematics_AllVars_1D_1000PE_MCStatFixed/TotalStatCovMatrix_AllVarNorm_rebinnedA";// Modified by andre 01.02.23
  //const std::string covMatrix ="Systematics_AllVars_1D_1000PE_MCStatFixed/TotalStatSystCovMatrix_AllVarNorm_rebinnedA";// Modified by andre 01.02.23
    std::string covMatrixPath, covMatrixFile, covMatrixBranch, covMatrixBranch1,covMatrixBranch2;
  if (dataset == "2016preUL"){
    inDir = "/afs/desy.de/user/z/zimermma/work/EFTFitter/wbern_0314/root/pre_UL_";
    covMatrixPath ="/afs/desy.de/user/z/zimermma/work/EFTFitter/inputs/covariance_matrix/";
    covMatrixFile = "covmat_190114";
    if (covariance == "stat"){
      covMatrixBranch = "TotalStat_shape_a";
     // covMatrixTag = "";
    }
    else if (covariance == "total"){
      covMatrixBranch = "TotalStatSyst_shape_a";
     // covMatrixTag = "";
    }
    else{
      std::cout << "EFTFitter ERROR: covariance positional argument (3) must be either \"stat\" or \"total\" " << std::endl;
      return 0;
    }
  }
  else if (dataset == "2016UL"){
      inDir = "/afs/desy.de/user/z/zimermma/work/EFTFitter/wbern_0314/root/";
      covMatrixPath ="/afs/desy.de/user/z/zimermma/work/EFTFitter/inputs/covariance_matrix/";
      covMatrixFile = "Systematics_AllVars_1D_1000PE_MCStatFixed";
    if (covariance == "stat"){
      covMatrixBranch = "TotalStatCovMatrix_AllVarNorm_rebinnedA";
   //   covMatrixTag = "_1000PE_MCStatFixed";
    }
    else if (covariance == "total"){
      covMatrixBranch1 = "TotalStatCovMatrix_AllVarNorm_rebinnedA";
      covMatrixBranch2 = "TotalSystCovMatrix_AllVarNorm_rebinnedA";
     // covMatrixTag = "_1000PE_MCStatFixed";
    }
    else{
      std::cout << "EFTFitter ERROR: covariance positional argument (3) must be either \"stat\" or \"total\" " << std::endl;
      return 0;
    }
  }
  else{
    std::cout << "EFTFitter ERROR: dataset positional argument (2) must be either \"2016preUL\" or \"2016UL\" " << std::endl;
    return 0;
  }

  const std::string outDir = "/afs/desy.de/user/z/zimermma/work/EFTFitter/output-2016-1D/sigleObsFit/" + dataset + "_" + covariance + covMatrixTag + "/" + opName + "/";

  std::cout << "EFTFitter debuging: opName: " << opName << std::endl;
  std::cout << "EFTFitter debuging: dataset: " << dataset << std::endl;
  std::cout << "EFTFitter debuging: covariance: " << covariance << std::endl;
  std::cout << "EFTFitter debuging: covariance TAG: " << covMatrixTag << std::endl;


  gSystem->Exec(("mkdir -p "+ outDir).c_str() );

  // auto mkdir = std::experimental::filesystem::create_directories(outDir);
  // if (mkdir) {
  //       std::cout << "created output directory: " << std::endl;
  //       std::system("tree x_tmp");
  //   } else {
  //       std::cout << "create_directories() failed" << std::endl;
  //   }

  // goes forever due to 2D degeneracy: {cvv, c1} {cva, c3}, {cav, c123}
  const std::map<std::string, double> m_range = {{"ctG",0.2},{"ut", 0.5}, {"cvv", 1.0}, {"c1", 15.0},
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
  const int binToIgnore = 0;
  //const int binWidth = 3.; // adding the bin width in the theory prediction


  // make the range to interpolate over; in this case [min, max: step]
  const double eps = 0.0001;
  //const double eps = 0.00001;
  const std::vector<double> v_opPoint = makeInterval(-opRange, opRange, opRange / 10000.);
  const std::vector<EFT::Sample> v_sample = {EFT::Sample::linear};
  const std::string sample = (useAll) ? "all" : "linear";

  const std::vector<std::string> v_hStr = {"b1k", "b2k", "b1r", "b2r", "b1n", "b2n", "b1j", "b2j", "b1q", "b2q",
                                            "ckk", "crr", "cnn",
                                            "cP_rk", "cM_rk", "cP_nr", "cM_nr", "cP_nk", "cM_nk",
                                            "ckj","crq",
                                            "chan","ctra","csca",
                                            "cP_rj","cM_rj",
                                            "c_kjL", "c_rqL",
                                            "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM",
                                            "cHel"};

  //const std::vector<std::string> v_hStr = {"cHel"};
                                            //need to add an input file that contains b1j
  const std::vector<int> v_all_var = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};

  //const std::vector<int> v_all_var = {34};

  // just in case of typos
  if (v_hStr.size() != v_all_var.size())
    throw std::range_error( "EFTFitter: List of variable names doesn't match list of variable index. Something is likely wrong." );

  // dummy rate as usual
  const std::array<double, 2> rate_zero = {0., 0.};

  std::vector<std::unique_ptr<TF1>> v_spinvar;
  v_spinvar.emplace_back( std::make_unique<TF1>("f_bli", "0.5 * (1. + ([0] * x))", -1., 1.) );
  v_spinvar.emplace_back( std::make_unique<TF1>("f_cii", "0.5 * (1. - ([0] * x)) * std::log(1. / std::abs(x))", -1., 1.) );
  v_spinvar.emplace_back( std::make_unique<TF1>("f_cPMij", "0.5 * (1. - (0.5 * [0] * x)) * std::acos(std::abs(x))", -1., 1.) );
  v_spinvar.emplace_back( std::make_unique<TF1>("f_cHel", "0.5 * (1. - ([0] * x))", -1., 1.) );

  // variable index
  auto evoFile = std::make_unique<TFile>((outDir + opName + "_dchi2.root").c_str(), "recreate");
  std::vector<int> v_fit_var;

    std::cout << "EFTFitter: starting..." << std::endl;
    std::vector<std::array<std::unique_ptr<TGraphAsymmErrors>, 2>> v_current_result;


    for (int iFit = 0; iFit < v_hStr.size(); ++iFit) {
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

        // const std::string histName = v_hStr.at(iFit);
        // eft.addRawInput(opName + "_0", EFT::Sample::all, inDir + "ttbareft_translationp_Breuther_dim6top_reweighting_13TeV_SM_hist.root",
        //                 "gen_" + histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);
        //
        // // 1D inputs - here c1 = 10 and c1 = 272 is chosen as raw inputs
        // eft.addRawInput(opName + "_-4", EFT::Sample::all, inDir + "ttbareft_translationp_Breuther_dim6top_reweighting_13TeV_" + opName +"_m4_hist.root",
        //                 "gen_" + histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);
        //
        // eft.addRawInput(opName + "_8", EFT::Sample::all, inDir + "ttbareft_translationp_Breuther_dim6top_reweighting_13TeV_" + opName +"_8_hist.root",
        //                 "gen_" + histName, "", nRebin, {(k_nnlo_lo / br_tt_2l), 0.}, EFT::Stat::xsec);

        // emulate coeff fits with arbitrary template selection
        // be very specific about the captures!
        auto f_emulate = [&v_iEB, &rate_zero, &binToIgnore, &v_spinvar, &dataset] (const auto &v_binC) {
          // convert the coeffs into shapes using the known parametrization
          std::vector<std::array<double, 2>> v_tmpC = {rate_zero, rate_zero};
          TF1 *f_ref = v_spinvar.at(0).get();

          for (int iC = 2; iC < v_binC.size(); ++iC) {
            if (!std::count(std::begin(v_iEB), std::end(v_iEB), iC - 2)) continue;
              if (dataset == "2016preUL"){
                // get the function ref
                if (iC - 2 < 10)
                  f_ref = v_spinvar.at(0).get();
                else if (iC - 2 < 13)
                  f_ref = v_spinvar.at(1).get();
                else if (iC - 2 < 19)
                  f_ref = v_spinvar.at(2).get();
                else if (iC - 2 < 20)
                    f_ref = v_spinvar.at(3).get();
              }
              else
             // get the function ref
            if (iC - 2 < 10)
              f_ref = v_spinvar.at(0).get();
            else if (iC - 2 < 13)
              f_ref = v_spinvar.at(1).get();
            else if (iC - 2 < 19)
              f_ref = v_spinvar.at(2).get();
            else if (iC - 2 < 21)
                f_ref = v_spinvar.at(1).get();
            else if (iC - 2 < 24)
              f_ref = v_spinvar.at(3).get();
            else if (iC - 2 < 35)
                f_ref = v_spinvar.at(3).get();

            //std::cout << "DEBUG: Using this fit function: " << f_ref->GetName() << std::endl;


            const bool is_cii = (std::string(f_ref->GetName()).find("f_cii") != std::string::npos);

            f_ref->SetParameter(0, v_binC.at(0).at(0) * v_binC.at(iC).at(0));
            const double integral = (is_cii) ? f_ref->Integral(-1., -DBL_MIN) + f_ref->Integral(DBL_MIN, 1.) : f_ref->Integral(-1., 1.);
            //for (int iB = 0; iB < 6; ++iB) {
            for (int iB = 0; iB < 6; ++iB) {
              if (iB == binToIgnore) continue;

              double min = -1. + (iB / 3.), max = -1. + ((iB + 1) / 3.);
              if (is_cii) {
                if (min == 0.)
                  min = DBL_MIN;
                if (max == 0.)
                  max = -DBL_MIN;
              }

              //v_tmpC.push_back({0.33*f_ref->Integral(min, max) / integral, 0.}); // here the 0.33 factor divides by the bin width
              //v_tmpC.push_back({3.*f_ref->Integral(min, max) / (integral), 0.}); // here the 0.33 factor divides by the bin width
              v_tmpC.push_back({f_ref->Integral(min, max) / integral, 0.}); // here the 0.33 factor divides by the bin width

            }
          }

          return v_tmpC;
        };
        eft.setHybridTransformation(std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> (f_emulate));

        eft.addRawInput(opName + "_0", EFT::Sample::all,
                        inDir + opName + "_coeff.root", "sm_nominal", "", 1, rate_zero, EFT::Stat::xsec);
        eft.addRawInput(opName + "_1", EFT::Sample::all,
                        inDir + opName + "_coeff.root", "sm_nominal_plus_" + opName + "_nominal", "", 1, rate_zero, EFT::Stat::xsec);
        eft.addRawInput(opName + "_-1", EFT::Sample::all, inDir + opName + "_coeff.root",
                      "sm_nominal_minus_" + opName + "_nominal", "", 1, rate_zero, EFT::Stat::xsec);

        // prepare the base for interpolation
        eft.prepareInterpolationBase();

        // data
        if (useData) {
          auto f_shape_data = [&v_iEB, &rate_zero, &binToIgnore] (const auto &v_binC) {
            // list of variables to include - ensure consistent with main fitter loop below
            const double shapeSum = 22.;

            const int nBin = v_binC.size();
            const int nBinEach = (nBin - 2) / int(shapeSum);

            // actual rate to be used in addHybridData()
            std::vector<std::array<double, 2>> v_tmpC = {{v_binC.at(0).at(0), v_binC.at(0).at(1)}, rate_zero};

            for (int iB = 2; iB < v_binC.size(); ++iB) {
              const int iIt = (iB - 2) / 6;
              if (!std::count(std::begin(v_iEB), std::end(v_iEB), iIt)) continue;
              if ((iB - 2) % nBinEach == binToIgnore) continue;

              // no shapeSum scaling unlike in MC since addHybridData doesn't normalize the binContent
              v_tmpC.push_back({v_binC.at(iB).at(0), v_binC.at(iB).at(1)});
            }

            return v_tmpC;
          };
          eft.addHybridData("/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/covariance_matrix/unfolded_data_190114.root", "TTbarSpinDensityMatrix/snake_spinCorr_shape_data_a",
                            rate_zero, EFT::Stat::xsec, f_shape_data);
        }
        // else
          eft.assignAsData(opName + "_0", EFT::Sample::linear, false);


        // grab the cov matrix
        std::vector<std::array<int, 2>> covMat_binRange;
        for (auto &var : v_iEB){
          // if (dataset == "preUL" && var > 33){
          // std::cout << "EFTFitter preUL: var reached: "  << var << std::endl;
          //   var=var-15;
          //  std::cout << "New assigment var: "  << var << std::endl;
          //   covMat_binRange.push_back({(2+(var) * 6), ((var + 1) * 6)}); //modified by Andre 01.02. (always ignore the first bin of each dib)
          //   break;
          // }  
         //changes the correct indice for cHel in preUL cov mat
          covMat_binRange.push_back({(2+(var) * 6), ((var + 1) * 6)}); //modified by Andre 01.02. (always ignore the first bin of each dib)

          //covMat_binRange.push_back({((var + 1) * 6) - 5, ((var + 1) * 6 -1)}); //modified by Andre 01.02. (always ignore the last bin of each dib)
          //covMat_binRange.push_back({((var + 1) * 5) - 4, ((var + 1) * 5)});
        };
        //std::cout << "EFTFitter: covMat_binRange debugging" << covMat_binRange << std::endl;
        //eft.readCovMatRoot("finalcov", "/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/covariance_matrix/covmat_190114.root", covMatrix, covMat_binRange);
        if (dataset != "2016preUL"){
          if (covariance == "stat")
            eft.readCovMatRoot("finalcov", covMatrixPath + covMatrixFile + ".root", covMatrixBranch, covMat_binRange);
          else if (covariance == "total"){
            eft.readCovMatRoot("totalStat", covMatrixPath + covMatrixFile + ".root", covMatrixBranch1, covMat_binRange);
            eft.readCovMatRoot("totalSyst", covMatrixPath + covMatrixFile + ".root", covMatrixBranch2, covMat_binRange);
            eft.makeFinalCovMat({"totalStat","totalSyst"});
            std::cout << "covMat_binRange :"  << covMat_binRange[0] << std::endl;
            std::cout << "covMat_path :"  << covMatrixPath + covMatrixFile + ".root"<< std::endl;
            std::cout << "covMatrixBranch1 :"  << covMatrixBranch1 << std::endl;
            std::cout << "covMatrixBranch2 :"  << covMatrixBranch2 << std::endl;   
       }
          else
            std::cout << "EFTFitter ERROR: covariance positional argument (2) must be either \"stat\" or \"total\" " << std::endl;
        }
        else if (dataset == "2016preUL"){
        std::cout << "covMat_binRange :"  << covMat_binRange[0] << std::endl;
        std::cout << "covMat_path :"  << covMatrixPath + covMatrixFile + ".root"<< std::endl;
        eft.readCovMatRoot("finalcov", covMatrixPath + covMatrixFile + ".root", covMatrixBranch, covMat_binRange);
        }


        //eft.readCovMatRoot("totalStat", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_132x132.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
        //eft.readCovMatRoot("finalcov", "/afs/desy.de/user/z/zimermma/work/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_228x228_1000PE.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
        //eft.readCovMatRoot("totalStat", "/afs/desy.de/user/z/zimermma/work/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_1000PE_MCStatFixed.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
        //eft.readCovMatRoot("totalSyst", "/afs/desy.de/user/z/zimermma/work/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_1000PE_MCStatFixed.root", "TotalSystCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
        //eft.readCovMatRoot("finalcov", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/covMatSyst_purdue_full2016.root", covMatrix, covMat_binRange);
        //eft.makeFinalCovMat({"totalStat","totalSyst"});

        eft.drawCovMat(outDir,{},true); //if true drwas correlation matrix instead
        //gSystem->Exec( ("mv " + outDir + "cov_finalcov.txt " + outDir + "covMat_" + v_hStr.at(iFit) + "_iter_" + toStr(iIt) + ".txt").c_str() );
        //gSystem->Exec( ("mv " + outDir + "cov_finalcov.pdf " + outDir + "covMat_" + v_hStr.at(iFit) + "_iter_" + toStr(iIt) + ".pdf").c_str() );
        //gSystem->Exec( ("rm " + outDir + "cov_all.root").c_str() );

        std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
        vt_keySampleLegend.push_back({opName + "_-1", EFT::Sample::all, opName + " -1"});
        vt_keySampleLegend.push_back({opName + "_1", EFT::Sample::all, opName + " 1"});

        vt_keySampleLegend.push_back({opName + "_0", EFT::Sample::all, "SM"});
        vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});

        eft.drawHistogram(vt_keySampleLegend,
                            outDir + opName + "_var_" + v_hStr.at(iFit) + "_shape",
                            "Fraction", "Index", -0.4999, 0.4999, 0.0001, 1.9999,
                            false, "", false, "none");

        eft.listKeyToFit({ {opName, v_opPoint}});
        eft.computeFitChi2(v_sample);

        eft.draw1DChi2({ {opName, {opName, {/* op range in min, max */}, {0., 9.999}, {-opRange + eps, opRange - eps} }} }, outDir, v_sample);

        eft.clearContent();

        v_current_result.emplace_back(fitResult(opName, outDir + opName + "_dChi2.root", 0, iFit, useAll));
        if (v_current_result.back().at(0) == nullptr)
          std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " doesn't result in a constraint!" << std::endl;
        else
          std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " saved." << std::endl;

        gSystem->Exec( ("rm " + outDir + opName + "_dChi2.root").c_str() );
  
    }

    evoFile->cd();
    for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
      int iVar = iRes;
      if (dataset == "2016preUL" && iRes> 18 && iRes < 34)
      iVar = iVar + 15;
      if (v_current_result.at(iRes).at(0) == nullptr) continue;

        v_current_result.at(iRes).at(0)->SetName( (sample + "_" + toStr(v_hStr.at(iVar)) + "_sigma1").c_str() );
        v_current_result.at(iRes).at(0)->Write( (sample + "_" + toStr(v_hStr.at(iVar)) + "_sigma1").c_str() );

        v_current_result.at(iRes).at(1)->SetName( (sample + "_" + toStr(v_hStr.at(iVar)) + "_sigma2").c_str() );
        v_current_result.at(iRes).at(1)->Write( (sample + "_" + toStr(v_hStr.at(iVar)) + "_sigma2").c_str() );

    }

    std::cout << "EFTFitter: all obsevables fit completed!" << std::endl;

  return 0;
}
