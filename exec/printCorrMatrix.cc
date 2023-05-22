// -*- C++ -*-
// example code using the EFTFitter plugin
// plugins rely on ROOT6 and C++14 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh ./exec/printCorrMatrix.cc
// version for doing regular fits ie absolute or shape

#include "../src/EFTFitter.h"

int main() {
  using EFT = EFTFitter;
  // construction done with the key we want to treat as data, lambda, fit mode and stat mode
  EFT eft("data", 1., EFT::Fit::hybrid, EFT::Stat::xsec,3.);

  // these are just to avoid writing them out repeatedly
  const std::string covMat_file = "./inputs/covariance_matrix/Systematics_AllVars_1D_228x228_1000PE.root";
  const std::string output_dir = "./fit_output";
  // const std::string histName = "gen_cLab";
  // const std::string sumName = ""; //Modified by Andre 02.02..23
  //     //const std::string sumName = "TTbarSpinDensityMatrix/sumWgt_noCut";
  // const int histIndex = 20; // index of Hist in the order it appears on the covMatt, ranging from 1 to 22
  // const int nRebin = 1;
  // const int nCovMatBins = 6;


  // grab the total stat error matrix - as usual matrix name is such that file->Get("some_matrix_name") works
  // can also partially extract along the diagonal, pass a vector of bin index range eg {{1, 6}, {115, 120}} as last arg
  //eft.readCovMatRoot("finalcov", covMat_file, covMat ,{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
  //eft.readCovMatRoot("finalcov", input_dir + "covariance_matrix/covMatSyst_purdue_full2016.root", "TotalSystCovMatrix_AllVarNorm_rebinnedA",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
  //eft.readCovMatRoot("finalcov", "./inputs/covariance_matrix/covmat_190114.root", "TotalStat_shape_a",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
  // eft.readCovMatRoot("finalcov", "./inputs/covariance_matrix/covmat_190114.root", "TotalSyst_shape_a",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
  // eft.readCovMatRoot("totalStat", "./inputs/covariance_matrix/covMatStat_purdue_full2016.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
  // eft.readCovMatRoot("totalSyst", "./inputs/covariance_matrix/covMatSyst_purdue_full2016.root", "TotalSystCovMatrix_AllVarNorm_rebinnedA",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
  //eft.readCovMatRoot("totalStat", "./inputs/covariance_matrix/Systematics_AllVars_1D_132x132.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});
//  eft.readCovMatRoot("totalSyst", "./inputs/covariance_matrix/Systematics_AllVars_1D_132x132.root", "TotalSystCovMatrix_AllVarNorm_rebinnedA",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});

  //eft.readCovMatRoot("finalcov", "./inputs/covariance_matrix/Systematics_AllVars_1D_132x132.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA",{{(histIndex-1)*nCovMatBins+1, (histIndex-1)*nCovMatBins+(nCovMatBins)}});

  //eft.makeFinalCovMat({"totalStat","totalSyst"});


  eft.readCovMatRoot("finalcov", covMat_file, "TotalStatCovMatrix_AllVarNorm_rebinnedA",{{1, 228}});

  eft.drawCovMat(output_dir);

  return 0;
}
