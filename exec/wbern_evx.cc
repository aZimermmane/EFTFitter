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

//g++ `root-config --cflags --evelibs` -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I ./src -o ./exec/wbern_evx_total ./src/EFTFitter.cc ./exec/wbern_evx.cc

//g++ `root-config --cflags --evelibs` -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I ./src -o ./exec/wbern_evx_test ./src/EFTFitter.cc ./exec/wbern_evx.cc

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

bool independentSet(const std::vector<int> &v_index, const int &index) {
  // simple check if the set constitutes an independent measurement or not
  // eg set of b1k, b1r and b1n is not, since b1k^2 + b1r^2 + b1n^2 = 1
  // ie weed out anything that can be written as an identity
  static const std::vector<std::vector<int>> v_dependent = {
    {19},
    {20},
    {24},
    {25},
    {26},
    {27},
    {10, 11, 12, 34}, // cHel = cii / 3
    {11, 12, 15, 30, 34}, //ckk -> cPnr cnrP
    {11, 12, 16, 31, 34}, //ckk -> cMnr cnrM
    {11, 12, 30, 31, 34}, //ckk -> cnrP cnrM
    {12, 15, 21, 30, 34}, //ckk -> cPnr cnrP, crr -> chan
    {12, 15, 22 ,30, 34}, //ckk -> cPnr cnrP, crr -> ctra
    {12, 15, 23, 30, 34}, //ckk -> cPnr cnrP, crr -> csca
    {12, 16, 21, 31, 34}, //ckk -> cMnr cnrM, crr -> chan
    {12, 16, 22, 31, 34}, //ckk -> cMnr cnrM, crr -> ctra
    {12, 16, 23, 31, 34}, //ckk -> cMnr cnrM, crr -> csca
    {10, 13, 21, 28, 34}, // cnn -> cPrk crkP, crr -> chan
    {10, 13, 22, 28, 34}, // cnn -> cPrk crkP, crr -> ctra
    {10, 13, 23, 28, 34}, // cnn -> cPrk crkP, crr -> csca
    {10, 14, 21, 29, 34}, // cnn -> cMrk crkM, crr -> chan
    {10, 14, 22, 29, 34}, // cnn -> cMrk crkM, crr -> ctra
    {10, 14, 23, 29, 34}, // cnn -> cMrk crkM, crr -> csca

    // {11, 12, 15, 30, 34}, //ckk -> cPnr + cnrP
    // {11, 15, 21, 30, 34}, //ckk -> cPnr + cnrP , crr -> chan
    // {11, 15, 22, 30, 34}, //ckk = cPnr + cnrP , crr -> ctra
    // {11, 15, 23 , 30, 34}, //ckk -> cPnr + cnrP , crr -> csca
    // {10, 11, 14, 29, 34}, // cnn -> cPrk + crkP
    // {10, 11, 14, 29, 34}, // cnn -> cPrk + crkP
    // {10, 11, 14, 29, 34}, // cnn -> cPrk + crkP

    {11, 12, 19, 34}, // cHel = cii / 3
    {10, 12, 20, 34}, // cHel = cii / 3
    {12, 19, 20, 34}, // cHel = cii / 3


    {10, 11, 12, 21}, // cHan = cii / 3
    {11, 12, 19, 21}, // cHan = cii / 3
    {10, 12, 20, 21}, // cHan = cii / 3
    {19, 12, 20, 21}, // cHan = cii / 3

    {10, 11, 12, 22}, // cSca = cii / 3
    {11, 12, 19, 22}, // cSca = cii / 3
    {10, 12, 20, 22}, // cSca = cii / 3
    {19, 12, 20, 22}, // cSca = cii / 3

    {10, 11, 12, 23}, // cTra = cii / 3
    {11, 12, 19, 23}, // cTra = cii / 3
    {10, 12, 20, 23}, // cTra = cii / 3
    {19, 12, 20, 23}, // cTra = cii / 3

    {10, 11, 12, 24}, // ckjL = cii / 3
    {11, 12, 19, 24}, // ckjL = cii / 3
    {10, 12, 20, 24}, // ckjL = cii / 3
    {19, 12, 20, 24}, // ckjL = cii / 3

    {10, 11, 12, 25}, // cqrL = cii / 3
    {19, 12, 20, 25}, // cqrL = cii / 3
    {19, 12, 20, 25}, // cqrL = cii / 3
    {19, 12, 20, 25}, // cqrL = cii / 3

    {21, 22, 23, 34},  // CHan, cTra, csca, Chels alike
    {22, 23, 24, 34},  // ckjL,  cTra,  csca,  Chels alike
    {21, 23, 25, 34},  // CHan,  cqjL,  csca,  chel
    {23, 24, 25, 34},  // ckjL,  cqjL,  csca,  chel

    {10, 21, 34},  // ckk chan chel
    {19, 21, 34},  // ckq chan chel
    {10, 26, 34},  // ckk ckjL chel
    {19, 26, 23},  // ckq ckjL chel
    {11, 23, 34},  // crr csca chel
    {20, 23, 34},  // crj csca chel
    {11, 27, 34},  // crr crjL chel
    {20, 27, 34},  // crj crjL chel
    {12, 22, 34},  // cnn ctra chel

    {10, 22, 23}, // ckk csca ctra
    {15, 22, 23, 30}, //ckk -> cPnr cnrP, csca ctra
    {16, 22, 23, 31}, //ckk -> cMnr cnrM, csca ctra

    {11, 21, 22}, //crr chan ctra
    {17, 21, 22, 32}, // crr -> cPnk cnkP, chan ctra
    {18, 21, 22, 33}, // crr -> cMnk cnkM, chan ctra

    {12, 21, 23}, //cnn chan csca
    {13, 28, 21, 23}, //cnn->cPrk crkP chan csca
    {14, 29, 21, 23}, //cnn->cMrk crkM chan csca

    {11, 12, 21, 34},  //crr cnn chan chel
    {20, 12, 21, 34},  //crj cnn chan chel
    {11, 12, 26, 34},  //crr cnn ckjL chel
    {20, 12, 26, 34},  //crj cnn ckjL chel
    {10, 12, 23, 34},  // ckk cnn csca chel
    {10, 12, 23, 34},  // ckk cnn csca chel
    {19, 12, 23, 34},  // ckq cnn csca chel
    {10, 12, 27, 34},  // ckk cnn crjL chel
    {19, 12, 27, 34},  // ckq cnn crjL chel
    {10, 11, 22, 34},  // ckk crr ctra chel
    {10, 20, 22, 34},  // ckk crq ctra chel

    {10, 12, 21, 22, 34},  //combinations chan, ctra, cij giving chel
    {12, 19, 21, 22, 34},
    // {12, 15, 21, 22, 30, 34}, //combinations chan, ctra, cijPM giving chel
    // {12, 16, 21, 22, 31, 34},
    // {10, 13, 21, 22, 28, 34},
    // {13, 19, 21, 22, 28, 34},
    // {10, 14, 21, 22, 29, 34},
    // {14, 19, 21, 22, 29, 34},

    {10, 12, 26, 22, 34},  //combinations ckjL, ctra, cij giving chel
    {12, 19, 26, 22, 34},
    // {12, 15, 26, 22, 30, 34}, //combinations ckjL, ctra, cijPM giving chel
    // {12, 16, 26, 22, 31, 34},
    // {10, 13, 26, 22, 28, 34},
    // {13, 19, 26, 22, 28, 34},
    // {10, 14, 26, 22, 29, 34},
    // {14, 19, 26, 22, 29, 34},

    {10, 11, 21, 23, 34},                  //combinations chan, csca, cij giving chel
    {10, 20, 21, 23, 34},
    {11, 19, 21, 23, 34},
    {19, 20, 21, 23, 34},
    // {11, 15, 21, 23, 30, 34}, //combinations chan, csca, cijPM giving chel
    // {15, 20, 21, 23, 30, 34},
    // {11, 16, 21, 23, 31, 34},
    // {16, 20, 21, 23, 31, 34},
    // {10, 17, 21, 23, 32, 34},
    // {17, 19, 21, 23, 32, 34},
    // {10, 18, 21, 23, 33, 34},
    // {18, 19, 21, 23, 33, 34},

    {10, 11, 26, 23, 34},//combinations ckjl, csca, cij giving chel
    {10, 20, 26, 23, 34},
    {11, 19, 26, 23, 34},
    {19, 20, 26, 23, 34},
    // {11, 15, 26, 23, 30, 34},//combinations ckjl, csca, cijPM giving chel
    // {20, 15, 26, 23, 30, 34},
    // {11, 16, 26, 23, 31, 34},
    // {16, 10, 26, 23, 31, 34},
    // {10, 17,26, 23, 32, 34},
    // {17, 19, 26, 23, 32, 34},
    // {10, 18, 26, 23, 33, 34},
    // {18, 19, 26, 23, 33, 34},

    {10, 11, 21, 27, 34},   //combinations chan, crqL, cij giving chel
    {10, 20, 21, 27, 34},
    {11, 19, 21, 27, 34},
    {19, 20, 21, 27, 34},
    // {11, 15, 21, 27, 30, 34},//combinations chan, crqL, cijPM giving chel
    // {15, 20, 21, 27, 30, 34},
    // {11, 16, 21, 27, 31, 34},
    // {16, 10, 21, 27, 31, 34},
    // {10, 17, 21, 27, 32, 34},
    // {19, 17, 21, 27, 32, 34},
    // {10, 18, 21, 27, 33, 34},
    // {18, 19, 21, 27, 33, 34},

    {10, 11, 26, 27, 34},  //combinations crjL, crqL, cij giving chel
    {10, 20, 26, 27, 34},
    {11, 19, 26, 27, 34},
    {19, 20, 26, 27, 34},
    // {11, 15, 26, 27, 30, 34}, //combinations crjL, crqL, cijPM giving chel
    // {15, 10, 26, 27, 30, 34},
    // {11, 16, 26, 27, 31, 34},
    // {16, 20, 26, 27, 31, 34},
    // {10, 17, 26, 27, 32, 34},
    // {17, 19, 26, 27, 32, 34},
    // {10, 18, 26, 27, 33, 34},
    // {18, 19, 26, 27, 33, 34},

    {11, 12, 22, 23, 34}, //combinations ctra, csca, cij giving chel
    {12, 20, 22, 23, 34},
    // {12, 17, 22, 23, 32, 34}, //combinations ctra, csca, cijPM giving chel
    // {12, 18, 22, 23, 33, 34},
    // {11, 22, 23, 24, 28, 34},
    // {11, 22, 23, 25, 29, 34},
    // {11, 13, 22, 23, 28, 34},
    // {11, 14, 22, 23, 29, 24},
    // {20, 22, 23, 24, 28, 34},
    // {20, 22, 23, 25, 29, 34},
    // {13, 20, 22, 23, 28, 34},
    // {14, 20, 22, 23, 29, 24},

    {11, 12, 22, 27, 34}, //combinations ctra, crqL, cij giving chel
    {12, 20, 22, 27, 34},
    // {12, 17, 22, 27, 32, 34}, //combinations ctra, crqL, cijPM giving chel
    // {12, 18, 22, 27, 33, 34},
    // {11, 22, 24, 27, 28, 34},
    // {11, 22, 25, 27, 29, 34},
    // {11, 13, 22, 27, 28, 34},
    // {11, 14, 22, 27, 29, 24},
    // {20, 22, 24, 27, 28, 34},
    // {20, 22, 25, 27, 29, 34},
    // {13, 20, 22, 27, 28, 34},
    // {14, 20, 22, 27, 29, 24},

    {0, 2, 4}, // as above, krn
    {1, 3, 5},
    {4, 6, 8}, // njq
    {5, 7, 9},
    {0, 4, 8}, // knq
    {1, 5, 9},
    {2, 4, 6}, // rnj
    {3, 5, 7},
    {0, 6}, // le stars
    {1, 7},
    {2, 8},
    {3, 9},
    {10,19},
    {11,20},
    {14, 16, 18, 34}, // sHel^2 = cM_ij^2 ??
    {14, 16, 18, 21}, // sHel^2 = cM_ij^2 ??
    {14, 16, 18, 26}, // sHel^2 = cM_ij^2 ??
    {14, 16, 18, 22}, // sHel^2 = cM_ij^2 ??
    {14, 16, 18, 23}, // sHel^2 = cM_ij^2 ??
    {14, 16, 18, 27}, // sHel^2 = cM_ij^2 ??
    {14, 16, 18, 24}, // sHel^2 = cM_ij^2 ??
    {13, 15, 17, 34}, // maybe related ??
    {13, 15, 17, 21}, // maybe related ??
    {13, 15, 17, 26}, // maybe related ??
    {13, 15, 17, 22}, // maybe related ??
    {13, 15, 17, 23}, // maybe related ??
    {13, 15, 17, 27}, // maybe related ??
    {13, 15, 17, 24}, // maybe related ??
    {13, 14}, // for any ij pair, forbid both + and - being together
    {15, 16},
    {17, 18},
    {24, 25},
    {0, 1, 10}, // cii = b1i * b2i
    {0, 1, 19},
    {2, 3, 11},
    {2, 3, 20},
    {4, 5, 12},
    {6, 7, 10},
    {6, 7, 19},
    {8, 9, 11},
    {8, 9, 20},
    {0, 7, 10}, // as above up to a sign
    {0, 7, 19},
    {1, 6, 10},
    {1, 6, 19},
    {2, 9, 11},
    {2, 9, 20},
    {3, 8, 11},
    {3, 8, 20},
    {0, 1, 11, 12, 34}, // cHel = cii / 3 but with one bli pair
    {0, 1, 11, 12, 21},
    {0, 1, 11, 12, 22},
    {0, 1, 11, 12, 23},
    {0, 1, 11, 12, 27},
    {0, 1, 11, 12, 24},
    {2, 3, 10, 12, 34},
    {2, 3, 10, 12, 21},
    {2, 3, 10, 12, 22},
    {2, 3, 10, 12, 23},
    {2, 3, 10, 12, 24},
    {4, 5, 10, 11, 34},
    {4, 5, 10, 11, 21},
    {4, 5, 10, 11, 22},
    {4, 5, 10, 11, 23},
    {4, 5, 10, 11, 24},
    {6, 7, 11, 12, 34},
    {6, 7, 11, 12, 21},
    {6, 7, 11, 12, 22},
    {6, 7, 11, 12, 23},
    {6, 7, 11, 12, 24},
    {8, 9, 10, 12, 34},
    {8, 9, 10, 12, 21},
    {8, 9, 10, 12, 22},
    {8, 9, 10, 12, 23},
    {8, 9, 10, 12, 24},
    {0, 7, 11, 12, 34},
    {0, 7, 11, 12, 21},
    {0, 7, 11, 12, 22},
    {0, 7, 11, 12, 23},
    {0, 7, 11, 12, 24},
    {1, 6, 11, 12, 34},
    {1, 6, 11, 12, 21},
    {1, 6, 11, 12, 22},
    {1, 6, 11, 12, 23},
    {1, 6, 11, 12, 24},
    {2, 9, 11, 12, 34},
    {2, 9, 11, 12, 21},
    {2, 9, 11, 12, 22},
    {2, 9, 11, 12, 23},
    {2, 9, 11, 12, 24},
    {3, 8, 11, 12, 34},
    {3, 8, 11, 12, 21},
    {3, 8, 11, 12, 22},
    {3, 8, 11, 12, 23},
    {3, 8, 11, 12, 24},
    {10, 11, 13}, // cii + cjj related to cij +- cji? cva doesnt like it...
    {11, 12, 15},
    {10, 12, 17},
    {10, 11, 14},
    {11, 12, 16},
    {10, 12, 18},
    {0, 1, 11, 13}, // as above but with one bli pair
    {2, 3, 10, 13},
    {6, 7, 11, 13},
    {8, 9, 10, 13},
    {0, 7, 11, 13},
    {1, 6, 11, 13},
    {2, 9, 10, 13},
    {3, 8, 10, 13},
    {2, 3, 12, 15},
    {4, 5, 11, 15},
    {8, 9, 12, 15},
    {2, 9, 12, 15},
    {3, 8, 12, 15},
    {0, 1, 12, 17},
    {4, 5, 10, 17},
    {6, 7, 12, 17},
    {0, 7, 12, 17},
    {1, 6, 12, 17},
    {0, 1, 11, 14},
    {2, 3, 10, 14},
    {6, 7, 11, 14},
    {8, 9, 10, 14},
    {0, 7, 11, 14},
    {1, 6, 11, 14},
    {2, 9, 10, 14},
    {3, 8, 10, 14},
    {2, 3, 12, 16},
    {4, 5, 11, 16},
    {8, 9, 12, 16},
    {2, 9, 12, 16},
    {3, 8, 12, 16},
    {0, 1, 12, 18},
    {4, 5, 10, 18},
    {6, 7, 12, 18},
    {0, 7, 12, 18},
    {1, 6, 12, 18},
    {12, 13, 28}, //c_rkP
    {4, 5, 13, 28}, // as above but with one bli pair
    {0, 3, 17, 28},
    {1, 2, 17, 28},
    {12, 14, 29}, //crkM
    {4, 5, 13, 29}, // as above but with one bli pair
    {0, 3, 17, 29},
    {1, 2, 17, 29},

    {10, 15, 30}, //c_nrP
    {15, 19, 30},
    {0, 1, 15, 30}, // as above but with one bli pair
    {6, 7, 15, 30},
    {0, 7, 15, 30},
    {1, 6, 15, 30},
    // {2, 3, 4, 5, 10, 30}, // as above but with two/threre bli pair
    // {2, 3, 4, 5, 19, 30},
    // {0, 1, 2, 3, 4, 5, 30},
    // {0, 7, 2, 3, 4, 5, 30},
    // {1, 6, 2, 3, 4, 5, 30},
    // {6, 7, 2, 3, 4, 5, 30},
    //
    // {2, 4, 5, 9, 10, 30},
    // {2, 4, 5, 9, 19, 30},
    // {0, 1, 2, 4, 5, 9, 30},
    // {0, 7, 2, 4, 5, 9, 30},
    // {1, 6, 2, 4, 5, 9, 30},
    // {6, 7, 2, 4, 5, 9, 30},
    //
    // {3, 4, 5, 8, 10, 30},
    // {3, 4, 5, 8, 19, 30},
    // {0, 1, 3, 4, 5, 8, 30},
    // {0, 7, 3, 4, 5, 8, 30},
    // {1, 6, 3, 4, 5, 8, 30},
    // {6, 7, 3, 4, 5, 8, 30},
    //
    // {2, 5, 6, 9, 10, 30},
    // {2, 5, 6, 9, 19, 30},
    // {0, 1, 2, 5, 6, 9, 30},
    // {0, 7, 2, 5, 6, 9, 30},
    // {1, 6, 2, 5, 6, 9, 30},
    // {6, 7, 2, 5, 6, 9, 30},

    {10, 16, 31}, //c_nrM
    {16, 19, 31},
    {0, 1, 16, 31}, // as above but with one bli pair
    {6, 7, 16, 31},
    {0, 7, 16, 31},
    {1, 6, 16, 31},
    // {2, 3, 4, 5, 10, 31}, // as above but with two/threre bli pair
    // {2, 3, 4, 5, 19, 31},
    // {0, 1, 2, 3, 4, 5, 31},
    // {0, 7, 2, 3, 4, 5, 31},
    // {1, 6, 2, 3, 4, 5, 31},
    // {6, 7, 2, 3, 4, 5, 31},

    // {2, 4, 5, 9, 10, 31},
    // {2, 4, 5, 9, 19, 31},
    // {0, 1, 2, 4, 5, 9, 31},
    // {0, 7, 2, 4, 5, 9, 31},
    // {1, 6, 2, 4, 5, 9, 31},
    // {6, 7, 2, 4, 5, 9, 31},

    // {3, 4, 5, 8, 10, 31},
    // {3, 4, 5, 8, 19, 31},
    // {0, 1, 3, 4, 5, 8, 31},
    // {0, 7, 3, 4, 5, 8, 31},
    // {1, 6, 3, 4, 5, 8, 31},
    // {6, 7, 3, 4, 5, 8, 31},
    //
    // {2, 5, 6, 9, 10, 31},
    // {2, 5, 6, 9, 19, 31},
    // {0, 1, 2, 5, 6, 9, 31},
    // {0, 7, 2, 5, 6, 9, 31},
    // {1, 6, 2, 5, 6, 9, 31},
    // {6, 7, 2, 5, 6, 9, 31},
    {11, 17, 32}, //c_nkP
    {2, 3, 17, 32}, // as above but with one bli pair
    {8, 9, 17, 32},
    {2, 9, 17, 32},
    {3, 8, 17, 32},
    {0, 5, 11, 32},
    {1, 4, 11, 32},
    {4, 7, 11, 32},
    {5, 6, 11, 32},
    {11, 18, 33}, //c_nkM
    {2, 3, 17, 33}, // as above but with one bli pair
    {8, 9, 17, 33},
    {2, 9, 17, 33},
    {3, 8, 17, 33},
    {0, 5, 11, 33},
    {1, 4, 11, 33},
    {4, 7, 11, 33},
    {5, 6, 11, 33},
  };

  std::vector<int> v_tmp(v_index);
  v_tmp.push_back(index);
  std::sort(std::begin(v_tmp), std::end(v_tmp));

  // hack to disable the check
  //return true;

  for (auto &v_dep : v_dependent) {
    if (std::includes(std::begin(v_tmp), std::end(v_tmp), std::begin(v_dep), std::end(v_dep)))
      return false;
  }

  return true;
}

int main(int argc, char** argv) {
  if (argc < 2)
    return 0;

  using EFT = EFTFitter;

  // common flags for an evolution test
  // fit with data or SM, check evolution of all or linear part, 2 sigma or 1 sigma interval width as figure of merit
  const bool useData = false, useAll = false, useSig2 = false;

  // construction done with the key we want to treat as data, lambda, fit mode and stat mode (optionally sum of shape templates)
  //EFT eft("data", 1., EFT::Fit::shape, EFT::Stat::xsec); // Modified by andre 01.02.23
  EFT eft("data", 1., EFT::Fit::hybrid, EFT::Stat::xsec);
  //const std::string inDir = "/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/wbern_0519/root/", opName = argv[1];
  const std::string inDir = "/nfs/dust/cms/user/zimermma/EFTFitter/wbern_0314/root/", opName = argv[1];
  //const std::string inDir = "../histograms/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/", opName = argv[1];
  const std::string covMatrix ="Systematics_AllVars_1D_228x228_1000PE/TotalStatCovMatrix_AllVarNorm_rebinnedA";// Modified by andre 01.02.23
  //const std::string covMatrix ="TotalStatCovMatrix_AllVarNorm_rebinnedA";
  //const std::string covMatrix ="TotalSystCovMatrix_AllVarNorm_rebinnedA";
  //const std::string covMatrix ="TotalSystStatCovMatrix_AllVarNorm_rebinnedA";
  //const std::string covMatrix = "TotalStat_shape_a_drop_bin_1";
  const std::string outDir = "/nfs/dust/cms/user/zimermma/EFTFitter/fit_output/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV/" + covMatrix + "/" + opName + "/";

  gSystem->Exec(("mkdir -p "+ outDir).c_str() );

  // auto mkdir = std::experimental::filesystem::create_directories(outDir);
  // if (mkdir) {
  //       std::cout << "created output directory: " << std::endl;
  //       std::system("tree x_tmp");
  //   } else {
  //       std::cout << "create_directories() failed" << std::endl;
  //   }

  // goes forever due to 2D degeneracy: {cvv, c1} {cva, c3}, {cav, c123}
  const std::map<std::string, double> m_range = {{"ctG",0.2},{"ut", 0.04}, {"cvv", 0.06}, {"c1", 0.5},
                                                 {"dt", 0.04}, {"cmm", 0.08},
                                                 {"cva", 0.5}, {"c3", 0.5},
                                                 {"cav", 0.1}, {"c123", 0.5},
                                                 {"cmp", 0.05}};

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
                                           "cPrj","cM_rj",
                                           "c_kjL", "c_rqL",
                                           "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM",
                                           "cHel"};



  //const std::vector<std::string> v_hStr = {"cP_rk"};
                                           //need to add an input file that contains b1j
  const std::vector<int> v_all_var = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
  //const std::vector<int> v_all_var = {13};

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
  auto evoFile = std::make_unique<TFile>((outDir + opName + "_evolution.root").c_str(), "recreate");
  std::vector<int> v_fit_var;

  //for (int iIt = 0; iIt < v_hStr.size(); ++iIt) {
  for (int iIt = 0; iIt < 5; ++iIt) {
  //for (int iIt = 0; iIt < 1; ++iIt) {
    std::cout << "EFTFitter: iteration " << iIt << " starting..." << std::endl;
    std::vector<std::array<std::unique_ptr<TGraphAsymmErrors>, 2>> v_current_result;


    for (int iFit = 0; iFit < v_hStr.size(); ++iFit) {
      const int iVar = v_all_var.at(iFit);

      if ( std::count(std::begin(v_fit_var), std::end(v_fit_var), iVar) or !independentSet(v_fit_var, iVar)) {
        std::cout << "Entering the IF " << std::endl;
        //std::cout << "  std::count(std::begin(v_fit_var), std::end(v_fit_var): " << std::count(std::begin(v_fit_var), std::end(v_fit_var) << std::endl;
        //std::cout << "  std::end(v_fit_var)" << std::end(v_fit_var) << std::endl;
        std::cout << "  iVar" << iVar << std::endl;
        v_current_result.emplace_back(fitResult("", ""));
        continue;
      }

      std::vector<int> v_iEB(v_fit_var);
      v_iEB.push_back( iVar );
      std::sort(std::begin(v_iEB), std::end(v_iEB));

      std::cout << "EFTFitter: starting fit on variable " << v_hStr.at(iFit) << " in iteration " << iIt << std::endl;

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
      auto f_emulate = [&v_iEB, &rate_zero, &binToIgnore, &v_spinvar] (const auto &v_binC) {
        // convert the coeffs into shapes using the known parametrization
        std::vector<std::array<double, 2>> v_tmpC = {rate_zero, rate_zero};
        TF1 *f_ref = v_spinvar.at(0).get();

        for (int iC = 2; iC < v_binC.size(); ++iC) {
          if (!std::count(std::begin(v_iEB), std::end(v_iEB), iC - 2)) continue;

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
      for (auto &var : v_iEB)
      {
        covMat_binRange.push_back({(2+(var) * 6), ((var + 1) * 6)}); //modified by Andre 01.02. (always ignore the first bin of each dib)
        //covMat_binRange.push_back({((var + 1) * 6) - 5, ((var + 1) * 6 -1)}); //modified by Andre 01.02. (always ignore the last bin of each dib)

	std::cout << "var: " << var << std::endl;
        //covMat_binRange.push_back({((var + 1) * 5) - 4, ((var + 1) * 5)});
      };
      //eft.readCovMatRoot("finalcov", "/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/covariance_matrix/covmat_190114.root", covMatrix, covMat_binRange);
      //eft.readCovMatRoot("finalcov", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/covmat_190114.root", "TotalStat_shape_a", covMat_binRange);
      //eft.readCovMatRoot("totalStat", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_132x132.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
      eft.readCovMatRoot("finalcov", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_228x228_1000PE.root", "TotalStatCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
    //  eft.readCovMatRoot("totalSyst", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/Systematics_AllVars_1D_228x228_1000PE.root", "TotalSystCovMatrix_AllVarNorm_rebinnedA", covMat_binRange);
      //eft.readCovMatRoot("finalcov", "/nfs/dust/cms/user/zimermma/EFTFitter/inputs/covariance_matrix/covMatSyst_purdue_full2016.root", covMatrix, covMat_binRange);
      //eft.makeFinalCovMat({"totalStat","totalSyst"});

      eft.drawCovMat(outDir,{},true);
      gSystem->Exec( ("mv " + outDir + "cov_finalcov.txt " + outDir + "covMat_" + v_hStr.at(iFit) + "_iter_" + toStr(iIt) + ".txt").c_str() );
      gSystem->Exec( ("mv " + outDir + "cov_finalcov.pdf " + outDir + "covMat_" + v_hStr.at(iFit) + "_iter_" + toStr(iIt) + ".pdf").c_str() );
      gSystem->Exec( ("rm " + outDir + "cov_all.root").c_str() );

      std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
      vt_keySampleLegend.push_back({opName + "_-1", EFT::Sample::all, opName + " -1"});
      vt_keySampleLegend.push_back({opName + "_1", EFT::Sample::all, opName + " 1"});

      vt_keySampleLegend.push_back({opName + "_0", EFT::Sample::all, "SM"});
      vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});

      eft.drawHistogram(vt_keySampleLegend,
                        outDir + opName + "_var_" + v_hStr.at(iFit) + "_iter_" + toStr(iIt) + "_shape",
                        "Fraction", "Index", -0.4999, 0.4999, 0.0001, 1.9999,
                        false, "", false, "none");

      eft.listKeyToFit({ {opName, v_opPoint}});
      eft.computeFitChi2(v_sample);

      eft.draw1DChi2({ {opName, {opName, {/* op range in min, max */}, {0., 9.999}, {-opRange + eps, opRange - eps} }} }, outDir, v_sample);

      eft.clearContent();

      v_current_result.emplace_back(fitResult(opName, outDir + opName + "_dChi2.root", iIt, iFit, useAll));
      if (v_current_result.back().at(0) == nullptr)
        std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " in iteration " << iIt << " doesn't result in a constraint!" << std::endl;
      else
        std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " in iteration " << iIt << " saved for evolution scan..." << std::endl;

      gSystem->Exec( ("rm " + outDir + opName + "_dChi2.root").c_str() );

    }

    // first find the best variable in a given iter by minimum sig2 or sig1 width
    int iMinF = -999, iMinV = -999;
    double width = 9999.;
    for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
      if (v_current_result.at(iRes).at(0) == nullptr) continue;
      auto &graph = (useSig2) ? v_current_result.at(iRes).at(1) : v_current_result.at(iRes).at(0);

      if (graph->GetErrorXlow(0) + graph->GetErrorXhigh(0) < width) {
        width = graph->GetErrorXlow(0) + graph->GetErrorXhigh(0);
        iMinF = iRes;
        iMinV = v_all_var.at(iRes);
      }
    }

    if (iMinF == -999) {
      std::cout << "EFTFitter: iteration " << iIt << " aborted; no further variables can make up an independent set." << std::endl;
      break;
    }

    v_fit_var.push_back(iMinV);
    std::cout << "EFTFitter: variable " << v_hStr.at(iMinF) << " with width " << width << " is the best variable in iteration " << iIt << std::endl;

    evoFile->cd();
    for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
      if (v_current_result.at(iRes).at(0) == nullptr) continue;

      if (iRes == iMinF) {
        v_current_result.at(iRes).at(0)->SetName( ("iter_" + toStr(iIt) + "_best_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma1").c_str() );
        v_current_result.at(iRes).at(0)->Write( ("iter_" + toStr(iIt) + "_best_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma1").c_str() );

        v_current_result.at(iRes).at(1)->SetName( ("iter_" + toStr(iIt) + "_best_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma2").c_str() );
        v_current_result.at(iRes).at(1)->Write( ("iter_" + toStr(iIt) + "_best_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma2").c_str() );
      }
      else {
        v_current_result.at(iRes).at(0)->SetName( ("iter_" + toStr(iIt) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma1").c_str() );
        v_current_result.at(iRes).at(0)->Write( ("iter_" + toStr(iIt) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma1").c_str() );

        v_current_result.at(iRes).at(1)->SetName( ("iter_" + toStr(iIt) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma2").c_str() );
        v_current_result.at(iRes).at(1)->Write( ("iter_" + toStr(iIt) + "_fit_" + sample + "_" + toStr(v_all_var.at(iRes)) + "_sigma2").c_str() );
      }
    }

    std::cout << "EFTFitter: iteration " << iIt << " completed!" << std::endl;
  }

  return 0;
}
