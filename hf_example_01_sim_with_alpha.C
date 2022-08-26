#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooRealProxy.h>
#include <RooAbsPdf.h>
#include <TROOT.h>
#include <TSystem.h>
#include <RooMinimizer.h>
#include <RooFitResult.h>
#include <RooListProxy.h>
#include <RooStats/HistFactory/Measurement.h>
#include <RooStats/HistFactory/MakeModelAndMeasurementsFast.h>
#include <RooStats/ModelConfig.h>

#include <TMath.h>

#ifdef BENCH
#include "benchmark/benchmark.h"
#endif

#include <clad/Differentiator/Differentiator.h>

#include "include/RooSimClasses.h"
#include "include/derivative.h"
#include <string>
#include <cmath>

std::unique_ptr<RooWorkspace> makeHistFactoryWorkspace(int nChannels)
{
   using namespace RooStats;
   using namespace HistFactory;

   // Create the measurement
   Measurement meas("meas", "meas");

   meas.SetOutputFilePrefix("jonas");
   meas.SetPOI("SigXsecOverSM");
   meas.AddConstantParam("alpha_chan0_syst1");
   meas.AddConstantParam("alpha_chan1_syst1");
   meas.AddConstantParam("Lumi");

   meas.SetLumi(1.0);
   meas.SetLumiRelErr(0.10);
   meas.SetExportOnly(false);
   meas.SetBinHigh(2);

   for (std::size_t iChannel = 0; iChannel < nChannels; ++iChannel) {
      // Create a channel
      std::string chanStr = std::to_string(iChannel);
      Channel chan(("channel" + chanStr).c_str());

      chan.SetStatErrorConfig(0.05, "Poisson");

      // set data
      auto dataHist = new TH1F{"data_hist", "data_hist", 2, 1.0, 2.0};
      dataHist->SetBinContent(1, 122);
      dataHist->SetBinContent(2, 112);
      HistFactory::Data data;
      data.SetHisto(dataHist);
      chan.SetData(data);

      // Now, create some samples

      // Create the signal sample
      Sample signal("signal");
      auto sigHist = new TH1F{"sig_hist", "sig_hist", 2, 1.0, 2.0};
      sigHist->SetBinContent(1, 20);
      sigHist->SetBinContent(2, 10);
      signal.SetHisto(sigHist);
      signal.AddOverallSys(("chan" + chanStr + "_syst1").c_str(), 0.95, 1.05);
      signal.AddNormFactor("SigXsecOverSM", 1, 0, 3);
      chan.AddSample(signal);

      // Background 1
      Sample background1("background1");
      auto bkg1Hist = new TH1F{"bkg1_hist", "bkg1_hist", 2, 1.0, 2.0};
      bkg1Hist->SetBinContent(1, 100);
      bkg1Hist->SetBinContent(2, 0);
      background1.SetHisto(bkg1Hist);
      auto bkg1UncertHist = new TH1F{"background1_statUncert", "background1_statUncert", 2, 1.0, 2.0};
      bkg1UncertHist->SetBinContent(1, 0.05);
      bkg1UncertHist->SetBinContent(2, 0.05);
      background1.GetStatError().Activate();
      background1.GetStatError().SetUseHisto();
      background1.GetStatError().SetErrorHist(bkg1UncertHist);
      background1.AddOverallSys(("chan" + chanStr + "_syst2").c_str(), 0.95, 1.05);
      chan.AddSample(background1);

      // Background 2
      Sample background2("background2");
      auto bkg2Hist = new TH1F{"bkg2_hist", "bkg2_hist", 2, 1.0, 2.0};
      bkg2Hist->SetBinContent(1, 0);
      bkg2Hist->SetBinContent(2, 100);
      background2.SetHisto(bkg2Hist);
      background2.ActivateStatError();
      background2.AddOverallSys(("chan" + chanStr + "_syst3").c_str(), 0.95, 1.05);
      chan.AddSample(background2);

      meas.AddChannel(chan);
   }

   return std::unique_ptr<RooWorkspace>{MakeModelAndMeasurementFast(meas)};
}

std::string generateNLLCode()
{

   contextManager ctx;

   ExRooRealVar X(ctx, "x", "{1.25, 1.75}", 2);
   ExRooConst B2Eps(ctx, 1);
   ExRooConst B1Eps(ctx, 1);
   ExRooRealVar NomLumi(ctx, "nominalLumi", "1");
   ExRooConst RandConst(ctx, 0.1);
   // ExRooRealVar AlphaSys(ctx, "alphaSys", "0");
   ExRooConst RandConst2(ctx, 1);
   // ExRooRealVar NomAlphaSys(ctx, "nomAlphaSys", "0");
   ExRooRealVar GammaB1(ctx, "gammaB1"); // Input
   ExRooRealVar NomGammaB1(ctx, "nomGammaB1", "400");
   ExRooRealVar GammaB2(ctx, "gammaB2");         // Input
   ExRooRealVar SigXOverM(ctx, "SigXsecOverSM"); // Input
   ExRooRealVar Lumi(ctx, "lumi");               // Input
   ExRooRealVar NomGammaB2(ctx, "nomGammaB2", "100");
   ExRooRealVar ChannelWeight(ctx, "cw", "1");
   ExRooRealVar ChannelWeight2(ctx, "cw2", "1");

   // Hist func part
   ExRooHistFunc Sig(ctx, &X, "sig", "{20, 10}", "{1, 1.5, 2}"), Bgk1(ctx, &X, "bgk1", "{100, 0}", "{1, 1.5, 2}"),
      Bgk2(ctx, &X, "bgk2", "{0, 100}", "{1, 1.5, 2}");
   ExRooParamHistFunc ParamHist(ctx, &X, {&GammaB1, &GammaB2});
   ExRooProduct constr1(ctx, {/* &AlphaSys, */ &SigXOverM, &Lumi});
   ExRooProduct constr2(ctx, {&Lumi, &B1Eps});
   ExRooProduct constr3(ctx, {&Lumi, &B2Eps});
   ExRooProduct Bgk1Shape(ctx, {&Bgk1, &ParamHist}), Bgk2Shape(ctx, {&Bgk2, &ParamHist});
   ExRooRealSum SumHF(ctx, "mu", {&Sig, &Bgk1Shape, &Bgk2Shape}, {&constr1, &constr2, &constr3});
   ExRooProduct ProdHF(ctx, {&SumHF});

   // NLL node
   ExRooNll2 Nll(ctx, &X, "2", &ProdHF, {122, 112});

   // Poisson B0 constraint
   ExRooProduct PoissProd1(ctx, {&NomGammaB1, &GammaB1});
   ExRooPoisson Poisson1(ctx, &NomGammaB1, &PoissProd1);
   // Poisson B1 constraint
   ExRooProduct PoissProd2(ctx, {&NomGammaB2, &GammaB2});
   ExRooPoisson Poisson2(ctx, &NomGammaB2, &PoissProd2);
   // Gaussian alpha constraint
   // ExRooGaussian GaussAlpha(ctx, &AlphaSys, &NomAlphaSys, &RandConst2);
   // Gaussian Lumi constraint
   ExRooGaussian GaussLumi(ctx, &Lumi, &NomLumi, &RandConst);

   // Constraint sum
   ExRooConstraintSum ConstrSum(ctx, {&Poisson1, &Poisson2, /* &GaussAlpha, */
                                      &GaussLumi});

   // Final root model node
   ExRooAddition Root(ctx, {&ConstrSum, &Nll});

   std::string code = Root.getCode();
   return "double nll(" + ctx.getParamList() + ") { \n" + code + " return " + Root.getResult() + ";\n}\n";
}

// Generated Code
// double nll(double lumi, double SigXsecOverSM, double gammaB1, double gammaB2)
// {
//    double nomGammaB1 = 400;
//    double nomGammaB2 = 100;
//    double nominalLumi = 1;
//    double constraint[3]{ExRooPoisson::poisson(nomGammaB1, (nomGammaB1 * gammaB1)),
//                         ExRooPoisson::poisson(nomGammaB2, (nomGammaB2 * gammaB2)),
//                         ExRooGaussian::gauss(lumi, nominalLumi, 0.100000)};
//    double cnstSum = 0;
//    double x[2]{1.25, 1.75};
//    double sig[2]{20, 10};
//    double binBoundaries1[3]{1, 1.5, 2};
//    double bgk1[2]{100, 0};
//    double binBoundaries2[3]{1, 1.5, 2};
//    double histVals[2]{gammaB1, gammaB2};
//    double bgk2[2]{0, 100};
//    double binBoundaries3[3]{1, 1.5, 2};
//    double weights[2]{122.000000, 112.000000};
//    for (int i = 0; i < 3; i++) {
//       cnstSum -= std::log(constraint[i]);
//    }
//    unsigned int b1, b2, b3;
//    double mu = 0;
//    double nllSum = 0;
//    double temp;
//    for (int iB = 0; iB < 2; iB++) {
//      b1 = ExRooHistFunc::getBin(binBoundaries1, x[iB]);
//      b2 = ExRooHistFunc::getBin(binBoundaries2, x[iB]);
//      b3 = ExRooHistFunc::getBin(binBoundaries3, x[iB]);
//      mu = 0;
//      mu += sig[b1] * (SigXsecOverSM * lumi);
//      mu += (bgk1[b2] * histVals[iB]) * (lumi * 1.000000);
//      mu += (bgk2[b3] * histVals[iB]) * (lumi * 1.000000);
//      temp = std::log((mu));
//      nllSum -= -(mu) + weights[iB] * temp;
//    }
//    return cnstSum + nllSum;
// }

double nll(double gammaB1, double gammaB2, double gamma1_B1, double gamma1_B2, double alpha_chan0_syst2,
           double alpha_chan0_syst3, double alpha_chan1_syst2, double alpha_chan1_syst3, double SigXsecOverSM)
{
   constexpr double lumi = 1.0;

   double nomGammaB1 = 400;
   double nomGammaB2 = 100;
   double nomGamma1_B1 = 400;
   double nomGamma1_B2 = 100;
   double nominalLumi = 1;
   double constraint[8]{
      ExRooPoisson::poisson(nomGammaB1, (nomGammaB1 * gammaB1)),
      ExRooPoisson::poisson(nomGammaB2, (nomGammaB2 * gammaB2)),
      ExRooPoisson::poisson(nomGamma1_B1, (nomGamma1_B1 * gamma1_B1)),
      ExRooPoisson::poisson(nomGamma1_B2, (nomGamma1_B2 * gamma1_B2)),
      // ExRooGaussian::gauss(lumi, nominalLumi, 0.100000),
      ExRooGaussian::gauss(alpha_chan0_syst2, 0.0, 1.000000),
      ExRooGaussian::gauss(alpha_chan0_syst3, 0.0, 1.000000),
      ExRooGaussian::gauss(alpha_chan1_syst2, 0.0, 1.000000),
      ExRooGaussian::gauss(alpha_chan1_syst3, 0.0, 1.000000),
   };
   double cnstSum = 0;
   double x[2]{1.25, 1.75};
   double sig0[2]{20, 10};
   double binBoundaries1[3]{1, 1.5, 2};
   double bgk0_1[2]{100, 0};
   double binBoundaries2[3]{1, 1.5, 2};
   double histVals0[2]{gammaB1, gammaB2};
   double bgk0_2[2]{0, 100};
   double binBoundaries3[3]{1, 1.5, 2};
   double nllSum0_weights[2]{122.000000, 112.000000};
   double sig1[2]{20, 10};
   double binBoundaries4[3]{1, 1.5, 2};
   double bgk1_1[2]{100, 0};
   double binBoundaries5[3]{1, 1.5, 2};
   double histVals1[2]{gamma1_B1, gamma1_B2};
   double bgk1_2[2]{0, 100};
   double binBoundaries6[3]{1, 1.5, 2};
   double nllSum1_weights[2]{122.000000, 112.000000};
   for (int i = 0; i < 8; i++) {
      cnstSum -= std::log(constraint[i]);
   }
   unsigned int b1, b2, b3, b4, b5, b6;
   double nllSum0 = 0;
   double mu0 = 0;
   double temp;
   for (int iB = 0; iB < 2; iB++) {
      mu0 = 0;
      b1 = ExRooHistFunc::getBin(binBoundaries1, x[iB]);
      b2 = ExRooHistFunc::getBin(binBoundaries2, x[iB]);
      b3 = ExRooHistFunc::getBin(binBoundaries3, x[iB]);
      mu0 += sig0[b1] * (SigXsecOverSM * lumi);
      mu0 += (bgk0_1[b2] * histVals0[iB]) * (lumi * 1.000000) * (1.0 + 0.05 * alpha_chan0_syst2);
      mu0 += (bgk0_2[b3] * histVals0[iB]) * (lumi * 1.000000) * (1.0 + 0.05 * alpha_chan0_syst3);
      temp = std::log((mu0));
      nllSum0 -= -(mu0) + nllSum0_weights[iB] * temp;
   }
   double nllSum1 = 0;
   double mu1 = 0;
   for (int iB = 0; iB < 2; iB++) {
      mu1 = 0;
      b4 = ExRooHistFunc::getBin(binBoundaries4, x[iB]);
      b5 = ExRooHistFunc::getBin(binBoundaries5, x[iB]);
      b6 = ExRooHistFunc::getBin(binBoundaries6, x[iB]);
      mu1 += sig1[b4] * (SigXsecOverSM * lumi);
      mu1 += (bgk1_1[b5] * histVals1[iB]) * (lumi * 1.000000) * (1.0 + 0.05 * alpha_chan1_syst2);
      mu1 += (bgk1_2[b6] * histVals1[iB]) * (lumi * 1.000000) * (1.0 + 0.05 * alpha_chan1_syst3);
      temp = std::log((mu1));
      nllSum1 -= -(mu1) + nllSum1_weights[iB] * temp;
   }
   return cnstSum + nllSum0 + nllSum1;
}

template <typename Func = void, typename Grad = void>
class RooGradFuncWrapper final : public RooAbsReal {
public:
   RooGradFuncWrapper(Func function, Grad gradient, RooWorkspace const &ws, std::vector<std::string> const &paramNames)
      : RooAbsReal{"RooGradFuncWrapper", "RooGradFuncWrapper"}, _paramProxies{"!params", "List of parameters", this},
        _func(function), _grad(gradient)
   {
      for (auto const &paramName : paramNames) {
         auto *var = ws.var(paramName.c_str());
         _paramProxies.add(*var);
         _params.emplace_back(var);
      }
   }
   RooGradFuncWrapper(const RooGradFuncWrapper &other, const char *name = nullptr)
      : RooAbsReal(other, name), _paramProxies("!params", this, other._paramProxies), _func(other._func),
        _grad(other._grad)
   {
   }

   TObject *clone(const char *newname) const override { return new RooGradFuncWrapper(*this, newname); }

   double defaultErrorLevel() const override { return 0.5; }

protected:
   double evaluate() const override
   {
      return _func(_params[0]->getValV(), _params[1]->getValV(), _params[2]->getValV(), _params[3]->getValV(),
                   _params[4]->getValV(), _params[5]->getValV(), _params[6]->getValV(), _params[7]->getValV(),
                   _params[8]->getValV());
   }

   void evaluateGradient(double *out) const override
   {
      // out[0] = 0;
      // out[1] = 0;
      // out[2] = 0;
      // out[3] = 0;
      // out[4] = 0;
      // out[5] = 0;
      // out[6] = 0;
      // out[7] = 0;
      // out[8] = 0;
      //_grad(_params[0]->getValV(), _params[1]->getValV(), _params[2]->getValV(), _params[3]->getValV(),
      //_params[4]->getValV(), _params[5]->getValV(),
      //_params[6]->getValV(), _params[7]->getValV(),
      //_params[8]->getValV(),
      //&out[0], &out[1], &out[2], &out[3], &out[4], &out[5],
      //&out[6], &out[7], &out[8],
      //);
   }

private:
   RooListProxy _paramProxies;
   std::vector<RooRealVar *> _params;
   Func _func;
   Grad _grad;
};

#ifdef BENCH
static void hf_example_01_sim(benchmark::State &state)
#else
int main()
#endif
{

   using namespace RooStats;
   using namespace RooFit;

   gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
   auto &msg = RooMsgService::instance();
   msg.setGlobalKillBelow(RooFit::WARNING);

   gInterpreter->Declare("#pragma cling optimize(2)");

   //   std::string func = generateNLLCode();
   //   std::cout << func.c_str();

   // clad::gradient(nll);

   // // FIXME :)
   // gInterpreter->ProcessLine(
   //     "#include "
   //     "\"/home/grimmyshini/ROOT/examples/roofit-clad-work/RooSimClasses.h\"");
   // gInterpreter->Declare(func.c_str());

   // gInterpreter->ProcessLine(
   //     "#include \"/home/grimmyshini/cern/src/tools/clad/include/clad/Differentiator/Differentiator.h\"");
   //   // calculate gradient
   //   gInterpreter->ProcessLine("clad::gradient(nll).dump();");

   //   // get correct gradient overload.
   //   gInterpreter->ProcessLine(
   //       "void (*ptr)(double, double, clad::array_ref<double>) = nll_grad;");
   //   // get the grad function pointer.
   //   auto gradObj =
   //       (void (*)(double, double,
   //                 clad::array_ref<double>))gInterpreter->ProcessLine("ptr");
   //   auto funcObj = (double (*)(doubl
   //  e,
   //   double))gInterpreter->ProcessLine("nll");

   std::unique_ptr<RooWorkspace> w = makeHistFactoryWorkspace(2);

   auto *mc = static_cast<ModelConfig *>(w->obj("ModelConfig"));

   RooArgSet params;
   if (mc->GetParametersOfInterest())
      params.add(*mc->GetParametersOfInterest());
   if (mc->GetNuisanceParameters())
      params.add(*mc->GetNuisanceParameters());

   // params.Print("t");

   std::unique_ptr<RooArgSet> originalParams{params.snapshot()};

   auto resetParameters = [&]() {
      params.assignFast(*originalParams);
      for (auto *param : static_range_cast<RooRealVar *>(params)) {
         param->setError(0.0);
      }
   };

   auto *pdf = w->pdf("simPdf");

   std::vector<std::string> paramNames{
      "gamma_stat_channel0_bin_0", "gamma_stat_channel0_bin_1", "gamma_stat_channel1_bin_0",
      "gamma_stat_channel1_bin_1", "alpha_chan0_syst2",         "alpha_chan0_syst3",
      "alpha_chan1_syst2",         "alpha_chan1_syst3",         "SigXsecOverSM",
   };

   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

   // Create the likelihood objects with RooFit and our generated code
   std::unique_ptr<RooAbsReal> nllRooFit{pdf->createNLL(*w->data("obsData"), Constrain(*mc->GetNuisanceParameters()),
                                                        GlobalObservables(*mc->GetGlobalObservables()))};
   std::unique_ptr<RooAbsReal> nllCodeGen{new RooGradFuncWrapper(nll, nll_grad, *w, paramNames)};

   // Create the minimizers for each configuration
   auto setupMinimizer = [](RooMinimizer &m) {
      m.setVerbose(false);
      m.setPrintLevel(-1);
      m.setStrategy(0);
   };

   RooMinimizer m1(*nllRooFit, RooMinimizer::FcnMode::classic);
   setupMinimizer(m1);
   RooMinimizer m2(*nllCodeGen, RooMinimizer::FcnMode::classic);
   setupMinimizer(m2);
   RooMinimizer m3(*nllCodeGen, RooMinimizer::FcnMode::clad_1);
   setupMinimizer(m3);
   RooMinimizer m4(*nllCodeGen, RooMinimizer::FcnMode::clad_2);
   setupMinimizer(m4);

#ifdef BENCH

   std::vector<RooMinimizer *> minimizers{&m1, &m2, &m3, &m4};
   std::size_t iMinimizer = state.range(0);
   auto &m = *minimizers[iMinimizer];

   // For validation
   resetParameters();
   m.minimize("");
   std::unique_ptr<RooFitResult> result{m.save()};
   result->Print();

   // If we don't use the actual RooFit model, we don't need to keep the dirty
   // flag propagation enabled.
   const bool disableDirtyFlagPropagation = &m != &m1;
   if (disableDirtyFlagPropagation) {
      RooAbsReal::setDirtyInhibit(true);
   }

   for (auto _ : state) {
      resetParameters();
      m.minimize("");
   }

   if (disableDirtyFlagPropagation) {
      RooAbsReal::setDirtyInhibit(false);
   }
#else
   auto &m = m4;
   m.minimize("");
   std::unique_ptr<RooFitResult> result{m.save()};
   result->Print();
#endif
}

#ifdef BENCH
auto unit = benchmark::kMicrosecond;

// For validating results
BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(0)->Iterations(1)->Name("RooFit_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(1)->Iterations(1)->Name("CodeGen_Numeric");
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(2)->Iterations(1)->Name("CodeGen_Clad_1");
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(3)->Iterations(1)->Name("CodeGen_Cald_2");

// For benchmarking
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(0)->Name("RooFit_Numeric");
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(1)->Name("CodeGen_Numeric");
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(2)->Name("CodeGen_Clad_1");
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(3)->Name("CodeGen_Cald_2");

// For profiling
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(3)->Iterations(100000)->Name("CodeGen_Cald_2");

// Define our main.
BENCHMARK_MAIN();
#endif