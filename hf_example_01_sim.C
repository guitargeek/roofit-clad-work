#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooRealProxy.h>
#include <RooAbsPdf.h>
#include <RooFunctor.h>
#include <TROOT.h>
#include <TSystem.h>
#include <RooMinimizer.h>
#include <RooFitResult.h>
#include <RooListProxy.h>
#include <RooStats/HistFactory/Measurement.h>
#include <RooStats/HistFactory/MakeModelAndMeasurementsFast.h>
#include <RooStats/ModelConfig.h>

#include <Minuit2/FCNGradientBase.h>
#include <Minuit2/VariableMetricMinimizer.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnStrategy.h>

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
   // meas.AddConstantParam("alpha_syst1");
   // meas.AddConstantParam("Lumi");

   meas.SetLumi(1.0);
   meas.SetLumiRelErr(0.10);
   meas.SetExportOnly(false);
   meas.SetBinHigh(2);

   for (std::size_t iChannel = 0; iChannel < nChannels; ++iChannel) {
      // Create a channel
      Channel chan(("channel" + std::to_string(iChannel)).c_str());

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
      // signal.AddOverallSys( "syst1",  0.95, 1.05 );
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
      // background1.AddOverallSys( "syst2", 0.95, 1.05  );
      chan.AddSample(background1);

      // Background 2
      Sample background2("background2");
      auto bkg2Hist = new TH1F{"bkg2_hist", "bkg2_hist", 2, 1.0, 2.0};
      bkg2Hist->SetBinContent(1, 0);
      bkg2Hist->SetBinContent(2, 100);
      background2.SetHisto(bkg2Hist);
      background2.ActivateStatError();
      // background2.AddOverallSys( "syst3", 0.95, 1.05  );
      chan.AddSample(background2);

      meas.AddChannel(chan);
   }

   return std::unique_ptr<RooWorkspace>{MakeModelAndMeasurementFast(meas)};
}

std::string generateNLLCode(contextManager &ctx)
{

   ExRooRealVar X(ctx, "x", "{1.25, 1.75}", 2);
   ExRooConst B2Eps(ctx, 1);
   ExRooConst B1Eps(ctx, 1);
   ExRooRealVar NomLumi(ctx, "nominalLumi", "1");
   ExRooConst RandConst(ctx, 0.1);
   // ExRooRealVar AlphaSys(ctx, "alphaSys", "0");
   ExRooConst RandConst2(ctx, 1);
   // ExRooRealVar NomAlphaSys(ctx, "nomAlphaSys", "0");
   ExRooRealVar GammaB1(ctx, "gamma_stat_channel0_bin_0"); // Input
   ExRooRealVar NomGammaB1(ctx, "nomGammaB1", "400");
   ExRooRealVar GammaB2(ctx, "gamma_stat_channel0_bin_1"); // Input
   ExRooRealVar SigXOverM(ctx, "SigXsecOverSM");           // Input
   ExRooRealVar Lumi(ctx, "Lumi");                         // Input
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
   return "double nll(double* in) { \n" + code + " return " + Root.getResult() + ";\n}\n";
   // return "double nll(" + ctx.getParamList() + ") { \n" + code + " return " + Root.getResult() + ";\n}\n";
}

// Generated Code
double nll(double *in)
{
   double nomGammaB1 = 400;
   double nomGammaB2 = 100;
   double nominalLumi = 1;
   double constraint[3]{ExRooPoisson::poisson(nomGammaB1, (nomGammaB1 * in[0])),
                        ExRooPoisson::poisson(nomGammaB2, (nomGammaB2 * in[1])),
                        ExRooGaussian::gauss(in[2], nominalLumi, 0.100000)};
   double cnstSum = 0;
   double x[2]{1.25, 1.75};
   double sig[2]{20, 10};
   double binBoundaries1[3]{1, 1.5, 2};
   double bgk1[2]{100, 0};
   double binBoundaries2[3]{1, 1.5, 2};
   double histVals[2]{in[0], in[1]};
   double bgk2[2]{0, 100};
   double binBoundaries3[3]{1, 1.5, 2};
   double weights[2]{122.000000, 112.000000};
   for (int i = 0; i < 3; i++) {
      cnstSum -= std::log(constraint[i]);
   }
   double mu = 0;
   double temp;
   double nllSum = 0;
   unsigned int b1, b2, b3;
   for (int iB = 0; iB < 2; iB++) {
      b1 = ExRooHistFunc::getBin(binBoundaries1, x[iB]);
      b2 = ExRooHistFunc::getBin(binBoundaries2, x[iB]);
      b3 = ExRooHistFunc::getBin(binBoundaries3, x[iB]);
      mu = 0;
      mu += sig[b1] * (in[3] * in[2]);
      mu += (bgk1[b2] * histVals[iB]) * (in[2] * 1.000000);
      mu += (bgk2[b3] * histVals[iB]) * (in[2] * 1.000000);
      temp = std::log((mu));
      nllSum -= -(mu) + weights[iB] * temp;
   }
   return cnstSum + nllSum;
}

template <typename Func = void>
class MinuitFuncWrapper final : public ROOT::Minuit2::FCNBase {
public:
   MinuitFuncWrapper(Func function) : _func(function) {}

   double Up() const override { return 0.5; }

   double operator()(const std::vector<double> &x) const override { return _func(const_cast<double *>(x.data())); }

private:
   Func _func;
};

template <typename Func = void, typename Grad = void>
class MinuitGradFuncWrapper final : public ROOT::Minuit2::FCNGradientBase {
public:
   MinuitGradFuncWrapper(Func function, Grad gradient, std::size_t n) : _func(function), _grad(gradient), _n{n} {}

   double Up() const override { return 0.5; }

   double operator()(const std::vector<double> &x) const override { return _func(const_cast<double *>(x.data())); }

   std::vector<double> Gradient(const std::vector<double> &x) const override
   {
      std::vector<double> out(_n);
      clad::array_ref<double> cladout(out.data(), _n);
      _grad(const_cast<double *>(x.data()), cladout);
      return out;
   }

private:
   Func _func;
   Grad _grad;
   std::size_t _n = 0;
};

constexpr bool verbose = true;

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

   contextManager ctx;
   std::string func = generateNLLCode(ctx);
   // std::cout << func.c_str();

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

   std::unique_ptr<RooWorkspace> w = makeHistFactoryWorkspace(1);

   auto *mc = static_cast<ModelConfig *>(w->obj("ModelConfig"));

   RooArgSet params;
   if (mc->GetParametersOfInterest())
      params.add(*mc->GetParametersOfInterest());
   if (mc->GetNuisanceParameters())
      params.add(*mc->GetNuisanceParameters());

   std::unique_ptr<RooArgSet> originalParams{params.snapshot()};

   auto resetParameters = [&]() {
      for (auto *param : static_range_cast<RooRealVar *>(params)) {
         param->setError(0.0);
      }
   };

   auto *pdf = w->pdf("simPdf");

   // Create the likelihood objects with RooFit
   std::unique_ptr<RooAbsReal> nllRooFit{pdf->createNLL(*w->data("obsData"), Constrain(*mc->GetNuisanceParameters()),
                                                        GlobalObservables(*mc->GetGlobalObservables()))};
   RooFunctor nllFunctor{*nllRooFit, {}, params, *w->data("obsData")->get()};
   std::unique_ptr<RooAbsReal> nllRooFitBatchMode{
      pdf->createNLL(*w->data("obsData"), Constrain(*mc->GetNuisanceParameters()),
                     GlobalObservables(*mc->GetGlobalObservables()), BatchMode(true))};
   RooFunctor nllFunctorBatchMode{*nllRooFitBatchMode, {}, params, *w->data("obsData")->get()};

   MinuitFuncWrapper minuitRooFunc(nllFunctor);
   MinuitFuncWrapper minuitRooFuncBatchMode(nllFunctorBatchMode);
   MinuitFuncWrapper minuitFunc(nll);
   MinuitGradFuncWrapper minuitGradFunc(nll, nll_grad, params.size());

   // We set all parameters away from the minimum as the initial state of the fit, such that the fit is not trivial
   // TODO: also set limits of parameters
   ROOT::Minuit2::MnUserParameters mnParamsRooFit;
   ROOT::Minuit2::MnUserParameters mnParamsCodeGen;
   for (std::size_t i = 0; i < params.size(); ++i) {
      mnParamsRooFit.Add(params[i]->GetName(), 0.9, 0.0);
      mnParamsCodeGen.Add(ctx.inputParams[i], 0.9, 0.0);
   }
   ROOT::Minuit2::MnStrategy mnStrategy{0};

   auto reorderFitParams = [&](ROOT::Minuit2::MnUserParameters const &fitPars) {
      ROOT::Minuit2::MnUserParameters orderedPars;

      for (std::size_t i = 0; i < params.size(); ++i) {
         const char *name = mnParamsCodeGen.GetName(i).c_str();
         orderedPars.Add(name, fitPars.Value(name), fitPars.Error(name));
      }

      return orderedPars;
   };

#ifdef BENCH

   std::size_t iBenchmark = state.range(0);

   if (iBenchmark == 0) {
      if (verbose) {
         auto minimum = ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitRooFunc, mnParamsRooFit, mnStrategy);
         std::cout << reorderFitParams(minimum.UserParameters()) << std::endl;
      }

      for (auto _ : state) {
         ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitRooFunc, mnParamsRooFit, mnStrategy);
      }
   } else if (iBenchmark == 1) {
      if (verbose) {
         auto minimum =
            ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitRooFuncBatchMode, mnParamsRooFit, mnStrategy);
         std::cout << reorderFitParams(minimum.UserParameters()) << std::endl;
      }

      for (auto _ : state) {
         ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitRooFuncBatchMode, mnParamsRooFit, mnStrategy);
      }
   } else if (iBenchmark == 2) {
      if (verbose) {
         auto minimum = ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitFunc, mnParamsCodeGen, mnStrategy);
         std::cout << minimum.UserParameters() << std::endl;
      }

      for (auto _ : state) {
         ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitFunc, mnParamsCodeGen, mnStrategy);
      }
   } else if (iBenchmark == 3) {
      if (verbose) {
         auto minimum = ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitGradFunc, mnParamsCodeGen, mnStrategy);
         std::cout << minimum.UserParameters() << std::endl;
      }

      for (auto _ : state) {
         ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitGradFunc, mnParamsCodeGen, mnStrategy);
      }
   }

#else
#endif
}

#ifdef BENCH
auto unit = benchmark::kMicrosecond;

const auto nIter = verbose ? 1 : 0;

BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(0)->Iterations(nIter)->Name("RooFit_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(1)->Iterations(nIter)->Name("BatchMode_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(2)->Iterations(nIter)->Name("CodeGen_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(3)->Iterations(nIter)->Name("CodeGen_Clad");

// For profiling
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(3)->Iterations(100000)->Name("CodeGen_Cald");

// Define our main.
BENCHMARK_MAIN();
#endif
