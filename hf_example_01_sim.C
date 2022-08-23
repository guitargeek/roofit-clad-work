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
#include <string>
#include <cmath>

std::unique_ptr<RooWorkspace> makeHistFactoryWorkspace()
{

   std::string InputFile = "./data/example.root";
   // in case the file is not found
   bool bfile = gSystem->AccessPathName(InputFile.c_str());
   if (bfile) {
      std::cout << "Input file is not found - run prepareHistFactory script " << std::endl;
      gROOT->ProcessLine(".! prepareHistFactory .");
      bfile = gSystem->AccessPathName(InputFile.c_str());
      if (bfile) {
         std::cout << "Still no " << InputFile << ", giving up.\n";
         exit(1);
      }
   }

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

   // Create a channel

   Channel chan("channel1");
   chan.SetData("data", InputFile);
   chan.SetStatErrorConfig(0.05, "Poisson");

   // Now, create some samples

   // Create the signal sample
   Sample signal("signal", "signal", InputFile);
   // signal.AddOverallSys( "syst1",  0.95, 1.05 );
   signal.AddNormFactor("SigXsecOverSM", 1, 0, 3);
   chan.AddSample(signal);

   // Background 1
   Sample background1("background1", "background1", InputFile);
   background1.ActivateStatError("background1_statUncert", InputFile);
   // background1.AddOverallSys( "syst2", 0.95, 1.05  );
   chan.AddSample(background1);

   // Background 1
   Sample background2("background2", "background2", InputFile);
   background2.ActivateStatError();
   // background2.AddOverallSys( "syst3", 0.95, 1.05  );
   chan.AddSample(background2);

   meas.AddChannel(chan);
   meas.CollectHistograms();

   return std::unique_ptr<RooWorkspace>{MakeModelAndMeasurementFast(meas)};
}

using namespace std;

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

// Generated Code
double nll(double lumi, double SigXsecOverSM, double gammaB1, double gammaB2)
{
   double nomGammaB1 = 100;
   double nomGammaB2 = 400;
   double nominalLumi = 1;
   double constraint[3]{ExRooPoisson::poisson(nomGammaB1, (nomGammaB1 * gammaB1)),
                        ExRooPoisson::poisson(nomGammaB2, (nomGammaB2 * gammaB2)),
                        ExRooGaussian::gauss(lumi, nominalLumi, 0.100000)};
   double cnstSum = 0;
   double x[2]{1.25, 1.75};
   double sig[2]{20, 10};
   double binBoundaries1[3]{1, 1.5, 2};
   double bgk1[2]{100, 0};
   double binBoundaries2[3]{1, 1.5, 2};
   double histVals[2]{gammaB1, gammaB2};
   double bgk2[2]{0, 100};
   double binBoundaries3[3]{1, 1.5, 2};
   double weights[2]{122.000000, 112.000000};
   for (int i = 0; i < 3; i++) {
      cnstSum -= std::log(constraint[i]);
   }
   double nllSum = 0;
   for (int iB = 0; iB < 2; iB++) {
      unsigned int b1 = ExRooHistFunc::getBin(binBoundaries1, x[iB]);
      unsigned int b2 = ExRooHistFunc::getBin(binBoundaries2, x[iB]);
      unsigned int b3 = ExRooHistFunc::getBin(binBoundaries3, x[iB]);
      double mu = 0;
      mu += sig[b1] * (SigXsecOverSM * lumi);
      mu += (bgk1[b2] * histVals[iB]) * (lumi * 1.000000);
      mu += (bgk2[b3] * histVals[iB]) * (lumi * 1.000000);
      double temp;
      temp = std::log((mu));
      nllSum -= -(mu) + weights[iB] * temp;
   }
   return cnstSum + nllSum;
}

void nll_grad(double lumi, double SigXsecOverSM, double gammaB1, double gammaB2, clad::array_ref<double> _d_lumi,
              clad::array_ref<double> _d_SigXsecOverSM, clad::array_ref<double> _d_gammaB1,
              clad::array_ref<double> _d_gammaB2)
{
}

template <typename Func = void, typename Grad = void>
class RooGradFuncWrapper final : public RooAbsReal {
public:
   RooGradFuncWrapper(Func function, Grad gradient, RooWorkspace const &ws, std::vector<std::string> const &paramNames)
      : RooAbsReal{"RooGradFuncWrapper", "RooGradFuncWrapper"}, _params{"!params", "List of parameters", this},
        _func(function), _grad(gradient)
   {
      for (auto const &paramName : paramNames) {
         _params.add(*ws.var(paramName.c_str()));
      }
   }
   RooGradFuncWrapper(const RooGradFuncWrapper &other, const char *name = nullptr)
      : RooAbsReal(other, name), _params("!params", this, other._params), _func(other._func), _grad(other._grad)
   {
   }

   TObject *clone(const char *newname) const override { return new RooGradFuncWrapper(*this, newname); }

   double defaultErrorLevel() const override { return 0.5; }

protected:
   double evaluate() const override
   {
      return _func(
         static_cast<RooAbsReal const &>(_params[0]).getVal(), static_cast<RooAbsReal const &>(_params[1]).getVal(),
         static_cast<RooAbsReal const &>(_params[2]).getVal(), static_cast<RooAbsReal const &>(_params[3]).getVal());
   }

   void evaluateGradient(double *out) const override
   {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
      out[3] = 0;
      _grad(static_cast<RooAbsReal const &>(_params[0]).getVal(), static_cast<RooAbsReal const &>(_params[1]).getVal(),
            static_cast<RooAbsReal const &>(_params[2]).getVal(), static_cast<RooAbsReal const &>(_params[3]).getVal(),
            &out[0], &out[1], &out[2], &out[3]);
   }

private:
   RooListProxy _params;
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

   std::string code = Root.getCode();
   //   std::string func = "double nll(" + ctx.getParamList() + ") { \n" + code + " return " + Root.getResult() + ";
   //   \n}\n"; std::cout << func.c_str();

   //   auto df = clad::gradient(nll);
   //   df.dump();

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

   std::unique_ptr<RooWorkspace> w = makeHistFactoryWorkspace();

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

   std::unique_ptr<RooFitResult> res{pdf->fitTo(
      *w->data("obsData"), Constrain(params), GlobalObservables(*mc->GetGlobalObservables()), Save(), PrintLevel(-1))};
   // auto * res = pdf->fitTo(*w->data("obsData"), Save(), PrintLevel(-1));
   res->Print();

   resetParameters();

   RooGradFuncWrapper wrapper(nll, nll_grad, *w,
                              {"Lumi", "SigXsecOverSM", "gamma_stat_channel1_bin_0", "gamma_stat_channel1_bin_1"});

   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // RooMinimizer m(wrapper, RooMinimizer::FcnMode::clad_2); //720 us          720 us
   // RooMinimizer m(wrapper, RooMinimizer::FcnMode::clad_1);
   // RooMinimizer m(wrapper, RooMinimizer::FcnMode::classic); //713 us          709 us
   std::unique_ptr<RooAbsReal> nll_1{pdf->createNLL(*w->data("obsData"))}; // 1653 us         1654 us
   RooMinimizer m(*nll_1, RooMinimizer::FcnMode::classic);

   m.setVerbose(false);
   m.setPrintLevel(-1);

#ifdef BENCH
   for (auto _ : state) {
      resetParameters();
      m.minimize("");
   }
#else
   m.minimize("");
#endif

   std::unique_ptr<RooFitResult> result{m.save()};
   result->Print();
}

#ifdef BENCH
auto unit = benchmark::kMicrosecond;
BENCHMARK(hf_example_01_sim)->Unit(unit);
// Define our main.
BENCHMARK_MAIN();
#endif
