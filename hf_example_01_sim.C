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
#include <string>
#include <cmath>

struct HFData {
   unsigned int bins;
   unsigned int channels;
   std::vector<double> weights;
   std::vector<double> binVals;
   std::vector<double> sig;
   std::vector<double> bkg1;
   std::vector<double> bkg2;
   std::vector<double> nomGammaVals;
   std::vector<double> binBoundaries;

   static std::string toString(const std::vector<double>& in) {
     std::string vec = "{";
     for (auto &it : in) {
       vec += std::to_string(it) + ",";
     }
     vec.pop_back();
     vec += "};";
     return vec;
   }

   void fillData(unsigned int numChannels, unsigned int numBins){

   }
};

std::unique_ptr<RooWorkspace> makeHistFactoryWorkspace(HFData& data)
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

   for (std::size_t iChannel = 0; iChannel < data.channels; ++iChannel) {
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

std::string generateNLLCode(contextManager &ctx, HFData& data)
{

   // ---------- Constants ----------
   ExRooRealVar X(ctx, "x", HFData::toString(data.binVals), 2);
   ExRooConst B2Eps(ctx, 1);
   ExRooConst B1Eps(ctx, 1);
   ExRooRealVar NomLumi(ctx, "nominalLumi", "1");
   ExRooConst RandConst(ctx, 0.1);
   // ExRooRealVar AlphaSys(ctx, "alphaSys", "0");
   ExRooConst RandConst2(ctx, 1);
   // ExRooRealVar NomAlphaSys(ctx, "nomAlphaSys", "0");

   // ---------- Channel Independent Vars ----------
   ExRooRealVar Lumi(ctx, "Lumi");                         // Input
   ExRooRealVar SigXOverM(ctx, "SigXsecOverSM");           // Input
   // Channel weight
   ExRooRealVar ChannelWeight(ctx, "cw", "1");
   ExRooRealVar ChannelWeight2(ctx, "cw2", "1");
   ExRooGaussian GaussLumi(ctx, &Lumi, &NomLumi, &RandConst);

   // ---------- Channel Dependent Vars ----------
   std::vector<ExRooReal*> channelVars;
   std::vector<ExRooReal*> nll;
   std::vector<ExRooReal*> constraints;
   constraints.push_back(&GaussLumi);
   std::string binBoundaries = HFData::toString(data.binBoundaries);
   std::string sigVals = HFData::toString(data.binBoundaries);
   std::string bkg1Vals = HFData::toString(data.binBoundaries);
   std::string bkg2Vals = HFData::toString(data.binBoundaries);
   for (int i = 0; i < data.channels; i++) {
     std::string ichan = std::to_string(i);

     // ---------- bin Dependent Vars ----------
     std::vector<ExRooReal *> gammas;
     for (int j = 0; j < data.bins; j++) {
       std::string ibin = std::to_string(j);

       // Gammas
       channelVars.push_back(new ExRooRealVar(
           ctx, "gamma_stat_channel" + ichan + "_bin_" + ibin)); // Input
       gammas.push_back(channelVars.back());
       int inomgamma = channelVars.size();
       channelVars.push_back(
           new ExRooRealVar(ctx, "nom_gamma_chan" + ichan + "_bin_" + ibin,
                            std::to_string(data.nomGammaVals[j])));

       // Constraints
       channelVars.push_back(new ExRooProduct(
           ctx, {channelVars[inomgamma], gammas.back()}));
       channelVars.push_back(
           new ExRooPoisson(ctx, channelVars[inomgamma], channelVars.back()));
       constraints.push_back(channelVars.back());
       // channelVars.push_back(new ExRooGaussian(ctx, &AlphaSys, &NomAlphaSys,
       // &RandConst2)); 
       //constraints.push_back(channelVars.back());
     }

     // HistFunc components
     int isig = channelVars.size();
     channelVars.push_back(new ExRooHistFunc(ctx, &X, "sig" + ichan, sigVals, binBoundaries));
     int ibgk1 = channelVars.size();
     channelVars.push_back(new ExRooHistFunc(ctx, &X, "bgk1" + ichan, bkg1Vals, binBoundaries));
     int ibgk2 = channelVars.size();
     channelVars.push_back(new ExRooHistFunc(ctx, &X, "bgk2" + ichan, bkg2Vals, binBoundaries));
     int iparamHist = channelVars.size();
     channelVars.push_back(new ExRooParamHistFunc(ctx, &X, gammas, "histVals" + ichan));
     int iscale1 = channelVars.size();
     channelVars.push_back(new ExRooProduct(ctx, {/* &AlphaSys, */ &SigXOverM, &Lumi}));
     int iscale2 = channelVars.size();
     channelVars.push_back(new ExRooProduct(ctx, {&Lumi, &B1Eps}));
     int iscale3 = channelVars.size();
     channelVars.push_back(new ExRooProduct(ctx, {&Lumi, &B2Eps}));
     int ibkgShape1 = channelVars.size();
     channelVars.push_back(new ExRooProduct(ctx, {channelVars[ibgk1], channelVars[iparamHist]}));
     int ibkgShape2 = channelVars.size();
     channelVars.push_back(new ExRooProduct(ctx, {channelVars[ibgk2], channelVars[iparamHist]}));
     channelVars.push_back(new ExRooRealSum( 
         ctx, "mu" + ichan,
         {channelVars[isig], channelVars[ibkgShape1], channelVars[ibkgShape2]},
         {channelVars[iscale1], channelVars[iscale2], channelVars[iscale3]}));
     channelVars.push_back(new ExRooProduct(ctx, {channelVars.back()}));
     channelVars.push_back(new ExRooNll2(ctx, &X, std::to_string(data.bins), channelVars.back(), data.weights, "nllSum" + ichan, "weights" + ichan));

     // NLL
     nll.push_back(channelVars.back());
   }

   // Constraint sum
   ExRooConstraintSum ConstrSum(ctx, constraints);
   nll.push_back(&ConstrSum);

   // Final root model node
   ExRooAddition Root(ctx, nll);

   std::string code = Root.getCode();
   std::string retVal = Root.getResult();

   // ---------- Cleanup ----------
   for(auto it : channelVars)
     delete it;

   return "double " + ctx.funcName + "(double* in) { \n" + code + " return " + retVal + ";\n}\n";
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

constexpr bool verbose = false;
std::unordered_map<int,bool> nllDeclared;

#ifdef BENCH
static void hf_example_01_sim(benchmark::State &state)
#else
int main()
#endif
{
   using namespace RooStats;
   using namespace RooFit;
#ifdef BENCH
   int testChannels = state.range(1);
#else
   int testChannels = 50;
   int testBins = 2;
#endif

   gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
   auto &msg = RooMsgService::instance();
   msg.setGlobalKillBelow(RooFit::WARNING);

   gInterpreter->Declare("#pragma cling optimize(2)");

   contextManager ctx("nll" + std::to_string(testChannels));
   HFData hf_data;

   // Fill HF data 
   hf_data.fillData(testChannels, testBins);

   std::string funcName = ctx.funcName;
   std::string gradName = funcName + "_grad";
   std::string func = generateNLLCode(ctx, hf_data);

   if (!nllDeclared[testChannels]) {
      gInterpreter->ProcessLine(
         "#include \"/home/grimmyshini/ROOT/examples/roofit-clad-work/include/RooSimClasses.h\"");
      gInterpreter->Declare(func.c_str());

      gInterpreter->ProcessLine("#include \"clad/Differentiator/Differentiator.h\"");
      // calculate gradient
      gInterpreter->ProcessLine(("clad::gradient(" + funcName + ");").c_str());
      nllDeclared[testChannels] = true;
   }

   // get the grad function pointer.
   auto gradObj =
       (void (*)(double *, clad::array_ref<double>))gInterpreter->ProcessLine(
           (gradName + ";").c_str());
   auto funcObj = (double (*)(double *))gInterpreter->ProcessLine(
       (funcName + ";").c_str());

   std::unique_ptr<RooWorkspace> w = makeHistFactoryWorkspace(hf_data);

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
   MinuitFuncWrapper minuitFunc(funcObj);
   MinuitGradFuncWrapper minuitGradFunc(funcObj, gradObj, params.size());

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
   auto minimum = ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitRooFunc, mnParamsRooFit, mnStrategy);
   std::cout << reorderFitParams(minimum.UserParameters()) << std::endl;

   auto minimum2 = ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitFunc, mnParamsCodeGen, mnStrategy);
   std::cout << minimum2.UserParameters() << std::endl;   

   auto minimum3 = ROOT::Minuit2::VariableMetricMinimizer{}.Minimize(minuitGradFunc, mnParamsCodeGen, mnStrategy);
   std::cout << minimum3.UserParameters() << std::endl;
#endif
}

#ifdef BENCH
auto unit = benchmark::kMicrosecond;

const auto nIter = 10;

BENCHMARK(hf_example_01_sim)->Unit(unit)->ArgsProduct({{{0},benchmark::CreateDenseRange(1, 100, /*step=*/10)}})->Iterations(nIter)->Name("RooFit_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->ArgsProduct({{{1},benchmark::CreateDenseRange(1, 100, /*step=*/10)}})->Iterations(nIter)->Name("BatchMode_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->ArgsProduct({{{2},benchmark::CreateDenseRange(1, 100, /*step=*/10)}})->Iterations(nIter)->Name("CodeGen_Numeric");
BENCHMARK(hf_example_01_sim)->Unit(unit)->ArgsProduct({{{3},benchmark::CreateDenseRange(1, 100, /*step=*/10)}})->Iterations(nIter)->Name("CodeGen_Clad");

// For profiling
// BENCHMARK(hf_example_01_sim)->Unit(unit)->Arg(3)->Iterations(1000)->Name("CodeGen_Cald");

// Define our main.
BENCHMARK_MAIN();
#endif
