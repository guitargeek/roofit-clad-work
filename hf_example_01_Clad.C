#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooRealProxy.h>
#include <RooAbsPdf.h>
#include <TROOT.h>
#include <TSystem.h>
#include <RooMinimizer.h>
#include <RooFitResult.h>
#include <RooStats/HistFactory/Measurement.h>
#include <RooStats/HistFactory/MakeModelAndMeasurementsFast.h>
#include <RooStats/ModelConfig.h>

#include <TMath.h>

#include <clad/Differentiator/Differentiator.h>

#include <chrono>

double nll002(double lumi, double sigXsecOverSM, double gamma_stat_channel1_bin_0, double gamma_stat_channel1_bin_1)
{

   constexpr std::size_t nBins = 2;
   constexpr std::size_t nSamples = 3;

   double taus[nBins]{
      400, 100, // TODO: explain what this is
   };
   double gammas[nBins]{gamma_stat_channel1_bin_0, gamma_stat_channel1_bin_1};
   double data[nBins]{122, 112};

   double scaleFactors[nBins * nSamples]{
      // sig,                  bkg1,       bkg2
      sigXsecOverSM, gammas[0], gammas[0], // bin 1
      sigXsecOverSM, gammas[1], gammas[1], // bin 2
   };

   double shapes[nBins * nSamples]{
      // sig,  bkg1,   bkg2
      20, 100, 0,   // bin 1
      10, 0,   100, // bin 2
   };

   double nll = 0.0;

   // lumi constraint
   {
      double arg = (lumi - 1.0) / 0.1;
      nll -= -0.5 * arg * arg; // Gaussian NLL without constant terms, where sigma is a constant
   }

   for (std::size_t iBin = 0; iBin < nBins; ++iBin) {

      // constraints
      {
         double mu = taus[iBin] * gammas[iBin];
         double n = taus[iBin];
         nll -= -mu + n * std::log(mu); // Poisson NLL without constant terms
      }

      double mu = 0.0;
      for (std::size_t iSample = 0; iSample < nSamples; ++iSample) {
         double scaleFactor = scaleFactors[nSamples * iBin + iSample];
         double shape = shapes[nSamples * iBin + iSample];
         mu += shape * scaleFactor * lumi;
      }

      double n = data[iBin];

      if (mu <= 0 && n > 0) {
         // Catch error condition: data present where zero events are
         // predicted
      } else if (std::abs(mu) < 1e-10 && std::abs(n) < 1e-10) {
         // Special handling of this case since log(Poisson(0,0)=0 but can't
         // be calculated with usual log-formula since log(mu)=0. No update
         // of result is required since term=0.
      } else {
         nll -= -mu + n * std::log(mu); // Poisson NLL without constant terms
      }
   }

   return nll;
}

auto nll002_gradient = clad::gradient(nll002, "lumi,sigXsecOverSM,gamma_stat_channel1_bin_0,gamma_stat_channel1_bin_1");

class HFLikelihoodWrapper final : public RooAbsReal {
public:
   HFLikelihoodWrapper(RooAbsReal &lumi, RooAbsReal &mu, RooAbsReal &gamma0, RooAbsReal &gamma1)
      : RooAbsReal{"RooFitDriverWrapper", "RooFitDriverWrapper"}, _lumi("lumi", "lumi", this, lumi),
        _mu("mu", "mu", this, mu), _gamma0("mu", "mu", this, gamma0), _gamma1("mu", "mu", this, gamma1)
   {
   }
   HFLikelihoodWrapper(const HFLikelihoodWrapper &other, const char *name = 0)
      : RooAbsReal(other, name), _lumi("lumi", this, other._lumi), _mu("sigma", this, other._mu),
        _gamma0("gamma0", this, other._gamma0), _gamma1("gamma1", this, other._gamma1)
   {
   }

   TObject *clone(const char *newname) const override { return new HFLikelihoodWrapper(*this, newname); }

   double defaultErrorLevel() const override { return 0.5; }

protected:
   double evaluate() const override { return nll002(_lumi, _mu, _gamma0, _gamma1); }

   void evaluateGradient(double *out) const override
   {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
      out[3] = 0;
      nll002_gradient.execute(_lumi, _mu, _gamma0, _gamma1, &out[0], &out[1], &out[2], &out[3]);
   }

private:
   RooRealProxy _lumi;
   RooRealProxy _mu;
   RooRealProxy _gamma0;
   RooRealProxy _gamma1;
};

void hf_example_01_Clad()
{

   using namespace RooStats;
   using namespace RooFit;
   using namespace HistFactory;

   gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");

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

   auto *w = MakeModelAndMeasurementFast(meas);

   auto *mc = static_cast<ModelConfig *>(w->obj("ModelConfig"));

   RooArgSet params;
   if (mc->GetParametersOfInterest())
      params.add(*mc->GetParametersOfInterest());
   if (mc->GetNuisanceParameters())
      params.add(*mc->GetNuisanceParameters());

   auto *pdf = w->pdf("simPdf");
   auto *mu = w->var("SigXsecOverSM");

   auto *res = pdf->fitTo(*w->data("obsData"), Constrain(params), GlobalObservables(*mc->GetGlobalObservables()),
                          Save(), PrintLevel(-10));
   // auto * res = pdf->fitTo(*w->data("obsData"), Save(), PrintLevel(-1));
   // res->Print();

   w->var("Lumi")->setVal(1.0);
   w->var("Lumi")->setError(0.0);
   mu->setVal(1.0);
   mu->setError(0.0);
   w->var("gamma_stat_channel1_bin_0")->setVal(1.0);
   w->var("gamma_stat_channel1_bin_0")->setError(0.0);
   w->var("gamma_stat_channel1_bin_1")->setVal(1.0);
   w->var("gamma_stat_channel1_bin_1")->setError(0.0);

   HFLikelihoodWrapper wrapper(*w->var("Lumi"), *mu, *w->var("gamma_stat_channel1_bin_0"),
                               *w->var("gamma_stat_channel1_bin_1"));

   // Clad doesn't converge perfectly for some reason, but it works!
   RooMinimizer m_clad(wrapper, RooMinimizer::FcnMode::clad_1);
   m_clad.setPrintLevel(10);
   m_clad.minimize("Minuit2");
}
