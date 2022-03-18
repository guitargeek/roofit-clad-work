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


double nll001(double sigXsecOverSM)
{

   constexpr std::size_t nBins = 2;
   constexpr std::size_t nSamples = 3;

   double data[nBins]{122, 112};

   double scaleFactors[nBins * nSamples]{
      // sig,                  bkg1,       bkg2
      sigXsecOverSM, 1, 1, // bin 1
      sigXsecOverSM, 1, 1, // bin 2
   };

   double shapes[nBins * nSamples]{
      // sig,  bkg1,   bkg2
      20, 100, 0,   // bin 1
      10, 0,   100, // bin 2
   };

   double nll = 0.0;

   for (std::size_t iBin = 0; iBin < nBins; ++iBin) {

      double mu = 0.0;
      for (std::size_t iSample = 0; iSample < nSamples; ++iSample) {
         double scaleFactor = scaleFactors[nSamples * iBin + iSample];
         double shape = shapes[nSamples * iBin + iSample];
         mu += shape * scaleFactor;
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

auto nll001_gradient = clad::gradient(nll001, "sigXsecOverSM");

class HFLikelihoodWrapper final : public RooAbsReal {
public:
   HFLikelihoodWrapper(RooAbsReal &mu)
      : RooAbsReal{"RooFitDriverWrapper", "RooFitDriverWrapper"},
        _mu("mu", "mu", this, mu)
   {
   }
   HFLikelihoodWrapper(const HFLikelihoodWrapper &other, const char *name = 0)
      : RooAbsReal(other, name),
      _mu("sigma", this, other._mu)
   {
   }

   TObject *clone(const char *newname) const override { return new HFLikelihoodWrapper(*this, newname); }

   double defaultErrorLevel() const override { return 0.5; }

protected:
   double evaluate() const override { return nll001(_mu); }

   void evaluateGradient(double *out) const override
   {
      out[0] = 0;
      nll001_gradient.execute(_mu, &out[0]);
   }

private:
   RooRealProxy _mu;
};

void hf_example_02()
{

   using namespace RooStats;
   using namespace RooFit;
   using namespace HistFactory;

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
   meas.AddConstantParam("Lumi");

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
   signal.AddNormFactor("SigXsecOverSM", 1, 0, 3);
   chan.AddSample(signal);

   // Background 1
   Sample background1("background1", "background1", InputFile);
   chan.AddSample(background1);

   // Background 1
   Sample background2("background2", "background2", InputFile);
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

   pdf->graphVizTree("model.dot") ;

   auto * res = pdf->fitTo(*w->data("obsData"), Save(), PrintLevel(-1));
   res->Print();

   mu->setVal(1.0);
   mu->setError(0.0);

   HFLikelihoodWrapper wrapper(*mu);

   // Clad doesn't converge perfectly for some reason, but it works!
   RooMinimizer m(wrapper, RooMinimizer::FcnMode::clad_2);
   //RooMinimizer m(wrapper, RooMinimizer::FcnMode::classic);
   m.setPrintLevel(-1);
   m.minimize("Minuit2");
   m.save()->Print();
}
