#include <RooAbsDataStore.h>
#include <RooAbsReal.h>
#include <RooDataSet.h>
#include <RooDerivative.h>
#include <RooGaussian.h>
#include <RooMinimizer.h>
#include <RooRealVar.h>
#include <RooVectorDataStore.h>
#include <RooFitResult.h>

#include <cmath>
#include <iostream>

#include <clad/Differentiator/Differentiator.h>

double gaus(double x, double mu, double sigma) {
  const double c = 1. / std::sqrt(2 * M_PI);
  const double arg = (x - mu) / sigma;
  return c / sigma * std::exp(-0.5 * arg * arg);
}

auto gaus_grad = clad::gradient(gaus, "mu,sigma");
auto gaus_dmu = clad::differentiate(gaus, "mu");
auto gaus_dsigma = clad::differentiate(gaus, "sigma");

class RooAutodiffGaussianNLL : public RooAbsReal {
public:
  RooAutodiffGaussianNLL(const char *name, const char *title, RooDataSet &data,
                         RooAbsReal &x, RooAbsReal &mean, RooAbsReal &sigma)
      : // Hack to get a reference to a vector with the generated x values
        _data{&static_cast<RooVectorDataStore *>(data.store())
                   ->realStoreList()[0]
                   ->data()},
        _mean("mean", "Mean", this, mean),
        _sigma("sigma", "Width", this, sigma) {}
  RooAutodiffGaussianNLL(const RooAutodiffGaussianNLL &other,
                         const char *name = nullptr)
      : _data{other._data}, _mean("mean", this, other._mean),
        _sigma("sigma", this, other._sigma) {}
  TObject *clone(const char *newname) const override {
    return new RooAutodiffGaussianNLL(*this, newname);
  }

  double evaluate() const override {
    double nll = 0.0;
    for (auto x : *_data) {
      double proba = gaus(x, _mean, _sigma);
      nll += -std::log(proba);
    }
    return nll;
  }

  void evaluateGradient(double *out) const override {
    std::cout << "evaluateGradient" << std::endl;
    double mean = _mean;
    double sigma = _sigma;

    if(sigma == 0.1) {
        out[0] = 0.0;
        out[1] = 0.0;
        return;
    }

    for (auto x : *_data) {
      double proba = gaus(x, _mean, _sigma);
      double res[2];
      res[0] = 0;
      res[1] = 0;
      gaus_grad.execute(x, mean, sigma, &res[0], &res[1]);
      //std::cout << res[0] << std::endl;
      out[0] += -1. / proba * res[0];
      out[1] += -1. / proba * res[1];
    }
    std::cout << mean << "   " << sigma << std::endl;
    std::cout << out[0] << "   " << out[1] << std::endl;
  }

  /// Return default level for MINUIT error analysis.
  double defaultErrorLevel() const override { return 0.5; }

private:
  std::vector<double> *_data = nullptr;
  RooRealProxy _mean;
  RooRealProxy _sigma;
};

void example_1() {

  // True values of parameters
  double mu_val_true = 0.0;
  double sigma_val_true = 1.0;

  // The values we use for evaluating the likelihood
  double mu_val = mu_val_true;
  double sigma_val = sigma_val_true;
  // double mu_val = 0.0268045;
  // double sigma_val = 1.03785;

  std::size_t nEvents = 1000;

  // Create the RooFit model
  RooRealVar rf_x{"x", "x", -10, 10};
  RooRealVar rf_mu{"mu", "mu", mu_val_true, -10, 10};
  RooRealVar rf_sigma{"sigma", "sigma", sigma_val_true, 0.1, 100};
  RooGaussian rf_gaus{"gaus", "gaus", rf_x, rf_mu, rf_sigma};

  // Create a RooFit dataset
  std::unique_ptr<RooDataSet> rf_data{rf_gaus.generate(rf_x, nEvents)};

  // The values of the parameters are changed to what we want to use for
  // evaluation
  rf_mu.setVal(mu_val);
  rf_sigma.setVal(sigma_val);

  //rf_sigma.setConstant();

  std::unique_ptr<RooAbsReal> rf_nll{rf_gaus.createNLL(*rf_data)};


  // Hack to get a reference to a vector with the generated x values
  auto const &data_x = static_cast<RooVectorDataStore *>(rf_data->store())
                           ->realStoreList()[0]
                           ->data();

  RooAutodiffGaussianNLL ad_nll{"ad_nll", "ad_nll", *rf_data,
                                rf_x,     rf_mu,    rf_sigma};

  double eps = 1e-4;
  std::unique_ptr<RooDerivative> ad_nll_grad_mu{
      ad_nll.derivative(rf_mu, 1, eps)};
  std::unique_ptr<RooDerivative> ad_nll_grad_sigma{
      ad_nll.derivative(rf_sigma, 1, eps)};

  std::cout << "RooFit land" << std::endl;
  std::cout << ad_nll.evaluate() << std::endl;
  std::cout << ad_nll_grad_mu->getVal() << std::endl;
  std::cout << ad_nll_grad_sigma->getVal() << std::endl;

  rf_mu.setVal(mu_val);
  rf_sigma.setVal(sigma_val);
  rf_mu.setError(0.0);
  rf_sigma.setError(0.0);

  std::cout << std::endl << "Standalone clad land" << std::endl;

  std::vector<double> nll_grad{};
  nll_grad.resize(2);

  ad_nll.evaluateGradient(nll_grad.data());

  std::cout << ad_nll.evaluate() << std::endl;
  std::cout << nll_grad[0] << std::endl;
  std::cout << nll_grad[1] << std::endl;

  //return;
  rf_mu.setVal(mu_val);
  rf_sigma.setVal(sigma_val);
  rf_mu.setError(0.0);
  rf_sigma.setError(0.0);

  RooMinimizer m(ad_nll, RooMinimizer::FcnMode::gradient);
  //m.setVerbose(10);
  //m.setPrintLevel(3);
  m.setStrategy(0);
  m.minimize("Minuit2");
  m.save()->Print();

  rf_mu.setVal(mu_val);
  rf_sigma.setVal(sigma_val);
  rf_mu.setError(0.0);
  rf_sigma.setError(0.0);

  RooMinimizer m2(ad_nll, RooMinimizer::FcnMode::clad);
  //m2.setVerbose(10);
  //m2.setPrintLevel(3);
  m2.setStrategy(0);
  m2.minimize("Minuit2");
  m2.save()->Print();
}
