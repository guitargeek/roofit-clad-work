// can be compiled with the root-compile script in the same repository

#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooExponential.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooDerivative.h>

#include <clad/Differentiator/Differentiator.h>

namespace calc {

double gaus(double x, double mu, double sigma)
{
   const double c = 1. / std::sqrt(2 * M_PI);
   const double arg = (x - mu) / sigma;
   return c / sigma * std::exp(-0.5 * arg * arg);
}

double expo(double x, double c)
{
   return -c * std::exp(c * x);
}

double add(double a, double b, double f)
{
   return f * a + (1 - f) * b;
}

auto gaus_clad_grad = clad::gradient(gaus, "x,mu,sigma");
auto expo_clad_grad = clad::gradient(expo, "x,c");
auto add_clad_grad = clad::gradient(add, "a,b,f");

} // namespace calc

template <std::size_t... Ints>
using Indices = std::index_sequence<Ints...>;

struct ChainRuleDriverSimple {
public:
   ChainRuleDriverSimple(std::size_t nVariables)
   {
      _cladBuffer.resize(nVariables);
      values.resize(nVariables);
      gradients.resize(nVariables);
      for (auto &v : gradients) {
         v.resize(nVariables);
      }
      // we can initialize the gradients as a diagonal matrix, as dx/dx = 1
      for (std::size_t i = 0; i < gradients.size(); ++i) gradients[i][i] = 1.0;
   }

   template <typename FuncType, typename GradFuncType, size_t... VarIndices>
   void eval(std::size_t iVar, FuncType &f, GradFuncType &g, Indices<VarIndices...> varIndices)
   {
      values[iVar] = callFunc(f, varIndices);
      callGradientFunc(g, varIndices, std::make_index_sequence<varIndices.size()>{});
      getGradientFromBuffer(iVar, varIndices);
      applyChainRule(iVar, varIndices);
   }

   std::vector<double> values;
   std::vector<std::vector<double>> gradients;

private:
   template <typename FuncType, size_t... VarIndices>
   double callFunc(FuncType &f, Indices<VarIndices...>)
   {
      return f(values[VarIndices]...);
   }

   template <typename GradFuncType, size_t... VarIndices, size_t... I>
   void callGradientFunc(GradFuncType &g, Indices<VarIndices...>, Indices<I...>)
   {
      for (auto &x : _cladBuffer) {
         x = 0.0;
      }
      g.execute(values[VarIndices]..., &_cladBuffer[I]...);
   }

   template <size_t... VarIndices>
   void getGradientFromBuffer(std::size_t iVar, Indices<VarIndices...>)
   {
      std::size_t i = 0;
      for (std::size_t j : {VarIndices...}) {
         gradients[iVar][j] = _cladBuffer[i];
         ++i;
      }
   }

   template <size_t... VarIndices>
   void applyChainRule(std::size_t iVar, Indices<VarIndices...>)
   {
      // loop over direct variables
      for (std::size_t j : {VarIndices...}) {
         // loop over indirect variables
         for (std::size_t k = 0; k < values.size(); ++k) {
            // make sure to not double count the gradients of direct variables
            if(k != j) {
                gradients[iVar][k] += gradients[iVar][j] * gradients[j][k];
            }
         }
      }
   }

   std::vector<double> _cladBuffer;
};

void chain_rule_example_01()
{

   /////////////////
   // Free functions
   {
      std::cout << "Free functions\n";
      std::cout << "==============\n";

      // Variables
      double x = 1.0;
      double mu = 5.0;
      double sigma = 1.0;
      double c = -0.5;
      double f = 0.2;

      ChainRuleDriverSimple driver{8};

      auto &values = driver.values;
      values[0] = x;
      values[1] = mu;
      values[2] = sigma;
      values[3] = c;
      values[4] = f;
      values[5] = 0.0; // intermediate result from Gaussian
      values[6] = 0.0; // intermediate result from Exponential
      values[7] = 0.0; // final result from the addition

      using namespace calc;

      // The ChainRuleDriverSimple makes it possible to 

      driver.eval(5, gaus, gaus_clad_grad, Indices<0, 1, 2>{});
      driver.eval(6, expo, expo_clad_grad, Indices<0, 3>{});
      driver.eval(7, add, add_clad_grad, Indices<5, 6, 4>{});

      for (auto x : values)
         std::cout << x << "  ";
      std::cout << std::endl;
      std::cout << std::endl;

      for(std::size_t i = 0; i < 8; ++i) {
          for (auto x : driver.gradients[i])
             std::cout << x << "  ";
          std::cout << std::endl;
      }
   }

   ////////////////
   // RooFit model
   {
      std::cout << "RooFit model\n";
      std::cout << "============\n";

      // Observables
      RooRealVar x{"x", "x", 1.0, 0.0, 100};

      // Parameters
      RooRealVar mu{"mu", "mu", 5, 0.0, 10};
      RooRealVar sigma{"sigma", "sigma", 1.0, 0.1, 10};
      RooRealVar c{"c", "c", -0.5, -10, -0.1};
      RooRealVar f{"c", "c", 0.2, 0.0, 1.0};

      // PDFs
      RooGaussian gaus{"gaus", "gaus", x, mu, sigma};
      RooExponential expo{"expo", "expo", x, c};
      RooAddPdf model{"model", "model", {gaus, expo}, {f}};

      RooArgSet normSet{x};

      std::cout << model.getVal(normSet) << std::endl;
      std::cout << std::endl;

      using AbsRealPtr = std::unique_ptr<RooAbsReal>;
      AbsRealPtr{model.derivative(mu, normSet, 1)}->Print();
      AbsRealPtr{model.derivative(sigma, normSet, 1)}->Print();
      AbsRealPtr{model.derivative(c, normSet, 1)}->Print();
      AbsRealPtr{model.derivative(f, normSet, 1)}->Print();

      std::cout << std::endl;

      // For visualization of the computation graph,
      // create gif file with `dot -Tgif -o model.gif model.dot`
      model.graphVizTree("model.dot");
   }
}
