# RooFit AD Master Student Project

The full code for this introduction can be found in the [chain_rule_example_01.cpp](https://github.com/guitargeek/roofit-clad-work/blob/main/chain_rule_example_01.cpp).

## Introduction

RooFit allows you to build statistical models in form of **probability distribution functions (PDFs)**. These PDFs depend in general on many unknown parameters. The optimum parameters are found by minimizing the likelihood (a function of the PDF and a dataset), which is called "parameter fitting". To find the minimum it's important to know the gradient of the PDF.

Currently, the gradient is found numerically, which is expensive in the case of many parameters because the parameters are varied one-by-one and the PDF is re-evaluated each time. That's why we want to implemet functionality to evaluate also the gradient in one pass through the PDF's computation graph, which consists of multiple objects inheriting from `RooAbsArg` in a tree-like structure that is evaluated bottom to top.

## The RooFit computation graph

As explained on the [ROOT website](https://root.cern/manual/roofit/), a RooFit model is specified as a tree of `RooAbsArgs` that represents the mathematical pdf.

For example, we can build a PDF that is the sum of a Gaussian and an exponential probability density in one dimension:
$$
\text{model}(x|\mu,\sigma,c,f) = f \cdot\text{gaus}(x|\mu,\sigma) + (1-f) \cdot\text{expo}(x|c).
$$
Here, $x$ is a fixed observed value (*observable*), and $\mu,\sigma,c,f$ are the PDFs parameters.

RooFit doesn't make a difference between observables and parameters. In general they are called **variables** and repreesnted by the `RooRealVar` class. In RooFit, the model above can be built as follows:
```C++
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
```
To evaluate the model in a meaningful way, we also need to create a set of **normalization variables**. Any PDF is normalized to one when integrating over the observables by definition, but since the model doesn't know what are parameters and observables we need to pass the observables as a normalization set for the PDF to get normalized correctly.
```C++
      RooArgSet normSet{x}; // create the normalization set
      std::cout << model.getVal(normSet) << std::endl; // evaluate the PDF
```
The output will be:
```
0.242639
```
It is also useful to export the model to a dot file to visualize the ocmputation graph:
```C++
      // For visualization of the computation graph,
      // create gif file with `dot -Tgif -o model.gif model.dot`
      model.graphVizTree("model.dot");
```
The resulting GIF file looks like this (fundamental variables are blue and derived functions or PDFs are red):
![](https://codimd.web.cern.ch/uploads/upload_0f37496fc941b9eb49cf4a5985989336.gif)

## Numeric differentiation

Any RooFit function or PDF can be numerically differentiated using [RooAbsReal::derivative()](https://root.cern.ch/doc/master/classRooAbsReal.html#a912f32dc088deafab33317c9e87d3698), which is useful for validating our (automatically generated) analytical gradients in this project.

We can use it to get the full gradient, i.e., the vector of derivatives with respect to each parameter.

```C++
      using AbsRealPtr = std::unique_ptr<RooAbsReal>;
      AbsRealPtr{model.derivative(mu, normSet, 1)}->Print();
      AbsRealPtr{model.derivative(sigma, normSet, 1)}->Print();
      AbsRealPtr{model.derivative(c, normSet, 1)}->Print();
      AbsRealPtr{model.derivative(f, normSet, 1)}->Print();
```

*Note: the `derivative()` function returns a new RooFit function that is owned by the caller, so don't forget to use smart pointers here (`std::unique_ptr`).*

The output will be:

```
RooDerivative::model_DERIV_mu[ nset=(x) function=model x=mu ] = -0.000107064
RooDerivative::model_DERIV_sigma[ nset=(x) function=model x=sigma ] = 0.000408505
RooDerivative::model_DERIV_c[ nset=(x) function=model x=c ] = -0.242612
RooDerivative::model_DERIV_c[ nset=(x) function=model x=c ] = -0.304629
```
Now that we know how to get the gradient of a PDF numerically, let's talk about automatically-generated analytic gradients.

## Implementing this model outside RooFit with the Chain Rule

Let's take a step outside of RooFit for a moment and see how we could implement this PDF with its gradient using the chain rule.

Let's write the code for the Gaussian, the exponential, and the addition in some free functions analogous to the RooFit classes:

```C++
namespace calc {

double gaus(double x, double mu, double sigma) {
   const double c = 1. / std::sqrt(2 * M_PI);
   const double arg = (x - mu) / sigma;
   return c / sigma * std::exp(-0.5 * arg * arg);
}

double expo(double x, double c) {
   return -c * std::exp(c * x);
}

double add(double a, double b, double f) {
   return f * a + (1 - f) * b;
}

} // namespace calc
```
The gradients, we will get automatically with [clad](https://github.com/vgvassilev/clad):
```C++
#include <clad/Differentiator/Differentiator.h>

namespace calc {

auto gaus_clad_grad = clad::gradient(gaus, "x,mu,sigma");
auto expo_clad_grad = clad::gradient(expo, "x,c");
auto add_clad_grad = clad::gradient(add, "a,b,f");

} // namespace calc
```
Now we need to evaluate the free functions in the right order, tracking intermediate results and gradients for the chain rule.

To show how this could be done, I wrote a small ChainRuleDriverSimple class.
You tell it the number of computation graph nodes in the constructor, and then when evaluating the node you pass:
1. The index of the node that is evaluated
2. The free functions for the value and the gradient
3. An list with the input node indices (length must match the numer of computation function parameters)

With this approach, the Gaussian plus exponential example would look like:
```C++
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
```

The ChainRuleDriverSimple code is too verbose to reproduce here fully, but the key points are that it has members to track the values and gradients:
```C++
   std::vector<double> values;
   std::vector<std::vector<double>> gradients;
```
...and a function to apply the chain rule for the gradients (`iVar` is the currently-evaluated node):
```C++
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
```

For simplicity, the driver is computing the gradient with respect to *all variables* and not only the model parameters. But if you compare the derivatives for the model parameters only, you will see that the results match.

## Implementation in RooFit

We have now seen how one can evaluate the gradient of a computation graph that consists of free functions, using the chain rule in driver class that tracks the intermediate results and gradients.

How can we do the same for RooFit computation graphs?

In RooFit, we also have a driver class that can evaluate the computation graph, the [RooFitDriver](https://github.com/root-project/root/blob/master/roofit/roofitcore/src/RooFitDriver.cxx). It is a bit more compicated than the simple driver in the standalone example, because it also allows for the evaluation of nodes in the GPU and it figures out the node evaluation order (which we hardcoded in the example above). The driver class is still actively developed and will get easier to understand in the coming weeks.

The driver has a [RooFitDriver::getVal()](https://github.com/root-project/root/blob/master/roofit/roofitcore/src/RooFitDriver.cxx#L422) function that loops over the computation graph nodes (that inherit from RooAbsArg in RooFit), and evaluates them using `RooAbsArg::computeBatch()`. The intermediate results are stored in the class member `RooFitDriver::_dataMapCPU`.

To get gradient evaluation with the chain rule work in RooFit, we will have to:
* extend the **RooAbsArg** interface by a function that returns the gradient (which might be automatically generated with clad)
* extend the **RooFitDriver** to also track the intermediate gradients

The implementation needs to be more involved than the simple chain rule driver above, because calculating the gradient in all variables would be wasteful. It is also the responsability of the RooFitDriver to trigger the computation of the gradient only in the model parameters.

## First steps

Concretely, the steps to do this would be:
1. Try to compile the standalone example and understand the logic of evaluating a computational graph with the chain rule
2. Set up a RooFit developer environment
3. Add gradient interface to `RooAbsArg` and implement the gradient for the few classes we used in the first example (`RooGaussian`, `RooExponential`, and `RooAddPdf`) - return the gradient in all parameters for now
4. Also evaluate the gradient in the RooFitDriver using the chain rule, see if the result matches with the numeric derivatives
