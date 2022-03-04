# RooFit clad work

Repository for sharing code and information in the RooFit AD with clad effort.

## Compiling test codes

The test codes in this repo must be comiled with the [root-compile](./root-compile) script in this repository.
You can't run these macros in the ROOT interpreter.

## Interfacing the analytical gradient with Minuit2 via the RooMinimizer

A first prototype for the plumming between RooFit and Minuit2 is implemented in [this `RooCladMinimizerFcn` branch](https://github.com/guitargeek/root/tree/RooCladMinimizerFcn_1).
This implementation reuses much of the `RooGradMinimizerFcn` logic, implementing `IMultiGradFunction::DoDerivative` where the derivatives are filled from cached gradient results.
You can test the `RooCladMinimizerFcn` with the [example_1.cpp](./example_1.cpp) code, which implements a complete Gaussian likelihood in one RooFit object that can also return the gradient.

Several known aspects about the RooMinimizer interface need to be improved in the future:

* implement `IMultiGradFunction::Gradient` directly
* avoid code duplication with `RooGradMinimizerFcn`

## Passing the second derivative to Minuit

Work on this hasn't started yet, but it's important to do so for example to compute the initial step size. Right now, the initial gradient is calculated numerically, even if an analytical gradient calculator is provided (see [PR #7001](https://github.com/root-project/root/pull/7001)).
