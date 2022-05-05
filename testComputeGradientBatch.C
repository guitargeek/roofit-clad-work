#include <RooGaussian.h>
#include <RooRealVar.h>

#include <RooBatchCompute.h>

void testComputeGradientBatch() {

    RooRealVar x{"x", "x", 3.5, -10, 10};
    RooRealVar mu{"mu", "mu", 3.0, -10, 10};
    RooRealVar sigma{"sigma", "sigma", 2.0, 0.1, 10};

    RooGaussian gauss{"gauss", "gauss", x, mu, sigma};

    // Evaluate once with getValues to trigger the initialization of RooBatchCompute
    RooDataSet data{}; // empy dataset
    gauss.getValues(data);

    RooArgSet normSet{x};

    std::cout << gauss.getVal(normSet) << std::endl;

    std::map<RooAbsReal const*, double> owningDataMap{
        {&x, x.getVal()},
        {&mu, mu.getVal()},
        {&sigma, sigma.getVal()},
        {&gauss, 0.0}
    };

    RooBatchCompute::DataMap dataMap{
        {&x, {&owningDataMap[&x], 1}},
        {&mu, {&owningDataMap[&mu], 1}},
        {&sigma, {&owningDataMap[&sigma], 1}},
        {&gauss, {&owningDataMap[&gauss], 1}}
    };
    
    // hack because computeBatch is not public in RooGaussian (needs fix in ROOT!)
    auto& pdf = static_cast<RooAbsReal&>(gauss);

    pdf.computeBatch(nullptr, &owningDataMap.at(&gauss), 1, dataMap);

    for(auto const& item : dataMap) {
        std::cout << item.first->GetName() << "    " << item.second[0] << std::endl;
    }

    std::vector<double> gradientOutput;
    gradientOutput.resize(3);

    // This function needs to be implemented for RooGaussian
    // pdf.computeGradientBatch(nullptr, gradientOutput.data(), 1, dataMap);

    std::cout << std::endl;
    for(auto& x : gradientOutput) {
        std::cout << x << std::endl;
    }
}
