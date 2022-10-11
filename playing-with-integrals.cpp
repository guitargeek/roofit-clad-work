// clang++ test.cpp -S -emit-llvm -o input.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops --std=c++17
// opt input.ll --enable-new-pm=0 -load=Enzyme/build/Enzyme/LLDEnzyme-14.so  -enzyme -o output.ll -S
// opt output.ll -O2 -o output_opt.ll -S
// clang++ output_opt.ll -o a.exe
// clang++ output.ll -O3 -o a.exe

#include <cmath>
#include <iostream>

class ConstVar {
public:
    constexpr ConstVar(double val) : _val{val} {}

    constexpr double val() const { return _val; }

    constexpr static bool isObs() { return false; }
    constexpr static bool isVec() { return false; }

private:
    const double _val;
};


template<bool hasMin_t=false, bool hasMax_t=false>
class RealObs {
public:

    constexpr RealObs() {
    }

    void setVal(double val) { _val = val; }

    // should be disabled at compile time with enable_if in case this variable has no ranges
    void setMin(double min) { _min = min; }
    void setMax(double max) { _max = max; }

    constexpr double val() const { return _val; }
    constexpr double min() const { return _min; }
    constexpr double max() const { return _max; }

    constexpr static bool isObs() { return true; }
    constexpr static bool hasMin() { return hasMin_t; }
    constexpr static bool hasMax() { return hasMax_t; }
    constexpr static bool isVec() { return false; }

private:
    double _val = 0.0;
    double _min = 0.0;
    double _max = 0.0;
};


template<class X_t, class Mu_t, class Sigma_t>
class Gaussian {

public:

    constexpr Gaussian(X_t const& x, Mu_t const& mu, Sigma_t const& sigma)
        : _x{x}, _mu{mu}, _sigma{sigma}
    {
        static_assert(!Mu_t::isObs(), "\"mu\" can't be an RealObs!");
        static_assert(!Sigma_t::isObs(), "\"sigma\" can't be an RealObs!");
    }

    inline void normalize() {
        _invNorm = 1./gaussNorm();
    }


    inline double val() const
    {
        double val = gaussUnnormalized(_x.val(), _mu.val(), _sigma.val());

        if constexpr(X_t::isObs()) {
            val *= _invNorm;
        }

        return val;
    }

private:
    inline double gaussNorm() const
    {
        if constexpr(!X_t::isObs()) return 1.0;

        double norm = 0.0;
        if constexpr(X_t::hasMax()) {
            norm += gaussIntegral(_x.max(), _mu.val(), _sigma.val());
        } else {
            norm += gaussIntegralFull(_mu.val(), _sigma.val());
        }
        if constexpr(X_t::hasMin()) {
            norm -= gaussIntegral(_x.min(), _mu.val(), _sigma.val());
        }

        return norm;
    }

    static inline double gaussUnnormalized(double x, double mu, double sigma)
    {
        const double arg = x - mu;
        const double inv = 1./(sigma * sigma);
        return std::exp(-arg * arg * inv);
    }

    static inline double gaussIntegral(double x, double mu, double sigma)
    {
        const double inv = 1. / sigma / std::sqrt(2);

        return std::sqrt(2 * M_PI) * sigma * 0.5 * (1 + std::erf((x - mu) * inv));
    }

    static inline double gaussIntegralFull(double mu, double sigma)
    {
        return std::sqrt(2 * M_PI) * sigma;
    }

    X_t const& _x;
    Mu_t const& _mu;
    Sigma_t const& _sigma;
    double _invNorm = 1.0;
};


//class NLL {
//};

//extern double __enzyme_autodiff(void*, double);

double foobar(double muVal)
{
    RealObs<true,true> x;
    x.setMin(-10);
    x.setMax(10);

    ConstVar mu{muVal};
    ConstVar sigma{1.0};
    Gaussian pdf{x, mu, sigma};

    pdf.normalize();

    return pdf.val();
}

double dfoobar(double x) {
    // This returns the derivative of foobar or 2 * x
    //return __enzyme_autodiff((void*) foobar, x);
}

int main()
{

    std::cout << foobar(0.0) << "   " << dfoobar(0.0) << std::endl;
    std::cout << foobar(1.0) << "   " << dfoobar(1.0) << std::endl;
    std::cout << foobar(-1.0) << "   " << dfoobar(-1.0) << std::endl;

    //std::vector<double> xVals{-1.0, 0.0, 1.0};


    //double nll = 0.0;

    //pdf.normalize();

    //for(const double xVal : xVals) {

        //x.setVal(xVal);

        //std::cout << pdf.val() << std::endl;
    //}

}

void test() {
    main();
}
