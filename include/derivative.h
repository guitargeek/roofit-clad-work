#include <clad/Differentiator/Differentiator.h>

#include <iostream>

inline void Exp_pullback(Double_t x, Double_t _d_y, clad::array_ref<Double_t> _d_x) {
    Double_t _t0;
    _t0 = x;
    double Exp_return = exp(_t0);
    goto _label0;
  _label0:
    {
        Double_t _r0 = _d_y * clad::custom_derivatives::exp_pushforward(_t0, 1.).pushforward;
        * _d_x += _r0;
    }
}
static void poisson_pullback(double x, double par, double _d_y, clad::array_ref<double> _d_x, clad::array_ref<double> _d_par) {
    bool _cond0;
    bool _cond1;
    bool _cond2;
    double _t0;
    double _t1;
    double _t2;
    double _t3;
    double _t4;
    double _t5;
    _cond0 = par < 0;
    if (_cond0) {
        Double_t poisson_return = TMath::QuietNaN();
        goto _label0;
    }
    _cond1 = x < 0;
    if (_cond1) {
        int poisson_return = 0;
        goto _label1;
    } else {
        _cond2 = x == 0.;
        if (_cond2) {
            _t0 = -par;
            Double_t poisson_return0 = TMath::Exp(_t0);
            goto _label2;
        } else {
            _t2 = x;
            _t3 = par;
            _t1 = log(_t3);
            _t4 = x + 1.;
            _t5 = _t2 * _t1 - TMath::LnGamma(_t4) - par;
            Double_t poisson_return1 = TMath::Exp(_t5);
            goto _label3;
        }
    }
    if (_cond1)
      _label1:
        ;
    else if (_cond2)
      _label2:
        {
            Double_t _grad0 = 0.;
            Exp_pullback(_t0, _d_y, &_grad0);
            double _r0 = _grad0;
            * _d_par += -_r0;
        }
    else {
      _label3:
        {
            Double_t _grad2 = 0.;
            Exp_pullback(_t5, _d_y, &_grad2);
            double _r1 = _grad2;
            double _r2 = _r1 * _t1;
            * _d_x += _r2;
            double _r3 = _t2 * _r1;
            double _r4 = _r3 * clad::custom_derivatives::log_pushforward(_t3, 1.).pushforward;
            * _d_par += _r4;
            double _r5 = -_r1 * numerical_diff::forward_central_difference(TMath::LnGamma, _t4, 0, 0, _t4);
            * _d_x += _r5;
            * _d_par += -_r1;
        }
    }
    if (_cond0)
      _label0:
        ;
}
static void gauss_pullback(double x, double mean, double sigma, double _d_y, clad::array_ref<double> _d_x, clad::array_ref<double> _d_mean, clad::array_ref<double> _d_sigma) {
    double _t0;
    double _t1;
    double _t2;
    double _t3;
    double _t4;
    double _t5;
    double _t6;
    double _t7;
    double _t8;
    double _t9;
    double _t10;
    Double_t _t11;
    double _t12;
    double _t13;
    double _d_arg = 0;
    const double arg = x - mean;
    double _d_sig = 0;
    const double sig = sigma;
    _t2 = arg;
    _t3 = -0.5 * _t2;
    _t1 = arg;
    _t4 = _t3 * _t1;
    _t6 = sig;
    _t5 = sig;
    _t0 = (_t6 * _t5);
    _t7 = _t4 / _t0;
    double _d_out = 0;
    double out = std::exp(_t7);
    _t11 = TMath::TwoPi();
    _t12 = std::sqrt(_t11);
    _t10 = sigma;
    _t9 = (_t12 * _t10);
    _t13 = 1. / _t9;
    _t8 = out;
    double gauss_return = _t13 * _t8;
    goto _label0;
  _label0:
    {
        double _r9 = _d_y * _t8;
        double _r10 = _r9 / _t9;
        double _r11 = _r9 * -1. / (_t9 * _t9);
        double _r12 = _r11 * _t10;
        Double_t _r13 = _r12 * clad::custom_derivatives::sqrt_pushforward(_t11, 1.).pushforward;
        double _r14 = _t12 * _r11;
        * _d_sigma += _r14;
        double _r15 = _t13 * _d_y;
        _d_out += _r15;
    }
    {
        double _r0 = _d_out * clad::custom_derivatives::exp_pushforward(_t7, 1.).pushforward;
        double _r1 = _r0 / _t0;
        double _r2 = _r1 * _t1;
        double _r3 = _r2 * _t2;
        double _r4 = -0.5 * _r2;
        _d_arg += _r4;
        double _r5 = _t3 * _r1;
        _d_arg += _r5;
        double _r6 = _r0 * -_t4 / (_t0 * _t0);
        double _r7 = _r6 * _t5;
        _d_sig += _r7;
        double _r8 = _t6 * _r6;
        _d_sig += _r8;
    }
    * _d_sigma += _d_sig;
    {
        * _d_x += _d_arg;
        * _d_mean += -_d_arg;
    }
}

static void getBin_pullback(double *binBoundaries, double x, int _d_y, clad::array_ref<double> _d_binBoundaries, clad::array_ref<double> _d_x) {
    unsigned long _t0;
    int _d_ibin = 0;
    int ibin = 0;
    _t0 = 0;
    while (binBoundaries[ibin + 1] < x)
        {
            _t0++;
            ibin++;
        }
    int getBin_return = ibin;
    goto _label0;
  _label0:
    _d_ibin += _d_y;
    while (_t0)
        {
            {
            }
            _t0--;
        }
}

void nll_grad(double *in, clad::array_ref<double> _d_in) {
    double _t0;
    double _t1;
    double _t2;
    double _t3;
    double _t4;
    double _t5;
    double _t6;
    double _t7;
    double _t8;
    double _t9;
    unsigned long _t10;
    clad::tape<int> _t11 = {};
    clad::tape<double> _t13 = {};
    unsigned long _t14;
    clad::tape<clad::array<double> > _t15 = {};
    clad::tape<int> _t16 = {};
    clad::tape<double> _t18 = {};
    clad::tape<clad::array<double> > _t19 = {};
    clad::tape<int> _t20 = {};
    clad::tape<double> _t22 = {};
    clad::tape<clad::array<double> > _t23 = {};
    clad::tape<int> _t24 = {};
    clad::tape<double> _t26 = {};
    clad::tape<double> _t27 = {};
    clad::tape<unsigned int> _t28 = {};
    clad::tape<double> _t30 = {};
    clad::tape<double> _t31 = {};
    clad::tape<double> _t33 = {};
    clad::tape<double> _t35 = {};
    clad::tape<double> _t36 = {};
    clad::tape<unsigned int> _t37 = {};
    clad::tape<double> _t39 = {};
    clad::tape<int> _t40 = {};
    clad::tape<double> _t42 = {};
    clad::tape<double> _t44 = {};
    clad::tape<double> _t45 = {};
    clad::tape<unsigned int> _t46 = {};
    clad::tape<double> _t48 = {};
    clad::tape<int> _t49 = {};
    clad::tape<double> _t51 = {};
    clad::tape<double> _t53 = {};
    clad::tape<double> _t54 = {};
    clad::tape<int> _t55 = {};
    clad::tape<double> _t57 = {};
    double _d_nomGammaB1 = 0;
    double nomGammaB1 = 400;
    double _d_nomGammaB2 = 0;
    double nomGammaB2 = 100;
    double _d_nominalLumi = 0;
    double nominalLumi = 1;
    _t0 = nomGammaB1;
    _t2 = nomGammaB1;
    _t1 = in[0];
    _t3 = (_t2 * _t1);
    _t4 = nomGammaB2;
    _t6 = nomGammaB2;
    _t5 = in[1];
    _t7 = (_t6 * _t5);
    _t8 = in[2];
    _t9 = nominalLumi;
    clad::array<double> _d_constraint(3UL);
    double constraint[3] = {ExRooPoisson::poisson(_t0, _t3), ExRooPoisson::poisson(_t4, _t7), ExRooGaussian::gauss(_t8, _t9, 0.10000000000000001)};
    double _d_cnstSum = 0;
    double cnstSum = 0;
    clad::array<double> _d_x(2UL);
    double x[2] = {1.25, 1.75};
    clad::array<double> _d_sig(2UL);
    double sig[2] = {20, 10};
    clad::array<double> _d_binBoundaries1(3UL);
    double binBoundaries1[3] = {1, 1.5, 2};
    clad::array<double> _d_bgk1(2UL);
    double bgk1[2] = {100, 0};
    clad::array<double> _d_binBoundaries2(3UL);
    double binBoundaries2[3] = {1, 1.5, 2};
    clad::array<double> _d_histVals(2UL);
    double histVals[2] = {in[0], in[1]};
    clad::array<double> _d_bgk2(2UL);
    double bgk2[2] = {0, 100};
    clad::array<double> _d_binBoundaries3(3UL);
    double binBoundaries3[3] = {1, 1.5, 2};
    clad::array<double> _d_weights(2UL);
    double weights[2] = {122., 112.};
    _t10 = 0;
    {
        int _d_i = 0;
        for (int i = 0; i < 3; i++) {
            _t10++;
            cnstSum -= std::log(clad::push(_t13, constraint[clad::push(_t11, i)]));
        }
    }
    double _d_mu = 0;
    double mu = 0;
    double _d_temp = 0;
    double temp;
    double _d_nllSum = 0;
    double nllSum = 0;
    unsigned int _d_b1 = 0, _d_b2 = 0, _d_b3 = 0;
    unsigned int b1, b2, b3;
    _t14 = 0;
    {
        int _d_iB = 0;
        for (int iB = 0; iB < 2; iB++) {
            _t14++;
            clad::push(_t15, binBoundaries1 , 3UL);
            b1 = ExRooHistFunc::getBin(binBoundaries1, clad::push(_t18, x[clad::push(_t16, iB)]));
            clad::push(_t19, binBoundaries2 , 3UL);
            b2 = ExRooHistFunc::getBin(binBoundaries2, clad::push(_t22, x[clad::push(_t20, iB)]));
            clad::push(_t23, binBoundaries3 , 3UL);
            b3 = ExRooHistFunc::getBin(binBoundaries3, clad::push(_t26, x[clad::push(_t24, iB)]));
            mu = 0;
            mu += clad::push(_t30, sig[clad::push(_t28, b1)]) * clad::push(_t27, (clad::push(_t33, in[3]) * clad::push(_t31, in[2])));
            mu += clad::push(_t42, (clad::push(_t39, bgk1[clad::push(_t37, b2)]) * clad::push(_t36, histVals[clad::push(_t40, iB)]))) * clad::push(_t35, (in[2] * 1.));
            mu += clad::push(_t51, (clad::push(_t48, bgk2[clad::push(_t46, b3)]) * clad::push(_t45, histVals[clad::push(_t49, iB)]))) * clad::push(_t44, (in[2] * 1.));
            temp = std::log(clad::push(_t53, mu));
            nllSum -= -mu + clad::push(_t57, weights[clad::push(_t55, iB)]) * clad::push(_t54, temp);
        }
    }
    double nll_return = cnstSum + nllSum;
    goto _label0;
  _label0:
    {
        _d_cnstSum += 1;
        _d_nllSum += 1;
    }
    for (; _t14; _t14--) {
        {
            double _r_d9 = _d_nllSum;
            _d_nllSum += _r_d9;
            _d_mu += - -_r_d9;
            double _r36 = -_r_d9 * clad::pop(_t54);
            int _t56 = clad::pop(_t55);
            _d_weights[_t56] += _r36;
            double _r37 = clad::pop(_t57) * -_r_d9;
            _d_temp += _r37;
            _d_nllSum -= _r_d9;
        }
        {
            double _r_d8 = _d_temp;
            double _r35 = _r_d8 * clad::custom_derivatives::log_pushforward(clad::pop(_t53), 1.).pushforward;
            _d_mu += _r35;
            _d_temp -= _r_d8;
        }
        {
            double _r_d7 = _d_mu;
            _d_mu += _r_d7;
            double _r30 = _r_d7 * clad::pop(_t44);
            double _r31 = _r30 * clad::pop(_t45);
            unsigned int _t47 = clad::pop(_t46);
            _d_bgk2[_t47] += _r31;
            double _r32 = clad::pop(_t48) * _r30;
            int _t50 = clad::pop(_t49);
            _d_histVals[_t50] += _r32;
            double _r33 = clad::pop(_t51) * _r_d7;
            double _r34 = _r33 * 1.;
            int _t52 = 2;
            _d_in[_t52] += _r34;
            _d_mu -= _r_d7;
        }
        {
            double _r_d6 = _d_mu;
            _d_mu += _r_d6;
            double _r25 = _r_d6 * clad::pop(_t35);
            double _r26 = _r25 * clad::pop(_t36);
            unsigned int _t38 = clad::pop(_t37);
            _d_bgk1[_t38] += _r26;
            double _r27 = clad::pop(_t39) * _r25;
            int _t41 = clad::pop(_t40);
            _d_histVals[_t41] += _r27;
            double _r28 = clad::pop(_t42) * _r_d6;
            double _r29 = _r28 * 1.;
            int _t43 = 2;
            _d_in[_t43] += _r29;
            _d_mu -= _r_d6;
        }
        {
            double _r_d5 = _d_mu;
            _d_mu += _r_d5;
            double _r21 = _r_d5 * clad::pop(_t27);
            unsigned int _t29 = clad::pop(_t28);
            _d_sig[_t29] += _r21;
            double _r22 = clad::pop(_t30) * _r_d5;
            double _r23 = _r22 * clad::pop(_t31);
            int _t32 = 3;
            _d_in[_t32] += _r23;
            double _r24 = clad::pop(_t33) * _r22;
            int _t34 = 2;
            _d_in[_t34] += _r24;
            _d_mu -= _r_d5;
        }
        {
            double _r_d4 = _d_mu;
            _d_mu -= _r_d4;
        }
        {
            unsigned int _r_d3 = _d_b3;
            clad::array<double> _r19 = clad::pop(_t23);
            double _grad12 = 0.;
            getBin_pullback(_r19, clad::pop(_t26), _r_d3, _d_binBoundaries3, &_grad12);
            clad::array<double> _r18(_d_binBoundaries3);
            double _r20 = _grad12;
            int _t25 = clad::pop(_t24);
            _d_x[_t25] += _r20;
            _d_b3 -= _r_d3;
        }
        {
            unsigned int _r_d2 = _d_b2;
            clad::array<double> _r16 = clad::pop(_t19);
            double _grad10 = 0.;
            getBin_pullback(_r16, clad::pop(_t22), _r_d2, _d_binBoundaries2, &_grad10);
            clad::array<double> _r15(_d_binBoundaries2);
            double _r17 = _grad10;
            int _t21 = clad::pop(_t20);
            _d_x[_t21] += _r17;
            _d_b2 -= _r_d2;
        }
        {
            unsigned int _r_d1 = _d_b1;
            clad::array<double> _r13 = clad::pop(_t15);
            double _grad8 = 0.;
            getBin_pullback(_r13, clad::pop(_t18), _r_d1, _d_binBoundaries1, &_grad8);
            clad::array<double> _r12(_d_binBoundaries1);
            double _r14 = _grad8;
            int _t17 = clad::pop(_t16);
            _d_x[_t17] += _r14;
            _d_b1 -= _r_d1;
        }
    }
    for (; _t10; _t10--) {
        {
            double _r_d0 = _d_cnstSum;
            _d_cnstSum += _r_d0;
            double _r11 = -_r_d0 * clad::custom_derivatives::log_pushforward(clad::pop(_t13), 1.).pushforward;
            int _t12 = clad::pop(_t11);
            _d_constraint[_t12] += _r11;
            _d_cnstSum -= _r_d0;
        }
    }
    {
        _d_in[0] += _d_histVals[0];
        _d_in[1] += _d_histVals[1];
    }
    {
        double _grad0 = 0.;
        double _grad1 = 0.;
        poisson_pullback(_t0, _t3, _d_constraint[0], &_grad0, &_grad1);
        double _r0 = _grad0;
        _d_nomGammaB1 += _r0;
        double _r1 = _grad1;
        double _r2 = _r1 * _t1;
        _d_nomGammaB1 += _r2;
        double _r3 = _t2 * _r1;
        _d_in[0] += _r3;
        double _grad2 = 0.;
        double _grad3 = 0.;
        poisson_pullback(_t4, _t7, _d_constraint[1], &_grad2, &_grad3);
        double _r4 = _grad2;
        _d_nomGammaB2 += _r4;
        double _r5 = _grad3;
        double _r6 = _r5 * _t5;
        _d_nomGammaB2 += _r6;
        double _r7 = _t6 * _r5;
        _d_in[1] += _r7;
        double _grad4 = 0.;
        double _grad5 = 0.;
        double _grad6 = 0.;
        gauss_pullback(_t8, _t9, 0.10000000000000001, _d_constraint[2], &_grad4, &_grad5, &_grad6);
        double _r8 = _grad4;
        _d_in[2] += _r8;
        double _r9 = _grad5;
        _d_nominalLumi += _r9;
        double _r10 = _grad6;
    }
}
