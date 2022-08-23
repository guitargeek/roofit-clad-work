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

void nll_grad(double gammaB1, double gammaB2, double gamma1_B1, double gamma1_B2, double lumi, double SigXsecOverSM, clad::array_ref<double> _d_gammaB1, clad::array_ref<double> _d_gammaB2, clad::array_ref<double> _d_gamma1_B1, clad::array_ref<double> _d_gamma1_B2, clad::array_ref<double> _d_lumi, clad::array_ref<double> _d_SigXsecOverSM) {
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
    double _t11;
    double _t12;
    double _t13;
    double _t14;
    double _t15;
    double _t16;
    double _t17;
    unsigned long _t18;
    clad::tape<int> _t19 = {};
    clad::tape<double> _t21 = {};
    unsigned long _t22;
    clad::tape<clad::array<double> > _t23 = {};
    clad::tape<int> _t24 = {};
    clad::tape<double> _t26 = {};
    clad::tape<clad::array<double> > _t27 = {};
    clad::tape<int> _t28 = {};
    clad::tape<double> _t30 = {};
    clad::tape<clad::array<double> > _t31 = {};
    clad::tape<int> _t32 = {};
    clad::tape<double> _t34 = {};
    clad::tape<double> _t35 = {};
    clad::tape<unsigned int> _t36 = {};
    clad::tape<double> _t38 = {};
    clad::tape<double> _t39 = {};
    clad::tape<double> _t40 = {};
    clad::tape<double> _t41 = {};
    clad::tape<double> _t42 = {};
    clad::tape<unsigned int> _t43 = {};
    clad::tape<double> _t45 = {};
    clad::tape<int> _t46 = {};
    clad::tape<double> _t48 = {};
    clad::tape<double> _t49 = {};
    clad::tape<double> _t50 = {};
    clad::tape<unsigned int> _t51 = {};
    clad::tape<double> _t53 = {};
    clad::tape<int> _t54 = {};
    clad::tape<double> _t56 = {};
    clad::tape<double> _t57 = {};
    clad::tape<double> _t58 = {};
    clad::tape<int> _t59 = {};
    clad::tape<double> _t61 = {};
    unsigned long _t62;
    clad::tape<clad::array<double> > _t63 = {};
    clad::tape<int> _t64 = {};
    clad::tape<double> _t66 = {};
    clad::tape<clad::array<double> > _t67 = {};
    clad::tape<int> _t68 = {};
    clad::tape<double> _t70 = {};
    clad::tape<clad::array<double> > _t71 = {};
    clad::tape<int> _t72 = {};
    clad::tape<double> _t74 = {};
    clad::tape<double> _t75 = {};
    clad::tape<unsigned int> _t76 = {};
    clad::tape<double> _t78 = {};
    clad::tape<double> _t79 = {};
    clad::tape<double> _t80 = {};
    clad::tape<double> _t81 = {};
    clad::tape<double> _t82 = {};
    clad::tape<unsigned int> _t83 = {};
    clad::tape<double> _t85 = {};
    clad::tape<int> _t86 = {};
    clad::tape<double> _t88 = {};
    clad::tape<double> _t89 = {};
    clad::tape<double> _t90 = {};
    clad::tape<unsigned int> _t91 = {};
    clad::tape<double> _t93 = {};
    clad::tape<int> _t94 = {};
    clad::tape<double> _t96 = {};
    clad::tape<double> _t97 = {};
    clad::tape<double> _t98 = {};
    clad::tape<int> _t99 = {};
    clad::tape<double> _t101 = {};
    double _d_nomGammaB1 = 0;
    double nomGammaB1 = 400;
    double _d_nomGammaB2 = 0;
    double nomGammaB2 = 100;
    double _d_nomGamma1_B1 = 0;
    double nomGamma1_B1 = 400;
    double _d_nomGamma1_B2 = 0;
    double nomGamma1_B2 = 100;
    double _d_nominalLumi = 0;
    double nominalLumi = 1;
    _t0 = nomGammaB1;
    _t2 = nomGammaB1;
    _t1 = gammaB1;
    _t3 = (_t2 * _t1);
    _t4 = nomGammaB2;
    _t6 = nomGammaB2;
    _t5 = gammaB2;
    _t7 = (_t6 * _t5);
    _t8 = nomGamma1_B1;
    _t10 = nomGamma1_B1;
    _t9 = gamma1_B1;
    _t11 = (_t10 * _t9);
    _t12 = nomGamma1_B2;
    _t14 = nomGamma1_B2;
    _t13 = gamma1_B2;
    _t15 = (_t14 * _t13);
    _t16 = lumi;
    _t17 = nominalLumi;
    clad::array<double> _d_constraint(5UL);
    double constraint[5] = {ExRooPoisson::poisson(_t0, _t3), ExRooPoisson::poisson(_t4, _t7), ExRooPoisson::poisson(_t8, _t11), ExRooPoisson::poisson(_t12, _t15), ExRooGaussian::gauss(_t16, _t17, 0.10000000000000001)};
    double _d_cnstSum = 0;
    double cnstSum = 0;
    clad::array<double> _d_x(2UL);
    double x[2] = {1.25, 1.75};
    clad::array<double> _d_sig0(2UL);
    double sig0[2] = {20, 10};
    clad::array<double> _d_binBoundaries1(3UL);
    double binBoundaries1[3] = {1, 1.5, 2};
    clad::array<double> _d_bgk0_1(2UL);
    double bgk0_1[2] = {100, 0};
    clad::array<double> _d_binBoundaries2(3UL);
    double binBoundaries2[3] = {1, 1.5, 2};
    clad::array<double> _d_histVals0(2UL);
    double histVals0[2] = {gammaB1, gammaB2};
    clad::array<double> _d_bgk0_2(2UL);
    double bgk0_2[2] = {0, 100};
    clad::array<double> _d_binBoundaries3(3UL);
    double binBoundaries3[3] = {1, 1.5, 2};
    clad::array<double> _d_nllSum0_weights(2UL);
    double nllSum0_weights[2] = {122., 112.};
    clad::array<double> _d_sig1(2UL);
    double sig1[2] = {20, 10};
    clad::array<double> _d_binBoundaries4(3UL);
    double binBoundaries4[3] = {1, 1.5, 2};
    clad::array<double> _d_bgk1_1(2UL);
    double bgk1_1[2] = {100, 0};
    clad::array<double> _d_binBoundaries5(3UL);
    double binBoundaries5[3] = {1, 1.5, 2};
    clad::array<double> _d_histVals1(2UL);
    double histVals1[2] = {gamma1_B1, gamma1_B2};
    clad::array<double> _d_bgk1_2(2UL);
    double bgk1_2[2] = {0, 100};
    clad::array<double> _d_binBoundaries6(3UL);
    double binBoundaries6[3] = {1, 1.5, 2};
    clad::array<double> _d_nllSum1_weights(2UL);
    double nllSum1_weights[2] = {122., 112.};
    _t18 = 0;
    {
        int _d_i = 0;
        for (int i = 0; i < 5; i++) {
            _t18++;
            cnstSum -= std::log(clad::push(_t21, constraint[clad::push(_t19, i)]));
        }
    }
    unsigned int _d_b1 = 0, _d_b2 = 0, _d_b3 = 0, _d_b4 = 0, _d_b5 = 0, _d_b6 = 0;
    unsigned int b1, b2, b3, b4, b5, b6;
    double _d_nllSum0 = 0;
    double nllSum0 = 0;
    double _d_mu0 = 0;
    double mu0 = 0;
    double _d_temp = 0;
    double temp;
    _t22 = 0;
    {
        int _d_iB = 0;
        for (int iB = 0; iB < 2; iB++) {
            _t22++;
            clad::push(_t23, binBoundaries1 , 3UL);
            b1 = ExRooHistFunc::getBin(binBoundaries1, clad::push(_t26, x[clad::push(_t24, iB)]));
            clad::push(_t27, binBoundaries2 , 3UL);
            b2 = ExRooHistFunc::getBin(binBoundaries2, clad::push(_t30, x[clad::push(_t28, iB)]));
            clad::push(_t31, binBoundaries3 , 3UL);
            b3 = ExRooHistFunc::getBin(binBoundaries3, clad::push(_t34, x[clad::push(_t32, iB)]));
            mu0 += clad::push(_t38, sig0[clad::push(_t36, b1)]) * clad::push(_t35, (clad::push(_t40, SigXsecOverSM) * clad::push(_t39, lumi)));
            mu0 += clad::push(_t48, (clad::push(_t45, bgk0_1[clad::push(_t43, b2)]) * clad::push(_t42, histVals0[clad::push(_t46, iB)]))) * clad::push(_t41, (lumi * 1.));
            mu0 += clad::push(_t56, (clad::push(_t53, bgk0_2[clad::push(_t51, b3)]) * clad::push(_t50, histVals0[clad::push(_t54, iB)]))) * clad::push(_t49, (lumi * 1.));
            temp = std::log(clad::push(_t57, mu0));
            nllSum0 -= -mu0 + clad::push(_t61, nllSum0_weights[clad::push(_t59, iB)]) * clad::push(_t58, temp);
        }
    }
    double _d_nllSum1 = 0;
    double nllSum1 = 0;
    double _d_mu1 = 0;
    double mu1 = 0;
    _t62 = 0;
    {
        int _d_iB = 0;
        for (int iB = 0; iB < 2; iB++) {
            _t62++;
            clad::push(_t63, binBoundaries4 , 3UL);
            b4 = ExRooHistFunc::getBin(binBoundaries4, clad::push(_t66, x[clad::push(_t64, iB)]));
            clad::push(_t67, binBoundaries5 , 3UL);
            b5 = ExRooHistFunc::getBin(binBoundaries5, clad::push(_t70, x[clad::push(_t68, iB)]));
            clad::push(_t71, binBoundaries6 , 3UL);
            b6 = ExRooHistFunc::getBin(binBoundaries6, clad::push(_t74, x[clad::push(_t72, iB)]));
            mu1 += clad::push(_t78, sig1[clad::push(_t76, b4)]) * clad::push(_t75, (clad::push(_t80, SigXsecOverSM) * clad::push(_t79, lumi)));
            mu1 += clad::push(_t88, (clad::push(_t85, bgk1_1[clad::push(_t83, b5)]) * clad::push(_t82, histVals1[clad::push(_t86, iB)]))) * clad::push(_t81, (lumi * 1.));
            mu1 += clad::push(_t96, (clad::push(_t93, bgk1_2[clad::push(_t91, b6)]) * clad::push(_t90, histVals1[clad::push(_t94, iB)]))) * clad::push(_t89, (lumi * 1.));
            temp = std::log(clad::push(_t97, mu1));
            nllSum1 -= -mu1 + clad::push(_t101, nllSum1_weights[clad::push(_t99, iB)]) * clad::push(_t98, temp);
        }
    }
    double nll_return = cnstSum + nllSum0 + nllSum1;
    goto _label0;
  _label0:
    {
        _d_cnstSum += 1;
        _d_nllSum0 += 1;
        _d_nllSum1 += 1;
    }
    for (; _t62; _t62--) {
        {
            double _r_d16 = _d_nllSum1;
            _d_nllSum1 += _r_d16;
            _d_mu1 += - -_r_d16;
            double _r70 = -_r_d16 * clad::pop(_t98);
            int _t100 = clad::pop(_t99);
            _d_nllSum1_weights[_t100] += _r70;
            double _r71 = clad::pop(_t101) * -_r_d16;
            _d_temp += _r71;
            _d_nllSum1 -= _r_d16;
        }
        {
            double _r_d15 = _d_temp;
            double _r69 = _r_d15 * clad::custom_derivatives::log_pushforward(clad::pop(_t97), 1.).pushforward;
            _d_mu1 += _r69;
            _d_temp -= _r_d15;
        }
        {
            double _r_d14 = _d_mu1;
            _d_mu1 += _r_d14;
            double _r64 = _r_d14 * clad::pop(_t89);
            double _r65 = _r64 * clad::pop(_t90);
            unsigned int _t92 = clad::pop(_t91);
            _d_bgk1_2[_t92] += _r65;
            double _r66 = clad::pop(_t93) * _r64;
            int _t95 = clad::pop(_t94);
            _d_histVals1[_t95] += _r66;
            double _r67 = clad::pop(_t96) * _r_d14;
            double _r68 = _r67 * 1.;
            * _d_lumi += _r68;
            _d_mu1 -= _r_d14;
        }
        {
            double _r_d13 = _d_mu1;
            _d_mu1 += _r_d13;
            double _r59 = _r_d13 * clad::pop(_t81);
            double _r60 = _r59 * clad::pop(_t82);
            unsigned int _t84 = clad::pop(_t83);
            _d_bgk1_1[_t84] += _r60;
            double _r61 = clad::pop(_t85) * _r59;
            int _t87 = clad::pop(_t86);
            _d_histVals1[_t87] += _r61;
            double _r62 = clad::pop(_t88) * _r_d13;
            double _r63 = _r62 * 1.;
            * _d_lumi += _r63;
            _d_mu1 -= _r_d13;
        }
        {
            double _r_d12 = _d_mu1;
            _d_mu1 += _r_d12;
            double _r55 = _r_d12 * clad::pop(_t75);
            unsigned int _t77 = clad::pop(_t76);
            _d_sig1[_t77] += _r55;
            double _r56 = clad::pop(_t78) * _r_d12;
            double _r57 = _r56 * clad::pop(_t79);
            * _d_SigXsecOverSM += _r57;
            double _r58 = clad::pop(_t80) * _r56;
            * _d_lumi += _r58;
            _d_mu1 -= _r_d12;
        }
        {
            unsigned int _r_d11 = _d_b6;
            clad::array<double> _r53 = clad::pop(_t71);
            double _grad22 = 0.;
            getBin_pullback(_r53, clad::pop(_t74), _r_d11, _d_binBoundaries6, &_grad22);
            clad::array<double> _r52(_d_binBoundaries6);
            double _r54 = _grad22;
            int _t73 = clad::pop(_t72);
            _d_x[_t73] += _r54;
            _d_b6 -= _r_d11;
        }
        {
            unsigned int _r_d10 = _d_b5;
            clad::array<double> _r50 = clad::pop(_t67);
            double _grad20 = 0.;
            getBin_pullback(_r50, clad::pop(_t70), _r_d10, _d_binBoundaries5, &_grad20);
            clad::array<double> _r49(_d_binBoundaries5);
            double _r51 = _grad20;
            int _t69 = clad::pop(_t68);
            _d_x[_t69] += _r51;
            _d_b5 -= _r_d10;
        }
        {
            unsigned int _r_d9 = _d_b4;
            clad::array<double> _r47 = clad::pop(_t63);
            double _grad18 = 0.;
            getBin_pullback(_r47, clad::pop(_t66), _r_d9, _d_binBoundaries4, &_grad18);
            clad::array<double> _r46(_d_binBoundaries4);
            double _r48 = _grad18;
            int _t65 = clad::pop(_t64);
            _d_x[_t65] += _r48;
            _d_b4 -= _r_d9;
        }
    }
    for (; _t22; _t22--) {
        {
            double _r_d8 = _d_nllSum0;
            _d_nllSum0 += _r_d8;
            _d_mu0 += - -_r_d8;
            double _r44 = -_r_d8 * clad::pop(_t58);
            int _t60 = clad::pop(_t59);
            _d_nllSum0_weights[_t60] += _r44;
            double _r45 = clad::pop(_t61) * -_r_d8;
            _d_temp += _r45;
            _d_nllSum0 -= _r_d8;
        }
        {
            double _r_d7 = _d_temp;
            double _r43 = _r_d7 * clad::custom_derivatives::log_pushforward(clad::pop(_t57), 1.).pushforward;
            _d_mu0 += _r43;
            _d_temp -= _r_d7;
        }
        {
            double _r_d6 = _d_mu0;
            _d_mu0 += _r_d6;
            double _r38 = _r_d6 * clad::pop(_t49);
            double _r39 = _r38 * clad::pop(_t50);
            unsigned int _t52 = clad::pop(_t51);
            _d_bgk0_2[_t52] += _r39;
            double _r40 = clad::pop(_t53) * _r38;
            int _t55 = clad::pop(_t54);
            _d_histVals0[_t55] += _r40;
            double _r41 = clad::pop(_t56) * _r_d6;
            double _r42 = _r41 * 1.;
            * _d_lumi += _r42;
            _d_mu0 -= _r_d6;
        }
        {
            double _r_d5 = _d_mu0;
            _d_mu0 += _r_d5;
            double _r33 = _r_d5 * clad::pop(_t41);
            double _r34 = _r33 * clad::pop(_t42);
            unsigned int _t44 = clad::pop(_t43);
            _d_bgk0_1[_t44] += _r34;
            double _r35 = clad::pop(_t45) * _r33;
            int _t47 = clad::pop(_t46);
            _d_histVals0[_t47] += _r35;
            double _r36 = clad::pop(_t48) * _r_d5;
            double _r37 = _r36 * 1.;
            * _d_lumi += _r37;
            _d_mu0 -= _r_d5;
        }
        {
            double _r_d4 = _d_mu0;
            _d_mu0 += _r_d4;
            double _r29 = _r_d4 * clad::pop(_t35);
            unsigned int _t37 = clad::pop(_t36);
            _d_sig0[_t37] += _r29;
            double _r30 = clad::pop(_t38) * _r_d4;
            double _r31 = _r30 * clad::pop(_t39);
            * _d_SigXsecOverSM += _r31;
            double _r32 = clad::pop(_t40) * _r30;
            * _d_lumi += _r32;
            _d_mu0 -= _r_d4;
        }
        {
            unsigned int _r_d3 = _d_b3;
            clad::array<double> _r27 = clad::pop(_t31);
            double _grad16 = 0.;
            getBin_pullback(_r27, clad::pop(_t34), _r_d3, _d_binBoundaries3, &_grad16);
            clad::array<double> _r26(_d_binBoundaries3);
            double _r28 = _grad16;
            int _t33 = clad::pop(_t32);
            _d_x[_t33] += _r28;
            _d_b3 -= _r_d3;
        }
        {
            unsigned int _r_d2 = _d_b2;
            clad::array<double> _r24 = clad::pop(_t27);
            double _grad14 = 0.;
            getBin_pullback(_r24, clad::pop(_t30), _r_d2, _d_binBoundaries2, &_grad14);
            clad::array<double> _r23(_d_binBoundaries2);
            double _r25 = _grad14;
            int _t29 = clad::pop(_t28);
            _d_x[_t29] += _r25;
            _d_b2 -= _r_d2;
        }
        {
            unsigned int _r_d1 = _d_b1;
            clad::array<double> _r21 = clad::pop(_t23);
            double _grad12 = 0.;
            getBin_pullback(_r21, clad::pop(_t26), _r_d1, _d_binBoundaries1, &_grad12);
            clad::array<double> _r20(_d_binBoundaries1);
            double _r22 = _grad12;
            int _t25 = clad::pop(_t24);
            _d_x[_t25] += _r22;
            _d_b1 -= _r_d1;
        }
    }
    for (; _t18; _t18--) {
        {
            double _r_d0 = _d_cnstSum;
            _d_cnstSum += _r_d0;
            double _r19 = -_r_d0 * clad::custom_derivatives::log_pushforward(clad::pop(_t21), 1.).pushforward;
            int _t20 = clad::pop(_t19);
            _d_constraint[_t20] += _r19;
            _d_cnstSum -= _r_d0;
        }
    }
    {
        * _d_gamma1_B1 += _d_histVals1[0];
        * _d_gamma1_B2 += _d_histVals1[1];
    }
    {
        * _d_gammaB1 += _d_histVals0[0];
        * _d_gammaB2 += _d_histVals0[1];
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
        * _d_gammaB1 += _r3;
        double _grad2 = 0.;
        double _grad3 = 0.;
        poisson_pullback(_t4, _t7, _d_constraint[1], &_grad2, &_grad3);
        double _r4 = _grad2;
        _d_nomGammaB2 += _r4;
        double _r5 = _grad3;
        double _r6 = _r5 * _t5;
        _d_nomGammaB2 += _r6;
        double _r7 = _t6 * _r5;
        * _d_gammaB2 += _r7;
        double _grad4 = 0.;
        double _grad5 = 0.;
        poisson_pullback(_t8, _t11, _d_constraint[2], &_grad4, &_grad5);
        double _r8 = _grad4;
        _d_nomGamma1_B1 += _r8;
        double _r9 = _grad5;
        double _r10 = _r9 * _t9;
        _d_nomGamma1_B1 += _r10;
        double _r11 = _t10 * _r9;
        * _d_gamma1_B1 += _r11;
        double _grad6 = 0.;
        double _grad7 = 0.;
        poisson_pullback(_t12, _t15, _d_constraint[3], &_grad6, &_grad7);
        double _r12 = _grad6;
        _d_nomGamma1_B2 += _r12;
        double _r13 = _grad7;
        double _r14 = _r13 * _t13;
        _d_nomGamma1_B2 += _r14;
        double _r15 = _t14 * _r13;
        * _d_gamma1_B2 += _r15;
        double _grad8 = 0.;
        double _grad9 = 0.;
        double _grad10 = 0.;
        gauss_pullback(_t16, _t17, 0.10000000000000001, _d_constraint[4], &_grad8, &_grad9, &_grad10);
        double _r16 = _grad8;
        * _d_lumi += _r16;
        double _r17 = _grad9;
        _d_nominalLumi += _r17;
        double _r18 = _grad10;
    }
}
