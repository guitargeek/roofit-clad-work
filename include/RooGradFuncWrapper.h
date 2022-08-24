#ifndef RooGradFuncWrapper_h
#define RooGradFuncWrapper_h

template <typename Func = void, typename Grad = void>
class RooGradFuncWrapper final : public RooAbsReal {
public:
   RooGradFuncWrapper(Func function, Grad gradient, RooWorkspace const &ws, std::vector<std::string> const &paramNames)
      : RooAbsReal{"RooGradFuncWrapper", "RooGradFuncWrapper"}, _paramProxies{"!params", "List of parameters", this},
        _func(function), _grad(gradient)
   {
      for (auto const &paramName : paramNames) {
         auto *var = ws.var(paramName.c_str());
         _paramProxies.add(*var);
         _params.emplace_back(var);
      }
   }

   RooGradFuncWrapper(const RooGradFuncWrapper &other, const char *name = nullptr)
      : RooAbsReal(other, name), _paramProxies("!params", this, other._paramProxies), _func(other._func),
        _grad(other._grad)
   {
   }

   TObject *clone(const char *newname) const override { return new RooGradFuncWrapper(*this, newname); }

   double defaultErrorLevel() const override { return 0.5; }

protected:
   double evaluate() const override
   {
      unsigned int n = _paramProxies.size();
      unsigned int idx = 0;
      double funcParams[n];
      for (auto const &param : _params) {
         funcParams[idx++] = param->getValV();
      }

      return _func(funcParams);
   }

   void evaluateGradient(double *out) const override
   {
      unsigned int n = _paramProxies.size();
      double funcParams[n];
      unsigned int idx = 0;
      for (auto const &param : _params) {
         funcParams[idx] = param->getValV();
         out[idx] = 0;
         idx++;
      }

      clad::array_ref<double> dx(out, n);

      _grad(funcParams, dx);
   }

private:
   RooListProxy _paramProxies;
   std::vector<RooRealVar *> _params;
   Func _func;
   Grad _grad;
};

#endif
