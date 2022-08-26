#include <iostream>
#include <string>
#include <unordered_set>

class contextManager {
  // Stores a set of decl'd funcs for reuse
  std::unordered_set<std::string> decldFuncs;
  // Variable to keep track of loops over data.
  bool isInsideLoop = false;

public:
  std::vector<std::string> inputParams;

  bool isFuncDecld(std::string funcName) {
    return decldFuncs.find(funcName) != decldFuncs.end();
  }

  void declFunc(std::string funcName) { decldFuncs.insert(funcName); }

  void inLoop() { isInsideLoop = true; }
  bool isInLoop() { return isInsideLoop; }

  void addInputParam(std::string name) { inputParams.push_back(name); }

  std::string getParamList() { 
    std::string params = "";
    for(auto& it : inputParams) {
      params += "double " + it + ",";
    }
    params.pop_back();
    return params;
  }

  std::string getNextInputIdx() { return std::to_string(inputParams.size() - 1); }
};

class ExRooReal {
  std::vector<ExRooReal*> child;
  std::vector<std::string> preFuncDecls;
  
protected:
  contextManager &ctxM;
  std::string result;

public:

  ExRooReal() = delete;
  ExRooReal(contextManager &ctx, std::vector<ExRooReal*> ch = {}) : child(ch), ctxM(ctx) {}
  virtual ~ExRooReal() = default;
  
  void setChildren(std::vector<ExRooReal *> ch) { child = ch;  }
  
  std::string getCode() {
    std::string code = "", global = "";
    getCodeRecur(this, code, global);
    return global + code;
  }
  void getCodeRecur(ExRooReal *head, std::string &code, std::string &global) {
    bool ILP = head->isLoopProducing();
    if (ILP)
      code += head->buildLoopBegin(global);
    for (auto it : head->child) {
      getCodeRecur(it, code, global);
    }
    code += head->translate(global, preFuncDecls);
    if (ILP)
      code += head->buildLoopEnd(global);
  }
  virtual std::string translate(std::string &globalScope,
                                std::vector<std::string> &preFuncDecls) = 0;
  virtual std::string getResult() { return result; }
  void updateResults(std::string newRes) { result = newRes; }
  virtual bool isLoopProducing() { return false; }
  virtual std::string buildLoopBegin(std::string &globalScope) { return ""; };
  virtual std::string buildLoopEnd(std::string &globalScope) { return ""; };
};

class ExRooRealVar : public ExRooReal {
  std::string obj;
  std::string def;
  std::string binVar = "";
  std::string init = "";
  int ele = 0;
  bool initDef = true;

public:
  ExRooRealVar(contextManager &ctxMgr, std::string obj_name, std::string initialize = "", int elements = 0)
      : ExRooReal(ctxMgr), obj(obj_name), init(initialize), ele(elements) {}
  // This is obviously a leaf node.
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    if (!initDef) return "";

    if (ele) {
      globalScope += "double " + obj + "[" + std::to_string(ele) + "]" + init + ";\n";
    } else if (!ele && init != "") {
      globalScope += "double " + obj + " = " + init + ";\n";
    } else {
      ctxM.addInputParam(obj);
      obj = "in[" + ctxM.getNextInputIdx() + "]";
    }

    initDef = false;
    return "";
  }

  void setBinVar(std::string bv) { binVar = bv; }
  std::string getBinVar() { return binVar; }

  std::string getResult() override {
    if (ele)
      return obj + "[" + binVar + "]";
    else
      return obj;
  }
};

class ExRooConst : public ExRooReal {
  double const_val;

public:
  ExRooConst(contextManager &ctxMgr, double val)
      : ExRooReal(ctxMgr), const_val(val) {}
  // This is obviously a leaf node.
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    result = std::to_string(const_val);
    return "";
  }
};

// Class to maintain the histogram?
class ExRooHistFunc : public ExRooReal {
  ExRooRealVar *x;
  inline static unsigned int _ObjCount = 0;
  unsigned int mycount;
  std::size_t nBins = 2;
  std::string name;
  std::string initalizer;
  bool unroll;
  std::vector<double> binBoundaries;
  std::vector<double> shapes;
  std::string binInit;
  std::string idx;
  bool init = true;

public:
  ~ExRooHistFunc() { _ObjCount--; }
  ExRooHistFunc(contextManager &ctxMgr, ExRooRealVar *in, std::string shapeName,
                std::string val_init, std::string bin_init = "")
      : ExRooReal(ctxMgr, {in}), x(in), name(shapeName), 
        initalizer(val_init), binInit(bin_init) {
    _ObjCount++;
    mycount = _ObjCount;
  }
  // Possibly also leaf-like (i.e. returns a scalar)
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    globalScope += "double " + name + "[" + std::to_string(nBins) + "]" +
                   initalizer + ";\n";
    globalScope += "double binBoundaries" + getCount() + "[" + std::to_string(nBins + 1) + "]" + binInit + ";\n";
    idx = "b" + getCount();
    return "unsigned int " + idx +
          " = ExRooHistFunc::getBin(binBoundaries" + getCount() + ", " +
          x->getResult() + ");\n";
    ;
  }

  std::string getResult() override {
    return name + "[" + idx + "]";
  }

  static int getBin(double* binBoundaries, double x) {
    int ibin = 0;
    while (binBoundaries[ibin + 1] < x) {
      ibin++;
    }
    return ibin;
  }

  std::string getCount() { return std::to_string(mycount); }
};

// Returns values dependant on the bin of x
class ExRooParamHistFunc : public ExRooReal {
  ExRooRealVar *x;
  std::vector<ExRooReal *> values;
  bool init = true;

public:
  ExRooParamHistFunc(contextManager &ctx, ExRooRealVar *in, std::vector<ExRooReal *> val) : ExRooReal(ctx), x(in), values(val) {
    val.push_back(in);
    setChildren(val);
  }

  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    if (init) {
      std::string decl =
          "double histVals[" + std::to_string(values.size()) + "]{";
      int idx = 0;
      for (auto &it : values) {
        decl += it->getResult() + ",";
        it->updateResults("histVals[" + std::to_string(idx) + "]");
        idx++;
      }
      decl[decl.size() - 1] = '}';
      init = false;
      globalScope += decl + ";\n";
    }
    result = "histVals[" + x->getBinVar() + "]";

    return "";
  }
};

class ExRooRealSum : public ExRooReal {
  std::string sum_val;
  std::vector<ExRooReal *> funcs;
  std::vector<ExRooReal *> coeffs;

public:
  ExRooRealSum(contextManager &ctxMgr, std::string sum,
               std::vector<ExRooReal *> func, std::vector<ExRooReal *> coeff)
      : ExRooReal(ctxMgr), sum_val(sum), funcs(func), coeffs(coeff) {
    func.insert(func.end(), coeff.begin(), coeff.end());
    setChildren(func);
  }

  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    std::string code = "";
    // list<RooAbsReal> lt= {}
    code += "double " + sum_val + " = 0;\n";
    for (int i = 0; i < funcs.size(); i++) {
      code += sum_val + " += " + funcs[i]->getResult() + " * " +
              coeffs[i]->getResult() + ";\n";
    }
    result = sum_val;
    return code;
  }
  void AddFunc(ExRooReal *obj) { funcs.push_back(obj); }
  void AddCoeff(ExRooReal *obj) { coeffs.push_back(obj); }
};

class ExRooProduct : public ExRooReal {
  std::vector<ExRooReal *> prods;

public:
  ExRooProduct(contextManager &ctxMgr, std::vector<ExRooReal *> prod)
      : ExRooReal(ctxMgr, prod), prods(prod) {}
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    return "";
  }

  std::string getResult() override {
    if (result == "") {
      result = '(' + prods[0]->getResult();
      for (int i = 1; i < prods.size(); i++) {
        result += "*" + prods[i]->getResult();
      }
      result += ')';
    }
    return result;
  }
};

class ExRooPoisson : public ExRooReal {
  ExRooReal *x;
  ExRooReal *par;

public:
  ExRooPoisson(contextManager &ctxMgr, ExRooReal *x, ExRooReal *par)
      : ExRooReal(ctxMgr, {x, par}), x(x), par(par) {}
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    result = "ExRooPoisson::poisson(" + x->getResult() + "," +
             par->getResult() + ")";
    return "";
  }

  static double poisson(double x, double par) {
    if (par < 0)
      return TMath::QuietNaN();
    if (x < 0)
      return 0;
    else if (x == 0.0)
      return TMath::Exp(-par);
    else {
      return TMath::Exp(x * log(par) - TMath::LnGamma(x + 1.) - par);
    }
  }
};

class ExRooGaussian : public ExRooReal {
  ExRooReal *_x;
  ExRooReal *_mu;
  ExRooReal *_sigma;

public:
  ExRooGaussian(contextManager &ctx, ExRooReal *in, ExRooReal *mean,
                ExRooReal *sig)
      : ExRooReal(ctx, {in, mean, sig}), _x(in), _mu(mean), _sigma(sig) {}
      
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    result = "ExRooGaussian::gauss(" + _x->getResult() + " ," + _mu->getResult() +
             " ," + _sigma->getResult() + ")";
    return "";
  }

  static double gauss(double x, double mean, double sigma) {
    const double arg = x - mean;
    const double sig = sigma;
    double out = std::exp(-0.5 * arg * arg / (sig * sig));
    return 1. / (std::sqrt(TMath::TwoPi()) * sigma) * out;
  }
};

class ExRooConstraintSum : public ExRooReal {
  std::vector<ExRooReal *> constraints;
  bool init = true;

public:
  ExRooConstraintSum(contextManager &ctx, std::vector<ExRooReal*> elements) : ExRooReal(ctx, elements), constraints(elements) {}
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    if (init) {
      std::string decl =
          "double constraint[" + std::to_string(constraints.size()) + "]{";
      int idx = 0;
      for (auto &it : constraints) {
        decl += it->getResult() + ",";
        it->updateResults("constraints[" + std::to_string(idx) + "]");
        idx++;
      }
      decl[decl.size() - 1] = '}';
      init = false;
      globalScope += decl + ";\n";
    }

    globalScope += "double cnstSum = 0;\n";
    std::string code = "for(int i = 0; i < " +
                       std::to_string(constraints.size()) + "; i++) \n{ \
        cnstSum -= std::log(constraint[i]); \n}\n";
    result = "cnstSum";
    return code;
  }
};

class ExRooAddition : public ExRooReal {
  std::vector<ExRooReal *> elements;
  bool init = true;

public:
  ExRooAddition(contextManager &ctx, std::vector<ExRooReal*> in) : ExRooReal(ctx, in), elements(in) {}
  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    if (elements.size() > 3) {
      if (init) {
        std::string decl =
            "double elements[" + std::to_string(elements.size()) + "]{";
        int idx = 0;
        for (auto &it : elements) {
          decl += it->getResult() + ",";
          it->updateResults("elements[" + std::to_string(idx) + "]");
          idx++;
        }
        decl[decl.size() - 1] = '}';
        init = false;
        globalScope += decl + ";\n";
      }

      globalScope += "double eleSum = 0;";
      std::string code = "for(int i = 0; i < " +
                         std::to_string(elements.size()) + "; i++) {\n \
          eleSum += elements[i]; \n \
        }\n";
      result = "eleSum";
      return code;
    }

    for (auto &it : elements) {
      result += it->getResult() + '+';
    }
    result.pop_back();

    return "";
  }
};

class KSum {
  double c = 0;
  double y = 0;
  double t = 0;

public:
  double accumulate(double sum, double a) { 
    y = - a - c;
    t = sum + y;
    c = (t - sum) - y;
    sum = t;
    return sum;
  }
};

class ExRooNll : public ExRooReal {
  ExRooRealVar *x;
  ExRooReal *func;
  std::string res = "nllSum";
  std::string kahnSum = "ksum";
  std::string bins;

public:
  ExRooNll(contextManager &ctx, ExRooRealVar *in, std::string binNum,
           ExRooReal *funcProd)
      : ExRooReal(ctx, {in, funcProd}), bins(binNum), x(in), func(funcProd) {}
  bool isLoopProducing() override { return true; }
  std::string buildLoopBegin(std::string &globalScope) override {
    std::string code = "double " + res + " = 0;\n";
    result = res;
    // code += "KSum " + kahnSum + ";\n";
    code += "for(int iB = 0; iB < " + bins + "; iB++) {\n";
    x->setBinVar("iB");
    return code;
  }

  std::string buildLoopEnd(std::string &globalScope) override { return "}\n"; }

  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    std::string code = "double temp;\n";
    code += "temp = std::log(" + func->getResult() + ");\n";
    // code += res + " = " + kahnSum + ".accumulate(" + res + ", temp);\n";
    code += res + " -= temp;\n";
    return code;
  }
};

class ExRooNll2 : public ExRooReal {
  ExRooRealVar *x;
  ExRooReal *func;
  std::string res = "nllSum";
  std::string kahnSum = "ksum";
  std::string bins;
  std::string idx;
  std::vector<double> weights;

public:
  ExRooNll2(contextManager &ctx, ExRooRealVar *in, std::string binNum,
           ExRooReal *funcProd,std::vector<double> wgts)
      : ExRooReal(ctx, {in, funcProd}), bins(binNum), x(in), func(funcProd), weights(wgts) {}
  bool isLoopProducing() override { return true; }
  std::string buildLoopBegin(std::string &globalScope) override {
    std::string code = "double " + res + " = 0;\n";
    result = res;
    // code += "KSum " + kahnSum + ";\n";
    code += "for(int iB = 0; iB < " + bins + "; iB++) {\n";
    x->setBinVar("iB");
    idx = "iB";
    return code;
  }

  std::string buildLoopEnd(std::string &globalScope) override { return "}\n"; }

  std::string translate(std::string &globalScope,
                        std::vector<std::string> &preFuncDecls) override {
    std::string decl = "double weights[" + std::to_string(weights.size()) + "]{";
    for(auto& it : weights) {
      decl += std::to_string(it) + ",";
    }
    decl.pop_back();
    decl += "};\n";
    globalScope += decl;
    std::string code = "double temp;\n";
    code += "temp = std::log(" + func->getResult() + ");\n";
    // code += res + " = " + kahnSum + ".accumulate(" + res + ", temp);\n";
    code += res + " -= -" + func->getResult() + " + weights[" + idx + "] * temp;\n";
    return code;
  }
};
