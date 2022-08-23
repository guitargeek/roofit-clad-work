#include <iostream>
#include <vector>

class ExRooReal;

class NTree {
    public: 
    std::vector<ExRooReal*> child;
    ExRooReal* data;
    std::vector<std::string> preFuncDecls;

    NTree(ExRooReal* obj = nullptr, std::vector<ExRooReal*> ch = {}) : data(obj), child(ch) { }

    void setChildren(std::vector<ExRooReal *> ch) { child = ch;  }
    std::string getCode() {
      std::string code = "", global = "";
      getCodeRecur(data, code, global);
      return global + code;
    }

    void getCodeRecur(ExRooReal *head, std::string &code, std::string &global);
};
