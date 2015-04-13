#ifndef __PSO_EVALFUNC__
#define __PSO_EVALFUNC__

#include <vector>

namespace PSO {
  using namespace std;

  enum VarType { REAL, INTEGER };

  class EvalFunc {
    friend class Swarm;

    public:
      EvalFunc() : numParams(0) {}
      EvalFunc(int np) : numParams(np) {
        lb = vector<double>(np,0.0);
        ub = vector<double>(np,0.0);
        type = vector<VarType>(np,REAL);
        locked = vector<bool>(np,false);
      }
      EvalFunc(const EvalFunc& ef)
        : numParams(ef.numParams), lb(ef.lb), ub(ef.ub),
          type(ef.type), locked(ef.locked) {}

      int numParams;                            /* number of parameters */
      vector<double> lb;                        /* lower bound */
      vector<double> ub;                        /* upper bound */
      vector<VarType> type;                     /* variable types */
      vector<bool> locked;                      /* lock variables */

      virtual double eval(const vector<double>& x) const = 0;

      inline void lockVar(int i, double x) {
        locked[i] = true;
        lb[i] = x;
        ub[i] = x;
      }

      static void force_quit() { PSO::EvalFunc::_force_quit = true; }
      static bool get_force_quit() { return PSO::EvalFunc::_force_quit; }

    private:
      static bool _force_quit;
  };
}

#endif // __PSO_EVALFUNC__

