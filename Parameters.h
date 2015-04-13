#ifndef __PSO_PARAMETERS__
#define __PSO_PARAMETERS__

#include <vector>

namespace PSO {
  using namespace std;
  enum VarType { REAL, INTEGER };

  class Parameters {
    friend class Swarm;

    public:
      Parameters() : numParams(0) {}
      Parameters(int np, void* ep = NULL) 
        : numParams(np), evalParams(ep)
      {
        lb = vector<double>(np,0.0);
        ub = vector<double>(np,0.0);
        type = vector<VarType>(np,REAL);
        locked = vector<bool>(np,false);
      }

      int numParams;                            /* number of parameters */
      vector<double> lb;                        /* lower bound */
      vector<double> ub;                        /* upper bound */
      vector<VarType> type;                     /* variable types */
      vector<bool> locked;                      /* lock variables */

      double ssd;                               /* search space diameter */
      double stag_t;                            /* stagnation threshold */

      double (*evalFunc)(const Point& x, void* params);
      void* evalParams;

      inline void lockVar(int i, double x) {
        locked[i] = true;
        lb[i] = x;
        ub[i] = x;
      }
  };
}

#endif
