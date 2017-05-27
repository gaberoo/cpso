#ifndef __PSO_PARAMETERS__
#define __PSO_PARAMETERS__

#include <vector>
#include "rapidjson/document.h"

#include "Point.h"

namespace PSO {
  using namespace std;
  enum VarType { REAL, INTEGER };

  class Parameters {
    friend class Swarm;

    public:
      Parameters() : numParams(0) {}
      Parameters(int np, void* ep = NULL) 
        : numParams(np), evalParams(ep)
      { resize(np); }

      void from_json(rapidjson::Value& v);

      int numParams;                            /* number of parameters */
      vector<double> lb;                        /* lower bound */
      vector<double> ub;                        /* upper bound */
      vector<VarType> type;                     /* variable types */
      vector<bool> locked;                      /* lock variables */
      vector<char> scale;                       /* variable scaling */

      double ssd;                               /* search space diameter */
      double stag_t;                            /* stagnation threshold */

      double (*evalFunc)(int gen, int id, const Point& x, const void* params);
      void* evalParams;

      inline void resize(size_t n) {
        numParams = n;
        lb.resize(n,0.0);
        ub.resize(n,0.0);
        type.resize(n,REAL);
        locked.resize(n,false);
        scale.resize(n,'n');
      }

      inline void lockVar(int i, double x) {
        locked[i] = true;
        lb[i] = x;
        ub[i] = x;
      }

      inline double initVar(int i) const {
        if (i >= 0 && i < numParams) {
          return 0.5*(lb[i]+ub[i]);
        } else {
          return 0.0;
        }
      }
  };
}

inline void PSO::Parameters::from_json(rapidjson::Value& jpars) {
  rapidjson::Value::MemberIterator m1;
  rapidjson::Value::MemberIterator it = jpars.FindMember("pars");
  if (it == jpars.MemberEnd()) throw "pars";

  rapidjson::Value& d = it->value;

  resize(d.Size());

  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    type[i] = REAL;

    m1 = d[i].FindMember("limits"); 
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsArray()) throw "Bad variables limits.";
      if (! (m1->value.Size() == 2)) throw "Bad variables limits.";
      rapidjson::SizeType j;
      j = 0; lb[i] = m1->value[j].GetDouble();
      j = 1; ub[i] = m1->value[j].GetDouble();
    }

    m1 = d[i].FindMember("lock");
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsDouble()) {
        throw "Locked variable isn't a double.";
      } else {
        lockVar(i,m1->value.GetDouble());
      }
    }

    m1 = d[i].FindMember("scale");
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsString()) throw "Bad scale.";
      scale[i] = m1->value.GetString()[0];
    } else {
      scale[i] = 'n';
    }
  }
}

#endif
