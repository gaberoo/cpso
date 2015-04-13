#include <iostream>
using namespace std;

// #include "EvalFunc.h"
#include "Swarm.h"
using namespace PSO;

double eval(const Point& x, void* params) {
  double f;
  f = -(x[0]-10)*(x[0]-10) - (x[1]-5)*(x[1]-5) + 10;
  return f;
}

int main(int argc, char* argv[]) {
  Parameters pars(2);
  pars.lb[0] = 0.0;
  pars.lb[1] = 0.0;
  pars.ub[0] = 20.0;
  pars.ub[1] = 10.0;
  pars.type[0] = REAL;
  pars.type[1] = REAL;
  pars.evalFunc = &eval;

  Swarm s(20,2,&pars);

  Point phi_p(2,1.49445);
  Point phi_g(2,1.49445);
  Point omega(2,0.729);
  s.setVars(phi_p,phi_g,omega);

  Point init(2);
  init[0] = 4.0;
  init[1] = 9.0;

  s.begin(argc,argv);
  s.create();
  s.initialize(init,0.4);
//  s.evaluate();
//  s.display();
  s.run_mpi(1000,1);
  s.end();

  return 0;
}
