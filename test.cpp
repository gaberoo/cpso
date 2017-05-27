#include <iostream>
#include <fstream>
using namespace std;

// #include "EvalFunc.h"
#include "Swarm.h"
using namespace PSO;

double eval(int gen, int id, const Point& x, const void* params) {
  double f;
  f = -(x[0]-10)*(x[0]-10) - (x[1]-5)*(x[1]-5) + 10;
  return f;
}

int main(int argc, char* argv[]) {
  string json_input;
  ifstream in(argv[1]);
  in.seekg(0,ios::end);
  json_input.reserve(in.tellg());
  in.seekg(0,ios::beg);
  json_input.assign(istreambuf_iterator<char>(in), 
                    istreambuf_iterator<char>());

  /* Reading priors from JSON */
  rapidjson::Document jpars;
  try {
    jpars.Parse<0>(json_input.c_str());
    if (! jpars.IsObject()) throw 10;
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl;
    cerr << json_input << endl;
  }

  cerr << "Reading parameters..." << flush;

  Parameters pars;
  try {
    pars.from_json(jpars);
  } catch (char* str) {
    cerr << "JSON error: " << str << endl;
  }

  cerr << "done" << endl;

  pars.evalFunc = &eval;

  Swarm s(20,2,&pars);

  Point phi_p(2,1.49445);
  Point phi_g(2,1.49445);
  Point omega(2,0.729);
  s.setVars(phi_p,phi_g,omega);

  s.begin(argc,argv);
  s.initialize(1);
  s.evaluate();
//  s.display();
//  s.run_mpi(1000,1);
  s.run(1000000,0,1,NULL,NULL);
  s.end();

  return 0;
}
