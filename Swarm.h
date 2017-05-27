/* 
 * PARTICLE SWARM OPTIMIZATION
 * ===========================
 *
 * Author: Gabriel E Leventhal, 2011-12-02
 *
 * adapted from 
 *   Author: Tomas V. Arredondo
 *   SimPSOLib: A simple yet flexible PSO implementation in C++.
 *
 */

#ifndef __PSO_SWARM__
#define __PSO_SWARM__

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <ctime>
#include <queue>
using namespace std;

#include "GSLRng.h"
#include <gsl/gsl_math.h>

#include "Particle.h"
#include "Point.h"
#include "Parameters.h"
#include "Network.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace PSO {
  class Swarm {
    public:
      Swarm() {}

      Swarm(int size, int np, Parameters* pars) 
        : swarmSize(size), numParams(np), numEvals(0), curIt(0),
          numInform(3), 
          phiPi(np,2.5), phiPf(np,0.5), phiGi(np,0.5),phiGf(np,2.5),
          omegai(np,0.721), omegaf(np,0.721), p(pars),
          bestVal(-INFINITY), bestParticle(-1),
          mpi_type(0), mpi_rank(0), mpi_ntasks(1)
      {
        bestPos = Point(np,0.0);
        setInformants(numInform);
      }

      virtual ~Swarm() { destroy(); }

      inline void setMPI(int type) { mpi_type = (type >= 0) ? type : 0; }

      inline void setVars(const Point& pp, const Point& pg, const Point& o) {
        phiPi = pp; phiGi = pg; omegai = o;
        phiPf = pp; phiGf = pg; omegaf = o;
        phiP = pp; phiG = pg; omega = o;
      }

      inline void setInformants(int K) {
        numInform = (K > swarmSize-1) ? swarmSize-1 : K;
        links.resize(swarmSize);
        init_network();
      }

      double maximum_swarm_radius() const;
      double search_space_diameter() const;

      // get everything ready to go
      void begin(int argc, char* argv[]);
      // finish up stuff
      void end();

      void initialize(int type = 0);
      void display(ostream* out = &cout, string prefix = "");

      virtual void evaluate(int vflag = 0, ostream* out = NULL, ostream* hist = NULL);
      virtual double evaluateParticle(int j);

      void updateVelocity(int j, bool bestVals = true);
      void updateVelocity(bool bestVals = true);

      void updatePosition(int j);
      void updatePosition();

      void randomizeInf(int maxTries = 1000);
      void randomPosition(int j);

      virtual void run(int numIt, int slowdown = 0, int vflag = 0, 
                       ostream* out = NULL, ostream* hist = NULL);

      inline bool is_master() const { return (mpi_rank == 0); }

      inline void seed_rng(unsigned long seed) { rng.set(seed); }

      inline const Particle& at(size_t i) const { return *swarm[i]; }
      inline Particle best() const { return *swarm[bestParticle]; }

      inline double fitness(int i) const { return swarm.at(i)->fitness(); }

      inline int size() const { return swarmSize; }

      inline int nextIt() { return ++curIt; }
      inline void setIt(int i) { curIt = i; }

      double mean(int i) const;
      double var(int i) const;
      Point mean() const;
      Point var() const;

#ifdef USE_MPI
      virtual void evaluate_slave();
      void run_mpi(int numInt, int vflag = 0, ostream* out = NULL, ostream* hist = NULL);
      virtual void run_master(int numInt, int vflag, ostream* out, ostream* hist);
      inline void evaluate_master(int vflag = 0) { run_master(-1,vflag,NULL,NULL); }
#endif

    protected:
      void create();
      void destroy();

      void init_network();
      int find_best();

      int swarmSize;                            /* swarm size */
      vector<Particle*> swarm;                  /* swarm */
      int numParams;                            /* number of parameters */
      int numEvals;                             /* current number of function evaluations */
      int curIt;                                /* current swarm iteration */

      int numInform;                            /* number of informants */
      Network links;                            /* information network */

      Point phiP;
      Point phiG;
      Point omega;
      Point phiPi;                              /* phiP initial */
      Point phiPf;                              /* phiP final */
      Point phiGi;                              /* phiG initial */
      Point phiGf;                              /* phiG final */
      Point omegai;                             /* omega initial */
      Point omegaf;                             /* omega final */

      Parameters* p;

      double bestVal;                           /* current best value */
      Point bestPos;                            /* current best position */
      int bestParticle;

      myGSL::Rng rng;

      int mpi_type;
      int mpi_rank;
      int mpi_ntasks;
  };
}

#endif // __PSO_SWARM__
