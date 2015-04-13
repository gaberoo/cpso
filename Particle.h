#ifndef __PSO_PARTICLE__
#define __PSO_PARTICLE__

#include <iostream>
#include <vector>
#include <cmath>

#include "Point.h"

namespace PSO {
  using namespace std;

  class Particle {
    public:
      friend class Swarm;

      Particle() {}

      Particle(int np) 
        : numParams(np), value(log(0.0)), bestValue(-INFINITY)
      {
        position.resize(np,0.0);
        velocity.resize(np,0.0);
        bestPosition.resize(np,0.0);
      }

      virtual ~Particle() {}

      friend ostream& operator<<(ostream& output, const Particle& p);

      inline void setPos(const Point& x) { position = x; }
      inline void setPos(int i, double xi) { position[i] = xi; }

      inline double& x(int i) { return position[i]; }
      inline double& v(int i) { return velocity[i]; }

      inline double operator[](int i) const { return position[i]; }
      inline double fitness() const { return value; }

      inline const double* pos() const { return position.data(); }

      int numParams;
      double value;
      double bestValue;
      Point position;
      Point velocity;
      Point bestPosition;
  };

  ostream& operator<<(ostream& output, const Particle& p);
}

#endif // __PSO_PARTICLE__
