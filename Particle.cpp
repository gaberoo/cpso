#include "Particle.h"

std::ostream& PSO::operator<<(std::ostream& out, const PSO::Particle& p) {
  out << p.value << " " << p.position << " " << p.velocity << " ";
  return out;
}

