#include "Point.h"

std::ostream& PSO::operator<<(std::ostream& out, const Point& p) {
  for (size_t i(0); i < p.size(); ++i) {
    out << p[i] << " ";
  }
  return out;
}

PSO::Point operator+(const PSO::Point& p1, const PSO::Point& p2) {
  PSO::Point p(p1.size());
  if (p1.size() != p2.size()) return p;
  for (size_t i(0); i < p1.size(); ++i) p[i] = p1[i]+p2[i];
  return p;
}

PSO::Point operator-(const PSO::Point& p1, const PSO::Point& p2) {
  PSO::Point p(p1.size());
  if (p1.size() != p2.size()) return p;
  for (size_t i(0); i < p1.size(); ++i) p[i] = p1[i]-p2[i];
  return p;
}

PSO::Point PSO::Point::operator*(double a) {
  PSO::Point p(size());
  for (size_t i(0); i < size(); ++i) p[i] = a*at(i);
  return p;
}

