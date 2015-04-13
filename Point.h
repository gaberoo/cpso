#ifndef __PSO_POINT__
#define __PSO_POINT__

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

namespace PSO {
  using namespace std;

  class Point : public vector<double> {
    public:
      Point() {}
      Point(int dim, double val = 0.0) : vector<double>(dim,val) {}
      Point(const Point& p) : vector<double>(p) {}
      Point(const vector<double>& x) : vector<double>(x) {}

      inline double norm2() const {
        double n(0.0);
        for (size_t i(0); i < size(); ++i) n += at(i)*at(i);
        return sqrt(n);
      }

      inline double euclid_dist(const Point& p) const {
        double d(0.0);
        double x(0.0);
        for (size_t i(0); i < size(); ++i) {
          x = at(i)-p.at(i);
          d += x*x;
        }
        return sqrt(d);
      }

      inline Point& operator=(const vector<double>& x) {
        if (size() == x.size()) assign(x.begin(),x.end());
        return *this;
      }

      friend Point operator+(const Point& p1, const Point& p2);
      friend Point operator-(const Point& p1, const Point& p2);
      Point operator*(double a);
      friend ostream& operator<<(ostream& out, const Point& p);
  };
  ostream& operator<<(ostream& out, const Point& p);
  Point operator+(const Point& p1, const Point& p2);
  Point operator-(const Point& p1, const Point& p2);
}

#endif
