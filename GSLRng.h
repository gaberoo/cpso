#ifndef _GSLRNG_H
#define _GSLRNG_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <vector>
using namespace std;

namespace myGSL {
  class Rng;
  class PoissonGenerator;
};

typedef myGSL::Rng GSLRng;

class myGSL::Rng {
public:
  Rng(const gsl_rng_type * T = gsl_rng_taus2, unsigned long seed = time(NULL)) 
    : rng(gsl_rng_alloc(T)), discrete_t(NULL)
  { gsl_rng_set(rng,seed); }

  Rng(const Rng& R)
    : rng(gsl_rng_clone(R.constptr())), discrete_t(NULL)
  {}
  
  ~Rng() { 
    gsl_rng_free(rng);
    if (discrete_t != NULL) discrete_free();
  }

  inline void set(unsigned long seed) { gsl_rng_set(rng,seed); }
  inline unsigned long min() { return gsl_rng_min(rng); }
  inline unsigned long max() { return gsl_rng_max(rng); }

  // uniform distributions
  inline unsigned long get() { return gsl_rng_get(rng); }
  inline double uniform() { return gsl_rng_uniform(rng); }
  inline double uniform_pos() { return gsl_rng_uniform_pos(rng); }
  inline double uniform(double a, double b) { return a + (b-a)*uniform(); }
  inline unsigned long uniform_int(unsigned long n) { return gsl_rng_uniform_int(rng,n); }

  // other distributions
  inline double gaussian(double mu, double s) { return gsl_ran_gaussian(rng,s) + mu; }
  inline unsigned long binomial(double r, unsigned n) { return gsl_ran_binomial(rng,r,n); }
  inline unsigned int poisson(double mu) { return gsl_ran_poisson(rng,mu); }
  inline double exponential(double mu) { return gsl_ran_exponential(rng,mu); }

  // discrete distributions
  inline void discrete_preproc(int K, const double * P) { discrete_t = gsl_ran_discrete_preproc(K,P); }
  inline int  get_discrete() { return gsl_ran_discrete(rng,discrete_t); }
  inline void discrete_free() { gsl_ran_discrete_free(discrete_t); discrete_t = NULL; }

  template<typename T> inline void shuffle(vector<T>& src, size_t first = 0, 
      size_t last = -1) {
    if (first >= src.size()) return;
    if (last > src.size()) last = src.size();
    gsl_ran_shuffle(rng,&src[first],last-first,sizeof(T));
  }

  template<typename T> inline void choose(vector<T>& src, vector<T>& dest) {
    gsl_ran_choose(rng,&dest[0],dest.size(),&src[0],src.size(),sizeof(T));
  }

  inline gsl_rng * ptr() { return rng; }
  inline const gsl_rng * constptr() const { return rng; }

private:
  gsl_rng * rng;
  gsl_ran_discrete_t * discrete_t;
};

class myGSL::PoissonGenerator {
  public:
    double mu;
    Rng* rng;
    PoissonGenerator(double m, Rng* r) { mu = m; rng = r; }
    double operator()() { return rng->poisson(mu); }
};

#endif
