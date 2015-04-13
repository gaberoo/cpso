#include "Swarm.h"

void PSO::Swarm::begin(int argc, char* argv[]) {
#ifdef USE_MPI
  MPI::Init(argc,argv);
  mpi_rank = MPI::COMM_WORLD.Get_rank();
  mpi_ntasks = MPI::COMM_WORLD.Get_size();
  // fprintf(stderr,"Hello, I'm process %d.\n",mpi_rank);
#endif

  create();
}

// ===========================================================================

void PSO::Swarm::end() {
  destroy();
#ifdef USE_MPI
  MPI::Finalize();
#endif
}

// ===========================================================================

void PSO::Swarm::create() {
#ifdef USE_MPI
  if (mpi_rank != 0) return;
#endif
  if (swarm.size() > 0) destroy();
  swarm.resize(swarmSize,NULL);
  vector<Particle*>::iterator it;
  for (it = swarm.begin(); it != swarm.end(); ++it) {
    *it = new Particle(numParams);
  }
}

// ===========================================================================

void PSO::Swarm::evaluate(int vflag, ostream* out, ostream* hist) {
#ifdef USE_MPI
  if (mpi_rank == 0) run_master(-1,vflag,out,hist);
  else evaluate_slave();
  MPI::COMM_WORLD.Barrier();

#else
  // double f(-INFINITY);
  for (size_t j = 0; j < swarm.size(); ++j) {
    /*f = */ evaluateParticle(j);
  }
#endif
}

// ===========================================================================

double PSO::Swarm::evaluateParticle(int j) {
  double f(-INFINITY);
  f = p->evalFunc(swarm[j]->position,p->evalParams);
  swarm[j]->value = f;
  if (f >= swarm[j]->bestValue) {
    swarm[j]->bestPosition = swarm[j]->position;
    swarm[j]->bestValue = f;
  }
  ++numEvals;
  if (f >= bestVal) {
    bestPos = swarm[j]->position;
    bestVal = f;
    bestParticle = j;
  }
  return f;
}

// ===========================================================================

#ifdef USE_MPI
/*
void PSO::Swarm::evaluate_master(int vflag) {
  double f(log(0.0));
  int id(0);
  int j(0);
  double* particlePos(NULL);
  MPI::Status status;
  int flag;
  int src;
  int idle = 0;
  // initialize slaves
  for (int k(1); k < mpi_ntasks && j < swarm.size(); ++k) {
    particlePos = swarm[j]->position.data();
    MPI::COMM_WORLD.Send(&j,1,MPI::INT,k,1);
    MPI::COMM_WORLD.Send(particlePos,numParams,MPI::DOUBLE,k,1);
    if (vflag) fprintf(stderr,"Sending particle %d to process %d.\n",j,k);
    ++j;
  }
  while (1) {
//    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
    flag = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE,MPI::ANY_TAG,status);
    if (flag) {
      src = status.Get_source();
      MPI::COMM_WORLD.Recv(&id,1,MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
      MPI::COMM_WORLD.Recv(&f,1,MPI::DOUBLE,src,MPI::ANY_TAG,status);
      if (vflag) fprintf(stderr,"Receiving particle %d from process %d.\n",id,src);
      swarm[id]->value = f;
      if (f >= swarm[id]->bestValue) {
        swarm[id]->bestPosition = swarm[id]->position;
        swarm[id]->bestValue = f;
      }
      ++numEvals;
      if (f >= bestVal) {
        bestPos = swarm[id]->position;
        bestVal = f;
        bestParticle = id;
      }
      if (j < swarm.size()) {
        if (vflag) fprintf(stderr,"Sending particle %d to process %d.\n",j,src);
        particlePos = swarm[j]->position.data();
        MPI::COMM_WORLD.Send(&j,1,MPI::INT,src,1);
        MPI::COMM_WORLD.Send(particlePos,numParams,MPI::DOUBLE,src,1);
        ++j;
      } else {
        ++idle;
        if (vflag) fprintf(stderr,"Sending done signal to process %d.\n",src);
        MPI::COMM_WORLD.Send(0,0,MPI::INT,src,0);
      }
      if (idle == mpi_ntasks-1) break;
    }
  }
  if (vflag) fprintf(stderr,"done.\n");
}
*/

// ===========================================================================

void PSO::Swarm::evaluate_slave() {
  double f(log(0.0));
  int id(0);
  int flag(0);
  int tag(0);
  int dest(0);
  Point position(numParams);
  MPI::Status status;
//  fprintf(stderr,"Slave %d ready.\n",mpi_rank);
  while (1) {
//    flag = MPI::COMM_WORLD.Iprobe(0,MPI::ANY_TAG,status);
//    if (flag) {
//      tag = status.Get_tag();
    MPI::COMM_WORLD.Recv(&id,1,MPI::INT,0,MPI::ANY_TAG,status);
    if (status.Get_tag() == 0) break;
    MPI::COMM_WORLD.Recv(position.data(),numParams,MPI::DOUBLE,0,MPI::ANY_TAG,status);
    f = p->evalFunc(position,p->evalParams);
    MPI::COMM_WORLD.Send(&id,1,MPI::INT,0,2);
    MPI::COMM_WORLD.Send(&f,1,MPI::DOUBLE,0,2);
//    }
  }
//  fprintf(stderr,"Slave %d done.\n",mpi_rank);
}
#endif

// ===========================================================================

void PSO::Swarm::destroy() {
#ifdef USE_MPI
  if (mpi_rank != 0) return;
#endif
  vector<Particle*>::iterator it;
  for (it = swarm.begin(); it != swarm.end(); ++it) {
    delete *it;
    *it = NULL;
  }
  swarm.clear();
}

// ===========================================================================

void PSO::Swarm::initialize(int type) {
#ifdef USE_MPI
  if (mpi_rank != 0) return;
#endif
  int i, j;
  double randPos;
  phiP = phiPi;
  phiG = phiGi;
  omega = omegai;
  // Latin Hypercube sampling
  vector< vector<int> > samples;
  if (type == 1) {
    samples.resize(numParams,vector<int>(swarmSize,0));
    for (i = 0; i < numParams; ++i) {
      for (j = 0; j < swarmSize; ++j) samples[i][j] = j;
      rng.shuffle(samples[i]);
    }
  }
  for (j = 0; j < swarmSize; ++j) {
    for (i = 0; i < numParams; ++i) {
      if (p->locked[i]) randPos = p->lb[i];
      else {
        if (type == 1) {
          // Latin Hypercube sampling
          double intervalSize = (p->ub[i]-p->lb[i])/swarmSize;
          randPos = p->lb[i] + samples[i][j]*(p->ub[i]-p->lb[i])/swarmSize
                    + rng.uniform(0.0,intervalSize);
        } else {
          // Uniform sampling
          randPos = rng.uniform(p->lb[i],p->ub[i]);
        }
        if (p->type[i] == INTEGER) randPos = floor(randPos);
        if (randPos < p->lb[i]) randPos = p->lb[i];
        else if (randPos > p->ub[i]) randPos = p->ub[i];
      }
      swarm[j]->setPos(i,randPos);
    }
  }
}
// ===========================================================================

void PSO::Swarm::randomizeInf(int maxTries) {
  int cnt(0);
  double f(-INFINITY);
  for (int j = 0; j < swarmSize; ++j) {
    if (swarm[j]->value == -INFINITY) {
      cerr << "f(" << j << ") = Inf. Randomizing position..." << flush;
      cnt = 0;
      f = -INFINITY;
      while (cnt++ <= maxTries) {
        randomPosition(j);
        f = evaluateParticle(j);
        if (f > -INFINITY) break;
      }
      cerr << " new pos = " << *swarm[j] << endl;
    }
  }
}

// ===========================================================================

void PSO::Swarm::randomPosition(int j) {
  double randPos;
  for (int i(0); i < numParams; ++i) {
    if (p->locked[i]) randPos = p->lb[i];
    else {
      // Uniform sampling
      randPos = rng.uniform(p->lb[i],p->ub[i]);
      if (p->type[i] == INTEGER) randPos = floor(randPos);
      if (randPos < p->lb[i]) randPos = p->lb[i];
      else if (randPos > p->ub[i]) randPos = p->ub[i];
    }
    swarm[j]->setPos(i,randPos);
  }
}

// ===========================================================================

void PSO::Swarm::updateVelocity(int j) {
  double rp;
  double rg;
  double dg;
  double dp;
  int i;
  rp = rng.uniform();
  rg = rng.uniform();

  // get informants for particle
  double f(-INFINITY);
  int bestInfId(bestParticle);
  Particle* bestInformant(swarm[bestParticle]);

  if (numInform < swarmSize-1) {
    set<int>::const_iterator nn;
    for (nn = links.nn_begin(j); nn != links.nn_end(j); ++nn) {
      if (swarm[*nn]->bestValue > f) {
        bestInfId = *nn;
        bestInformant = swarm[*nn];
        f = bestInformant->bestValue;
      }
    }
  }


  // Shi, Eberhart (1998, 2001)
  for (i = 0; i < numParams; ++i) {
    if (! p->locked[i]) {
      switch (p->type[i]) {
        case INTEGER:
          dp = round(swarm[j]->bestPosition[i]) - round(swarm[j]->position[i]);
          // dg = round(bestPos[i]) - round(swarm[j]->position[i]);
          dg = round(bestInformant->bestPosition[i]) - round(swarm[j]->position[i]);
          break;
        case REAL:
        default:
          dp = swarm[j]->bestPosition[i] - swarm[j]->position[i];
          // dg = bestPos[i] - swarm[j]->position[i];
          dg = bestInformant->bestPosition[i] - swarm[j]->position[i];
          break;
      }

      swarm[j]->velocity[i] = omega[i]*(swarm[j]->velocity[i]);
      swarm[j]->velocity[i] += phiP[i]*rp*dp;
      if (bestInfId != j) swarm[j]->velocity[i] += phiG[i]*rg*dg;

      // Limit velocity in each dimension to full dynamic range 
      // in search space (Simplifying PSO, Pedersen 2009)
      if (fabs(swarm[j]->velocity[i]) > (p->ub[i] - p->lb[i])) {
        swarm[j]->velocity[i] = rng.uniform()*(p->ub[i] - p->lb[i]);
      }
    }
  }
}

// ===========================================================================

void PSO::Swarm::updateVelocity() {
  for (int j(0); j < swarmSize; ++j) {
    // if (j == bestParticle) continue;
    updateVelocity(j);
  }
}

// ===========================================================================

void PSO::Swarm::updatePosition(int j) {
  // Shi, Eberhart (1998, 2001)
  for (int i(0); i < numParams; ++i) {
    if (! p->locked[i]) {
      swarm[j]->position[i] = swarm[j]->position[i] + swarm[j]->velocity[i];
      // Limit velocity in each dimension to full dynamic range
      // in search space (Simplifying PSO, Pedersen 2009)
      if (swarm[j]->position[i] > (p->ub[i])) {
        swarm[j]->position[i] = p->ub[i];
        swarm[j]->velocity[i] = 0.0;
      }
      else if (swarm[j]->position[i] < (p->lb[i])) {
        swarm[j]->position[i] = p->lb[i];
        swarm[j]->velocity[i] = 0.0;
      }
    }
  }
}

// ===========================================================================

void PSO::Swarm::updatePosition() {
  for (int j(0); j < swarmSize; ++j) {
    // if (j == bestParticle) continue;
    updatePosition(j);
  }
}

// ===========================================================================

#ifdef USE_MPI
void PSO::Swarm::run_mpi(int numInt, int vflag, ostream* out, ostream* hist) {
  if (mpi_rank == 0) run_master(numInt,vflag,out,hist);
  else evaluate_slave();
  MPI::COMM_WORLD.Barrier();
}

// ===========================================================================

void PSO::Swarm::run_master(int numIt, int vflag, ostream* out, ostream* hist) {
  double f(-INFINITY);
  int id(0);
  int j(0);
  double* particlePos(NULL);
  MPI::Status status;
  int flag;
  int src;
  int idle(0);
  int iter(0);
  queue<int> evalQueue;

  for (int i(0); i < swarm.size(); ++i) evalQueue.push(i);

  if (vflag) cerr << "Sending particles to slaves..." << endl;
  // initialize slaves
  for (int k(1); k < mpi_ntasks && (iter < numIt || numIt < 0); ++k) {
    j = evalQueue.front();
    evalQueue.pop();
    if (numIt > 0) {
      // numIt < 0 => evaluate at current position
      updateVelocity(j);
      updatePosition(j);
    }
    if (vflag) cerr << j << " " << (*swarm[j]) << endl;
    if (vflag) fprintf(stderr,"Sending particle %d to process %d.\n",j,k);
    MPI::COMM_WORLD.Send(&j,1,MPI::INT,k,1);
    MPI::COMM_WORLD.Send(swarm[j]->position.data(),numParams,MPI::DOUBLE,k,1);
    ++iter;
  }

  while (1) {
    flag = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE,MPI::ANY_TAG,status);
    if (flag) {
      // get function value
      src = status.Get_source();
      MPI::COMM_WORLD.Recv(&id,1,MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
      if (vflag) fprintf(stderr,"Receiving particle %d from process %d.\n",id,src);
      MPI::COMM_WORLD.Recv(&f,1,MPI::DOUBLE,src,MPI::ANY_TAG,status);

      // update particle information
      swarm[id]->value = f;
      if (f >= swarm[id]->bestValue) {
        swarm[id]->bestPosition = swarm[id]->position;
        swarm[id]->bestValue = f;
      }
      ++numEvals;

      if (hist != NULL) {
         *hist << id << " " << (*swarm[id]) << endl;
      }

      // check for new best value
      if (f >= bestVal) {
        bestPos = swarm[id]->position;
        bestVal = f;
        bestParticle = id;
        if (out != NULL) {
          *out << numEvals << " " << bestVal << " ";
          for (int j(0); j < bestPos.size(); ++j) *out << bestPos[j] << " ";
          *out << endl;
        }
      }

      if (numIt > 0) {
        // update velocity and position
        updateVelocity(id);
        updatePosition(id);
        evalQueue.push(id);
      }

      // send new work to slave
      // if (iter < numIt) {
      if ((iter < numIt || numIt < 0) && ! evalQueue.empty()) {
        j = evalQueue.front();
        evalQueue.pop();
        if (vflag) fprintf(stderr,"Sending particle %d to process %d.\n",j,src);
        MPI::COMM_WORLD.Send(&j,1,MPI::INT,src,1);
        MPI::COMM_WORLD.Send(swarm[j]->position.data(),numParams,MPI::DOUBLE,src,1);
        ++iter;
      } else {
        ++idle;
        if (vflag) fprintf(stderr,"Sending done signal to process %d.\n",src);
        MPI::COMM_WORLD.Send(0,0,MPI::INT,src,0);
      }

      if (idle == mpi_ntasks-1) break;
    }
  }
}
#endif

// ===========================================================================

void PSO::Swarm::run(int numEvals, int slowdown, int vflag, 
                     ostream* out, ostream* hist) {
  int numIt = numEvals / swarm.size();
  for (int i(0); i < numIt; ++i) {
    if (slowdown) {
      for (int j(0); j < numParams; ++j) {
        if (phiPf[j] > phiPi[j]) 
          phiP[j] = phiPi[j] + (phiPf[j]-phiPi[j])*((1.*numIt-i)/numIt);
        if (phiGf[j] > phiGi[j]) 
          phiG[j] = phiGi[j] + (phiGf[j]-phiGi[j])*((1.*numIt-i)/numIt);
        if (omegaf[j] > omegai[j]) 
          omega[j] = omegai[j] + (omegaf[j]-omegai[j])*((1.*numIt-i)/numIt);
      }
    }
    updateVelocity();
    updatePosition();
    evaluate();
    if (hist != NULL) {
      for (int id(0); id < swarmSize; ++id) {
        *hist << setw(4) << id << " " << scientific << (*swarm[id]) << endl;
      }
    }
    if (out != NULL) {
      *out << setw(6) << i*swarm.size() << " ";
      // *out << scientific << *swarm[bestParticle] << endl;
      *out << scientific << bestVal << " " << bestPos << endl;
    } 
    if (vflag) {
      cerr << "[";
      int progress = i/(numIt/30);
      for (int a(0); a < progress; ++a) cerr << "=";
      for (int a(progress); a < 30; ++a) cerr << " ";
      cerr << "] " << i << "/" << numIt << " best = " << bestPos
           << ", fx = " << scientific << bestVal
           << fixed << "                    \r";
    }
    /*
    if (vflag) {
      cerr << "[" << i << "] best solution: " << bestVal;
      cerr << " at x = (";
      for (size_t j(0); j < bestPos.size(); ++j) {
        cerr << bestPos[j];
        if (j < bestPos.size()-1) cerr << ",";
      }
      cerr << ")" << endl;
    }
    */
  }
  if (vflag) cerr << endl;
}

// ===========================================================================

void PSO::Swarm::display(ostream* out, string prefix) {
  if (mpi_rank == 0) {
    for (int i(0); i < swarmSize; ++i) {
      if (prefix != "") *out << prefix << " ";
      *out << i << " " << (*swarm[i]) << endl;
    }
  }
}

// ===========================================================================

double PSO::Swarm::mean(int j) const {
  double x(0.0);
  for (int i(0); i < swarmSize; ++i) {
    switch (p->type[j]) {
      case INTEGER:
        x += round((*swarm[i])[j]);
        break;
      default:
        x += (*swarm[i])[j];
        break;
    }
  }
  return x/swarmSize;
}

// ===========================================================================

double PSO::Swarm::var(int j) const {
  double m(mean(j));
  double s(0.0);
  double x(0.0);
  for (int i(0); i < swarmSize; ++i) {
    switch (p->type[j]) {
      case INTEGER:
        x = round((*swarm[i])[j])-m;
        break;
      default:
        x = (*swarm[i])[j]-m;
        break;
    }
    s += x*x;
  }
  return s/(swarmSize-1);
}

// ===========================================================================

PSO::Point PSO::Swarm::mean() const {
  Point x(numParams);
  for (int i(0); i < numParams; ++i) x[i] = mean(i);
  return x;
}

// ===========================================================================

PSO::Point PSO::Swarm::var() const {
  Point x(numParams);
  for (int i(0); i < numParams; ++i) x[i] = var(i);
  return x;
}

// ===========================================================================

double PSO::Swarm::maximum_swarm_radius() const {
  double max_rad(0.0);
  double d(0.0);
  for (int i(0); i < swarmSize; ++i) {
    d = bestPos.euclid_dist(swarm[i]->position);
    if (d > max_rad) max_rad = d;
  }
  return d;
}

// ===========================================================================

double PSO::Swarm::search_space_diameter() const {
  Point lb(p->lb);
  Point ub(p->ub);
  return lb.euclid_dist(ub);
}

// ===========================================================================

void PSO::Swarm::init_network() {
  if (numInform < swarmSize-1) {
    double p = 1.-gsl_pow_int(1.-1./swarmSize,numInform);
    for (int i(0); i < swarmSize; ++i) {
      for (int j(0); j < swarmSize; ++j) {
        if (i == j) links.add_link(i,j);
        if (rng.uniform() < p) links.add_link(i,j);
      }
    }
  }
}


