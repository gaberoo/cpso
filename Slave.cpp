#include "Slave.h"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
int  PSO::Slave::_rank        = -1;
int  PSO::Slave::_np          = -1;
int  PSO::Slave::_data_sent   =  0;
int  PSO::Slave::_data_rcvd   =  0;
bool PSO::Slave::_are_running = false;
bool PSO::Slave::_stop_ok     = false;

/*----------------------------------------*/
/*        initializations (private)       */
/*----------------------------------------*/
void NOMAD::Slave::init() const {
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&PSO::Slave::_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&PSO::Slave::_np);
#else
  PSO::Slave::_rank = 0;
  PSO::Slave::_np   = 1;
#endif

  // Slave::force_quit() will be called if ctrl-c is pressed:
  if (! PSO::Slave::is_master()) {
    PSO::EvalFunc::force_quit();
    signal(SIGTERM,PSO::Slave::force_quit);
    signal(SIGINT,PSO::Slave::force_quit);
    signal(SIGPIPE,PSO::Slave::force_quit);
  }
}

/*----------------------------------------*/
/*          get the process rank          */
/*               (static)                 */
/*----------------------------------------*/
int PSO::Slave::get_rank() {
  if (PSO::Slave::_rank < 0) {
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&PSO::Slave::_rank);
#else
    PSO::Slave::_rank = 0;
#endif
  }
  return PSO::Slave::_rank;
}

/*----------------------------------------*/
/*       get the number of processes      */
/*               (static)                 */
/*----------------------------------------*/
int PSO::Slave::get_nb_processes() {
  if (PSO::Slave::_np < 0) {
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&PSO::Slave::_np);
#else
    PSO::Slave::_np = 1;
#endif
  }
  return PSO::Slave::_np;
}

/*----------------------*/
/*  run the slave code  */
/*----------------------*/
void NOMAD::Slave::run() const {
#ifdef USE_MPI
  MPI_Request req;
  char signal(0);
  PSO::Point* x = NULL;

  while (true) {
    PSO::Slave::receive_data(&signal,1,MPI_CHAR,0,&req);
    PSO::Slave::send_data(&PSO::READY_SIGNAL,1,MPI_CHAR,0,false);
    PSO::Slave::wait_request(req);

    if (signal == PSO::EVAL_SIGNAL) x = eval();
    else if (signal == PSO::RESULT_SIGNAL) {
      send_eval_result(x);
      delete x;
      x = NULL;
    } 
    else if (signal == PSO::STOP_SIGNAL) break;
    // else if ( signal == NOMAD::WAIT_SIGNAL ) {}
  }

  if (x) delete x;
#endif
}

/*-----------------------------------------*/
/*        initialize all the slaves        */
/*                (static)                 */
/*-----------------------------------------*/
void PSO::Slave::init_slaves() {
#ifdef USE_MPI
  // don't reinitialize
  if (! PSO::Slave::is_master() || PSO::Slave::_are_running) return;

  // out << std::endl << NOMAD::open_block ( "initializing slaves" );

  MPI_Status status;
  MPI_Request** req = new MPI_Request*[PSO::Slave::_np];
  int nb_initialized = 0;
  int nb_slaves = PSO::Slave::_np-1;
  int source;
  char signal;
  
  // 1. launch requests:
  for (source = 1; source < PSO::Slave::_np; ++source) {
    req[source] = new MPI_Request;
    PSO::Slave::receive_data(&signal,1,MPI_CHAR,source,req[source]);
    // out << "." << std::endl;
  }

  // 2. test requests (with a maximal delay of MAX_REQ_WAIT):
  int cnt = 0;
  int flag;
  while (nb_initialized < nb_slaves) {
    for (source = 1; source < NOMAD::Slave::_np; ++source) {
      if (req[source]) {
        MPI_Test(req[source],&flag,&status);
	      if (flag) {
	        MPI_Wait(req[source],&status);
          PSO::Slave::send_data(&PSO::WAIT_SIGNAL,1,MPI_CHAR,source,true);
	        delete req[source];
	        req[source] = NULL;
	        ++nb_initialized;
	      }
      }
    }
    // a constant is used in order to display only a few '.' :
//    if ( display_degree == NOMAD::FULL_DISPLAY && cnt%1000000==0 )
//      out << "." << std::endl;
//
    ++cnt;
  }
    
  // 3. delete requests:
  std::list<int> err_list;
  for (source = 1; source < PSO::Slave::_np; ++source) {
    if (req[source]) {
      err_list.push_back(source);
      MPI_Cancel(req[source]);
      delete req[source];
    }
  }
  delete[] req;

  PSO::Slave::_are_running = true;
  PSO::Slave::_stop_ok     = false;

  if (! err_list.empty()) {
    std::ostringstream oss;
    oss << "could not initialize slave";
    if ( err_list.size() > 1 ) {
      oss << "s";
      std::list<int>::const_iterator it, end = err_list.end();
      for (it = err_list.begin(); it != end; ++it) oss << " #" << *it;
    }
    else oss << " #" << *err_list.begin();
    cerr << oss.str() << endl;
    // throw NOMAD::Exception ("Slave.cpp" , __LINE__ , oss.str());
  }
#endif
}

/*-----------------------------------------*/
/*             stop the slaves             */
/*                (static)                 */
/*-----------------------------------------*/
void PSO::Slave::stop_slaves() {
#ifdef USE_MPI
  if (! PSO::Slave::is_master() || PSO::Slave::_stop_ok) return;

  // std::cerr << std::endl << "stopping slaves" << std::endl;

  int nb_stopped = 0;
  int nb_slaves = PSO::Slave::_np-1;
  int source;
  char signal;

  // NOMAD::Clock clk;

  MPI_Status status;
  MPI_Request** req = new MPI_Request*[PSO::Slave::_np];

  // 1. launch requests:
  for (source = 1; source < PSO::Slave::_np; ++source) {
    req[source] = new MPI_Request;
    PSO::Slave::receive_data(&signal,1,MPI_CHAR,source,req[source]);
  }

  // 2. test requests (with a maximal delay of MAX_REQ_WAIT):
  int cnt = 0;
  int flag;
  while (nb_stopped < nb_slaves) {
    for (source = 1; source < PSO::Slave::_np; ++source) {
      if (req[source]) {
        MPI_Test(req[source],&flag,&status);
        if (flag) { 
          MPI_Wait(req[source],&status);
          PSO::Slave::send_data(&PSO::STOP_SIGNAL,1,MPI_CHAR,source,true);
          delete req[source];
          req[source] = NULL;
          ++nb_stopped;
        }
      }
    }
    ++cnt;
  }

  PSO::Slave::_are_running = false;
  PSO::Slave::_stop_ok     = true;

  // 3. delete requests:  
  for (source = 1; source < PSO::Slave::_np ; ++source) {
    if (req[source]) {
      MPI_Cancel(req[source]);
      delete req[source];
      PSO::Slave::_stop_ok = false;
    }
  }
  delete[] req;
#endif

}

#ifdef USE_MPI

/*------------------------------------------------------*/
/*               receive data (static, private)         */
/*------------------------------------------------------*/
int PSO::Slave::receive_data (void* buf, int count, 
    MPI_Datatype datatype, int source, MPI_Request* req)
{
  int tag = (PSO::Slave::is_master()) ? source : PSO::Slave::get_rank();

  // immediate receive:
  if (req) {
    if (source == MPI_ANY_SOURCE)
	    cerr << "Slave::receive_data(): immediate receive with no source" << endl;
    MPI_Irecv(buf,count,datatype,source,tag,MPI_COMM_WORLD,req);
  }

  // normal receive:
  else {
    MPI_Status status;
    if (source == MPI_ANY_SOURCE) tag = MPI_ANY_TAG;
    MPI_Recv(buf,count,datatype,source,tag,MPI_COMM_WORLD,&status);
    source = status.MPI_SOURCE;
  }

  // stats:
  int size;
  MPI_Type_size(datatype,&size);
  PSO::Slave::_data_rcvd += count*size;
  
  return source;
}

/*------------------------------------------------------*/
/*              send data (static, private)             */
/*------------------------------------------------------*/
void NOMAD::Slave::send_data(const void* buf, int count, 
    MPI_Datatype datatype, int dest, bool ready_send)
{
  int tag = (PSO::Slave::is_master()) ? dest : PSO::Slave::get_rank();

  // ready send:
  if (ready_send)
    MPI_Rsend(const_cast<void*>(buf),count,datatype,dest,tag,MPI_COMM_WORLD);

  // normal send:
  else
    MPI_Send(const_cast<void*>(buf),count,datatype,dest,tag,MPI_COMM_WORLD);

  // stats:
  int size;
  MPI_Type_size(datatype,&size);
  PSO::Slave::_data_sent += count*size;
}

/*------------------------------------------------------*/
/*  receive and evaluate an Eval_Point from the master  */
/*  (private)                                           */
/*------------------------------------------------------*/
PSO::Point* PSO::Slave::eval() const {
  // 1. receive the point:
  int itab[3];
  MPI_Request req;
  PSO::Slave::receive_data(itab,3,MPI_INT,0,&req);
  PSO::Slave::send_data(&PSO::READY_SIGNAL,1,MPI_CHAR,0,false);
  PSO::Slave::wait_request(req);

  int n = itab[0];
  double* dtab = new double[n+1];
  PSO::Slave::receive_data(dtab,n+1,MPI_DOUBLE,0,&req);
  PSO::Slave::send_data(&PSO::READY_SIGNAL,1,MPI_CHAR,0,false);
  PSO::Slave::wait_request(req);

  // 2. create the Eval_Point:
  PSO::Point* x = new PSO::Point(n,_p->get_bb_nb_outputs());
  for ( int i = 0 ; i < n ; ++i )
    (*x)[i] = dtab[i];
  NOMAD::Double h_max = dtab[n];

  x->set_tag       ( itab[2]                            );
  x->set_eval_type ( ( itab[1] > 0 ) ? NOMAD::SGTE : NOMAD::TRUTH );

  delete [] dtab;
  
  // 3. evaluate the point:
  bool eval_ok;
  try {
    eval_ok = _ev->eval_x ( *x , h_max , count_eval );
  }
  catch ( ... ) {
    eval_ok = false;
  }

  x->set_eval_status ( ( eval_ok ) ? NOMAD::EVAL_OK : NOMAD::EVAL_FAIL );

  return x;
}

/*-----------------------------------------------------*/
/*  send an evaluation result to the master (private)  */
/*-----------------------------------------------------*/
void NOMAD::Slave::send_eval_result ( const NOMAD::Eval_Point * x          ,
				      bool                      count_eval   ) const
{
  // receive a signal from the master:
  char signal;
  NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , 0 , NULL );

  // send the evaluation result:
  int                  m    = _p->get_bb_nb_outputs();
  int                  s    = 2*m+2;
  double             * dtab = new double [s];
  const NOMAD::Point & bbo  = x->get_bb_outputs();

  // bb_outputs (m values):
  for ( int i = 0 ; i < m ; ++i ) {
    if ( bbo[i].is_defined() ) {
      dtab[i  ] = bbo[i].value();
      dtab[i+m] = 1.0;
    }
    else {
      dtab[i  ] = NOMAD::INF;
      dtab[i+m] = -1.0;
    }
  }

  // evaluation status:
  dtab[2*m] = ( x->get_eval_status() == NOMAD::EVAL_OK ) ? 1.0 : -1.0;

  // count_eval:
  dtab[s-1] = ( count_eval ) ? 1.0 : -1.0;

  // send the array:
  NOMAD::Slave::send_data ( dtab , s , MPI_DOUBLE , 0 , true );

  delete [] dtab;
}

/*---------------------------------------------*/
/*  receive an evaluation result from a slave  */
/*---------------------------------------------*/
void NOMAD::Slave::receive_eval_result ( int                 slave_rank ,
					 NOMAD::Eval_Point * x          ,
					 bool              & eval_ok    ,
					 bool              & count_eval    ) const
{
  // send the RESULT signal to the slave:
  NOMAD::Slave::send_data ( &NOMAD::RESULT_SIGNAL , 1 , MPI_CHAR , slave_rank , true );

  // receive the evaluation result as a double array:
  int      m    =  _p->get_bb_nb_outputs();
  int      s    = 2*m+2;
  double * dtab = new double [s];

  MPI_Request req;
  NOMAD::Slave::receive_data ( dtab              , s , MPI_DOUBLE , slave_rank , &req  );
  NOMAD::Slave::send_data ( &NOMAD::READY_SIGNAL , 1 , MPI_CHAR   , slave_rank , false );
  NOMAD::Slave::wait_request ( req );

  // interpret the array:
  for ( int i = 0 ; i < m ; ++i )
    x->set_bb_output ( i , ( dtab[i+m] > 0.0 ) ? dtab[i] : NOMAD::Double() );
  
  eval_ok = ( dtab[2*m] > 0.0 );

  x->set_eval_status ( eval_ok ? NOMAD::EVAL_OK : NOMAD::EVAL_FAIL );

  count_eval = ( dtab[s-1] > 0.0 );

  delete [] dtab;
}

/*-----------------------------------------------------*/
/*            send an Eval_Point to a slave            */
/*-----------------------------------------------------*/
void NOMAD::Slave::send_eval_point ( const NOMAD::Eval_Point * x          ,
				     int                       slave_rank ,
				     const NOMAD::Double     & h_max        ) const
{
  char signal;
  int  itab[3];
  int  n = x->size();
    
  // n:
  itab[0] = n;

  // evaluation type (+1: sgte eval; -1: true eval):
  itab[1] = ( x->get_eval_type() == NOMAD::SGTE ) ? 1 : -1;

  // tag of the point:
  itab[2] = x->get_tag();

  // point coordinates:
  double * dtab = new double[n+1];
  for ( int i = 0 ; i < n ; ++i )
    dtab[i] = (*x)[i].value();
  dtab[n] = h_max.value();

  // wait for the slave signal:
  NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , slave_rank , NULL );
  
  // send n and evaluation type:
  NOMAD::Slave::send_data ( itab , 3 , MPI_INT , slave_rank , true );

  // wait for the slave signal:
  NOMAD::Slave::receive_data ( &signal , 1 , MPI_CHAR , slave_rank , NULL );

  // send the point coordinates:
  NOMAD::Slave::send_data ( dtab , n+1 , MPI_DOUBLE , slave_rank , true );

  delete [] dtab;
}

/*-----------------------------------------*/
/*          wait for a MPI request         */
/*                (static)                 */
/*-----------------------------------------*/
void PSO::Slave::wait_request(MPI_Request& req)
{
  MPI_Status status;
  MPI_Wait(&req,&status);
}

#endif
