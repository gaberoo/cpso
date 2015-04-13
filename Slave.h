/*-------------------------------------------------------------------------------------*/
/*  
 *  ADAPTED FROM:
 *
 *  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.5        *
 *                                                                                     *
 *  This program is free software: you can redistribute it and/or modify it under the  *
 *  terms of the GNU Lesser General Public License as published by the Free Software   *
 *  Foundation, either version 3 of the License, or (at your option) any later         *
 *  version.                                                                           *
 *                                                                                     *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY    *
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   *
 *                                                                                     *
 *  You should have received a copy of the GNU Lesser General Public License along     *
 *  with this program. If not, see <http://www.gnu.org/licenses/>.                     *
 *
 */

#ifndef __PSO_SLAVE__
#define __PSO_SLAVE__

#include <list>
#include <iostream>

#include "EvalFunc.h"

namespace PSO {
#ifdef USE_MPI
  // MPI constants
  const int   MAX_REQ_WAIT =  3 ;  // Maximum time to wait for a request
  const char   STOP_SIGNAL = 'S';  // Stop signal
  const char   EVAL_SIGNAL = 'X';  // Evaluation signal
  const char  READY_SIGNAL = 'R';  // Ready signal
  const char RESULT_SIGNAL = 'O';  // Result signal
  const char   WAIT_SIGNAL = 'W';  // Wait signal
#endif

  // Slave process for the parallel version.
  class Slave {
  private:
    static int  _rank;        // Process rank.
    static int  _np;          // Number of processes.
    static int  _data_sent;   // Stats on the sent data.
    static int  _data_rcvd;   // Stats on the received data.
    static bool _are_running; //< \c true if the slaves are running.
    static bool _stop_ok;     //< \c true if the slaves stopped without error.
    PSO::EvalFunc* _ev;       // Evaluator  (may be NULL).

    /// Initializations.
    void init ( void ) const;

    /// Force quit.
    /**
       - Called by pressing Ctrl-C.
       - Does nothing: the slave will be stopped by the master.
       \param signalValue Signal value -- \b IN.
    */
    static void force_quit ( int signalValue ) {}

  public:
    Slave(PSO::EvalFunc* ev) : _ev (ev) { init(); } 
    virtual ~Slave() {}

    // Run the slave code.
    void run() const;

    inline static int get_data_sent() { return Slave::_data_sent; }
    inline static int get_data_rcvd() { return Slave::_data_rcvd; }
    static int get_rank();
    static int get_nb_processes();
    inline static bool are_running() { return _are_running; }
    inline static bool is_master() { return (Slave::get_rank() == 0); }
    static void init_slaves();
    static void stop_slaves();

#ifdef USE_MPI
  private:
    Point* eval() const;

    void send_eval_result(const Point* x) const;
    
    /// Send data.
    /**
       \param buf        Data to send  -- \b IN.
       \param count      Data quantity -- \b IN.
       \param datatype   Data type     -- \b IN.
       \param dest       Destination   -- \b IN.
       \param ready_send Flag equal to \c true if \c Rsend is used instead of \c Send
                                       -- \b IN.
    */
    static void send_data(const void* buf, int count, 
        MPI_Datatype datatype, int dest, bool ready_send);
    
    /// Receive data.
    /**
       \param buf      Data to receive                           -- \b OUT.
       \param count    Data quantity                             -- \b IN.
       \param datatype Data type                                 -- \b IN.
       \param source   Source (may be \c MPI_ANY_SOURCE)         -- \b IN.
       \param req      Pointer to a MPI request (may be \c NULL) -- \b IN/OUT.
       \return         The source.
    */
    static int receive_data(void* buf, int count, 
        MPI_Datatype datatype, int source, MPI_Request* req);

    /// Wait for a MPI request.
    /**
       \param req The request -- \b IN/OUT.
    */
    static void wait_request(MPI_Request& req);
    
  public:

    void send_eval_point(const Point* x, int slave_rank) const;
    void receive_eval_result(int slave_rank, Point* x, bool& eval_ok) const;
    
    static int receive_signal(char & signal) {
      return Slave::receive_data(&signal,1,MPI_CHAR,MPI_ANY_SOURCE,NULL);
    }

    static void send_signal (char signal, int slave_rank) {
      Slave::send_data(&signal,1,MPI_CHAR,slave_rank,true);
    }
#endif // USE_MPI

  };
}
#endif
