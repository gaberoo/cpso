#ifndef __PSO_NETWORK__
#define __PSO_NETWORK__

#include <vector>
#include <set>
#include <algorithm>
using namespace std;

namespace PSO {
  class Network {
    public:
      Network() {}
      Network(int size) : order(size), links(size) {}
      virtual ~Network() {}

      inline void resize(int new_size) {
        clear();
        links.resize(new_size);
      }

      inline void clear() { links.clear(); }

      inline set<int>::const_iterator get_link(int i, int j) const {
        return links[i].find(j);
      }

      inline void add_link(int i, int j) {
        if (get_link(i,j) != nn_end(i)) links[i].insert(j);
      }

      inline set<int>::const_iterator nn_begin(int i) const {
        return links[i].begin();
      }

      inline set<int>::const_iterator nn_end(int i) const {
        return links[i].end();
      }

    protected:
      int order;
      vector< set<int> > links;
  };
}

#endif
