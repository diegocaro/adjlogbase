#ifndef TGRAPHREADER_H_
#define TGRAPHREADER_H_

#include "btree_map.h"
#include "btree_set.h"

#include <algorithm>
#include <vector>
#include <sys/types.h>

using namespace btree;
using namespace std;

class Event {
public:
  uint node;
  uint time;
  
  Event() {}
  
  Event( const Event &obj) {
	  node = obj.node;
	  time = obj.time;
  }
  
	inline bool operator<(const Event& e) const {
		if (this->time == e.time) return (node < e.node);
		return  (this->time < e.time);
	}

	uint getnode() { return node;}
	uint gettime() { return this->time;}
};


class TGraphReader {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        //TGraphReaderEventList* tgraph;
		btree_map< uint, btree_set <Event> > tgraph;
		
        //TGraphReaderReverseList *revgraph;
		btree_map< uint, btree_set <uint> > revgraph;

        TGraphReader(uint n, uint e, uint c, uint t) {
                nodes = n;
                edges = e;
                changes = c;
                maxtime = t;
                
                
                //tgraph = new TGraphReaderEventList[nodes];

                //revgraph = new TGraphReaderReverseList[nodes];
        }
        
        void addChange(uint u, uint v, uint t) {
                Event e;
                e.node = v;
                e.time = t;
                
                tgraph[u].insert(e);
        }

        void addReverseEdge(uint v, uint u) {
        	revgraph[v].insert(u);
        }
        
        

};



#endif
