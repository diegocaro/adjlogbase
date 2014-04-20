#ifndef TGRAPHREADER_H_
#define TGRAPHREADER_H_

#include <algorithm>
#include <vector>
#include <sys/types.h>


using namespace std;

struct event {
  uint node;
  uint time;
};



class TGraphReaderEventList {
public:
  static bool myfunction (const struct event &i, const struct event &j) { return (i.time<j.time); }
  
  vector<struct event> events;
  void sort() {
    std::sort(events.begin(), events.end(),myfunction);
  }
  void neighbors(vector<uint> &n) {
    n.clear();
    for(vector<struct event>::iterator it=events.begin(); it!=events.end(); ++it) {
      n.push_back(it->node);
    }
  }
  void timepoints(vector<uint> &n) {
    n.clear();
    for(vector<struct event>::iterator it=events.begin(); it!=events.end(); ++it) {
      n.push_back(it->time);
    }
  }
  
  void purge() {
    vector<struct event>().swap(events);
    events.clear();
  }
  
  uint changes() {
    return events.size();
  }
  
};


class TGraphReaderReverseList {
public:
	vector <uint> neighbors;
};



class TGraphReader {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        TGraphReaderEventList* tgraph;
        TGraphReaderReverseList *revgraph;

        TGraphReader(uint n, uint e, uint c, uint t) {
                nodes = n;
                edges = e;
                changes = c;
                maxtime = t;
                
                
                tgraph = new TGraphReaderEventList[nodes];

                revgraph = new TGraphReaderReverseList[nodes];
        }
        
        void addChange(uint u, uint v, uint t) {
                struct event e;
                e.node = v;
                e.time = t;
                
                tgraph[u].events.push_back(e);
        }

        void addReverseEdge(uint v, uint u) {
        	      revgraph[v].neighbors.push_back(u);
        }
        
        

};



#endif
