#ifndef TGRAPHREADER_H_
#define TGRAPHREADER_H_

#include <vector>
#include <sys/types.h>


using namespace std;



class TGraphReaderEventList {
public:
        uint changes;
        vector<uint> neighbors;
        vector<uint> timepoints;
        
        void addEvent(uint v, uint t) {
                neighbors.push_back(v);
                timepoints.push_back(t);
                changes++;
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
                tgraph[u].addEvent(v,t);
        }

        void addReverseEdge(uint v, uint u) {
        	revgraph[v].neighbors.push_back(u);
        }

};



#endif
