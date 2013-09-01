#ifndef TGRAPH_H_
#define TGRAPH_H_ 

#include <fstream>
#include <sys/types.h>
#include "tgraphreader.h"
#include "etdc.h"
#include "coding_policy.h"

#define BLOCKSIZE 128
#define BUFFER 67108864 //256 megabytes

#define CODINGPOLICY "pfor:32:s16:32"

class TGraphEventList {
public:
        uint changes;
        
        uint csize_neighbors;
        uint csize_time;
        
        unsigned char *cneighbors;
        uint *ctime;
};

class TGraph {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        uint *etdctable;
        uint etdcsize;
        
        TGraphEventList* tgraph;
        CodingPolicy *cc;
        
        void loadpolicy() {
                cc = new CodingPolicy(CodingPolicy::kPosition);
                cc->LoadPolicy(CODINGPOLICY);
        }
        
        
        void save(ofstream &out);
        
        static TGraph* load(ifstream &in);
        
        void create(TGraphReader &tgr);
        
        
        void decodetime(uint v, uint *res);
        void decodeneigh(uint v, uint *res);
        
        
        uint snapshot(uint t);
        
        int edge_point(uint u, uint v, uint t);
        int edge_weak(uint u, uint v, uint tstart, uint tend);
        int edge_strong(uint u, uint v, uint tstart, uint tend);
        
        int edge_next(uint u, uint v, uint t);
        
	void direct_point(uint node, uint t, uint *res) ;
	void direct_weak(uint node, uint tstart, uint tend, uint *res) ;
	void direct_strong(uint node, uint tstart, uint tend, uint *res) ;	
	
	/*
        void reverse_point(uint node, uint t, uint *res) const;
	void reverse_weak(uint node, uint tstart, uint tend, uint *res) const;
	void reverse_strong(uint node, uint tstart, uint tend, uint *res) const;
        */
};


void decodediff(uint *k, uint size);
void encodediff(vector<uint> &t);

#endif