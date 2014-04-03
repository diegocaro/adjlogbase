#ifndef TGRAPH_H_
#define TGRAPH_H_ 

#include <fstream>
#include <sys/types.h>
#include <vector>
//#include "tgraphreader.h"
#include "etdc.h"
#include "coding_policy.h"

#define BLOCKSIZE 128
#define BUFFER 67108864 //256 megabytes



#define CP_S9 "s9"
#define CP_S16 "s16"
#define CP_VBYTE "vbyte"
#define CP_RICE "rice"
#define CP_PFOR "pfor:32:s16:32"

using namespace std;

enum CP_FORMAT {
	S9, S16, VBYTE, RICE, PFOR,
};

class TGraphReader;

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
        
	enum CP_FORMAT cp;
	
        TGraphEventList* tgraph;
        CodingPolicy *cc;
        
	void set_policy(enum CP_FORMAT c) {
		cp = c;
	}
	
        void loadpolicy() {
                cc = new CodingPolicy(CodingPolicy::kPosition);
		
		switch(cp) {
			case S9: cc->LoadPolicy(CP_S9); break;
			case S16: cc->LoadPolicy(CP_S16); break;
			case VBYTE: cc->LoadPolicy(CP_VBYTE); break;
			case RICE: cc->LoadPolicy(CP_RICE); break;
			case PFOR: cc->LoadPolicy(CP_PFOR); break;
			default: cc->LoadPolicy(CP_PFOR); break;
		}
		
                
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

        int edge_interval(uint u, uint v, uint tstart, uint tend, uint semantic);
        
        int edge_next(uint u, uint v, uint t);
        
	void direct_point(uint node, uint t, uint *res) ;
	void direct_weak(uint node, uint tstart, uint tend, uint *res) ;
	void direct_strong(uint node, uint tstart, uint tend, uint *res) ;	
	
	void direct_interval(uint node, uint tstart, uint tend, uint semantic, uint *res);
	void reverse_interval(uint node, uint tstart, uint tend, uint semantic, uint *res);
	

        
	/*
        void reverse_point(uint node, uint t, uint *res) const;
	void reverse_weak(uint node, uint tstart, uint tend, uint *res) const;
	void reverse_strong(uint node, uint tstart, uint tend, uint *res) const;
        */
};


void decodediff(uint *k, uint size);
void encodediff(vector<uint> &t);

#endif
