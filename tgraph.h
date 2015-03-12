#ifndef TGRAPH_H_
#define TGRAPH_H_ 

#include <fstream>
#include <sys/types.h>
#include <vector>
//#include "tgraphreader.h"
#include "etdc.h"
#include "coding_policy.h"

#define BLOCKSIZE 128
#define BUFFER 67108864*4 //256*4 megabytes



//#define CP_S9 "s9"
//#define CP_S16 "s16"
//#define CP_VBYTE "vbyte"
//#define CP_RICE "rice"
//#define CP_PFOR "pfor:32:rice:32"
//#define CP_PFOR32 "pfor:32:s16:32"
//#define CP_PFOR64 "pfor:64:s16:32"
//#define CP_PFOR128 "pfor:128:s16:32"  //codec:blocksize:coder:paddingsize

// Defines names for each coding method.
//static const char kRiceCoding[] = "rice";
//static const char kTurboRiceCoding[] = "turbo-rice";
//static const char kPForDeltaCoding[] = "pfor";
//static const char kS9Coding[] = "s9";
//static const char kS16Coding[] = "s16";
//static const char kVarByteCoding[] = "vbyte";
//static const char kNullCoding[] = "null";


using namespace std;

//enum CP_FORMAT {
//	S9, S16, VBYTE, RICE, PFOR32, PFOR64, PFOR128, PFOR
//};


enum TypeGraph {
	kInterval, kGrowth, kPoint
};

struct opts {
	char *outfile;
	char cp[100]; //coding policy
	enum TypeGraph typegraph;
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


class TGraphReverse {
public:
	uint size; //out degree, number of elements
	uint csize; //compressed size
	uint *clist; //pointer to compressed list
};



class TGraph {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        struct opts opts;

        uint *etdctable;
        uint etdcsize;
	
        TGraphEventList* tgraph;

        TGraphReverse *reverse;

        CodingPolicy *cc;
        
		
        void set_opts(struct opts o) {
        	opts = o;
        }
	
        void loadpolicy() {
                cc = new CodingPolicy(CodingPolicy::kPosition);
				cc->LoadPolicy(opts.cp);
        }
        
        
        void save(ofstream &out);
        
        static TGraph* load(ifstream &in);
        
        void create(TGraphReader &tgr);
        
        
        void decodetime(uint v, uint *res);
        void decodeneigh(uint v, uint *res);
        void decodereverse(uint v, uint *res);
        
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
	

        

        void reverse_point(uint node, uint t, uint *res) ;
	void reverse_weak(uint node, uint tstart, uint tend, uint *res) ;
	void reverse_strong(uint node, uint tstart, uint tend, uint *res) ;


    size_t change_point(uint t);
    size_t change_interval(uint ts, uint te);
    size_t change_node(uint u, uint ts, uint te);

    size_t actived_point(uint t);
    size_t actived_interval(uint ts, uint te);
    size_t actived_node(uint u, uint ts, uint te);


    size_t deactived_point(uint t);
    size_t deactived_interval(uint ts, uint te);
    size_t deactived_node(uint u, uint ts, uint te);
	

    // point contact graphs
    void direct_interval_pg(uint node, uint tstart, uint tend, uint *res) ;
    //void reverse_interval_pg(uint node, uint tstart, uint tend, uint *res) ; //captured by edge_interval_pg
    int edge_interval_pg(uint u, uint v, uint tstart, uint tend) ;
    int edge_next_pg(uint u, uint v, uint t);
    //size_t snapshot_pg(uint t) ; //captured by direct_point

    //size_t change_interval_pg(uint ts, uint te); //do not need an upgrade
    size_t actived_interval_pg(uint ts, uint te);
    size_t deactived_interval_pg(uint ts, uint te);


//	size_t getSize() {
//		size_t size_time=0;
//		size_t size_events=0;
//
//		for(uint i=0; i < nodes; i++) {
//			size_time += (tgraph[i]->csize_time )*sizeof(uint);
//			size_events += tgraph[i]->csize_neighbors;
//		}
//
//		printf("neighbors: %.2lf MBytes\n", size_events/1024/1024);
//		printf("time: %.2lf MBytes\n", size_time/1024/1024):
//
//		return size_time + size_events;
//	}
};


void decodediff(uint *k, uint size);
void encodediff(vector<uint> &t);

#endif
