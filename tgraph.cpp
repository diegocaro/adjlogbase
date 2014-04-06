#include "tgraphreader.h"
#include "tgraph.h"


#include "uthash.h"
#include "etdc.h"

#include "cppUtils.h"
#include "debug.h"

#include "arraysort.h"

#include <unordered_map>
#include <map>
using namespace std;


void encodediff(vector<uint> &t) {
        uint val, old;
        old = t[0];
        for(uint j=1; j < t.size(); j++) {                        
                val = t[j];
                t[j] -= old;
                old = val;
        }
}

void decodediff(uint *k, uint size) {
        if (size == 0) return;
        
        for(uint i=1; i < size; i++) {
                k[i] += k[i-1];
        }
}


void TGraph::save(ofstream &out) {
        cds_utils::saveValue<TGraph>(out, this, 1);
        
        cds_utils::saveValue<TGraphEventList>(out, tgraph, nodes);
        
        for(uint i=0; i < nodes; i++) {
                if (tgraph[i].changes > 0) {
                        cds_utils::saveValue<unsigned char>(out, tgraph[i].cneighbors, tgraph[i].csize_neighbors);
                        cds_utils::saveValue<uint>(out, tgraph[i].ctime, tgraph[i].csize_time);
                }
        }
        
        cds_utils::saveValue<uint>(out, etdctable, etdcsize);


        //save reverse graph
        cds_utils::saveValue<TGraphReverse>(out, reverse, nodes);
        for(uint i=0; i < nodes; i++) {
        	if (reverse[i].size > 0) {
        		cds_utils::saveValue<uint>(out, reverse[i].clist, reverse[i].csize);
        	}
        }
}

TGraph* TGraph::load(ifstream &in) {
        TGraph *tg;
        
        tg = cds_utils::loadValue<TGraph>(in, 1);
        
        LOG("nodes: %u", tg->nodes);
        LOG("edges: %u", tg->edges);
        LOG("changes: %u", tg->changes);
        LOG("maxtime: %u", tg->maxtime);
        LOG("etdcsize: %u", tg->etdcsize);
        
        tg->tgraph = cds_utils::loadValue<TGraphEventList>(in, tg->nodes);
        
        for(uint i=0; i < tg->nodes; i++) {
                if (tg->tgraph[i].changes > 0) {
                        tg->tgraph[i].cneighbors = cds_utils::loadValue<unsigned char>(in, tg->tgraph[i].csize_neighbors);
                        tg->tgraph[i].ctime = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_time);
                }
        }
        
        tg->etdctable = cds_utils::loadValue<uint>(in, tg->etdcsize);
        
        //read reverse graph
        tg->reverse = cds_utils::loadValue<TGraphReverse>(in, tg->nodes);
        for(uint i=0; i < tg->nodes; i++) {
                if (tg->reverse[i].size > 0) {
                	tg->reverse[i].clist =  cds_utils::loadValue<uint>(in, tg->reverse[i].csize);
                }
        }

	tg->set_policy(tg->cp);
        tg->loadpolicy();
        
        return tg;
}






void TGraph::create(TGraphReader &tgr) {
        //uint *nodesbuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        //uint *timebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        unsigned char *ucharbuffer = new unsigned char[ (BUFFER/BLOCKSIZE+1) * BLOCKSIZE ];
        uint *uintbuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        //uint *timebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        
        nodes = tgr.nodes;
        edges = tgr.edges;
        changes = tgr.changes;
        maxtime = tgr.maxtime;
        
        tgraph = new TGraphEventList[nodes];
        
        // First pass for ETDC
        struct etdc_table *table = NULL;
        for(uint i=0; i < tgr.nodes; i++) {
                for (uint j = 0; j < tgr.tgraph[i].changes; j++) {
                         etdc_add(&table, tgr.tgraph[i].neighbors[j]);
                }
        }
        etdc_sort(&table);
        etdc_gencodes(table);
        
        etdcsize = etdc_size(table);
        etdctable = new uint[etdcsize];
        etdc_voc2uint(table, etdctable);
        
        
        LOG("nodes: %u", nodes);
        LOG("edges: %u", edges);
        LOG("changes: %u", changes);
        LOG("maxtime: %u", maxtime);
        LOG("etdcsize: %u", etdcsize);
        
        this->loadpolicy(); //load compression policy for time
        
        // Second pass, to really compress :)
        uint csize_time, csize_neigh, node_changes;
        for(uint i=0; i < nodes; i++) {
                if (i%1000==0) fprintf(stderr, "Compressing: %0.2f%%\r", (float)i*100/nodes);
                node_changes = tgr.tgraph[i].changes;
                
                tgraph[i].changes = 0;
                tgraph[i].csize_time = 0;
                tgraph[i].csize_neighbors = 0;
                
                if (node_changes == 0) { LOG("node %u with zero changes", i) ;continue;}
                 
                csize_neigh = etdc_encode(&table, tgr.tgraph[i].neighbors.data(), node_changes, ucharbuffer);
                tgraph[i].cneighbors = new unsigned char [csize_neigh];
                memcpy(tgraph[i].cneighbors, ucharbuffer, csize_neigh);
                
                encodediff(tgr.tgraph[i].timepoints);
                csize_time = cc->Compress(tgr.tgraph[i].timepoints.data(), uintbuffer, node_changes);
                tgraph[i].ctime = new uint [csize_time];
                memcpy(tgraph[i].ctime, uintbuffer, csize_time * sizeof(uint));
                //printf("Compression time ratio: %f\n", (float)csize_time/node_changes);


                tgraph[i].changes = node_changes;
                tgraph[i].csize_time = csize_time;
                tgraph[i].csize_neighbors = csize_neigh;
		
		

                //force free memory
                vector<uint>().swap(tgr.tgraph[i].neighbors);
                vector<uint>().swap(tgr.tgraph[i].timepoints);
    
		tgr.tgraph[i].neighbors.clear();
		tgr.tgraph[i].timepoints.clear();
        }
        fprintf(stderr, "\n");
        
        // Creating reverse structure
        reverse = new TGraphReverse[nodes];
        uint csize,size;
        for(uint i=0; i < nodes; i++) {
                if (i%1000==0) fprintf(stderr, "Compressing reverse graph: %0.2f%%\r", (float)i*100/nodes);

                reverse[i].size = 0;
                reverse[i].csize = 0;
                reverse[i].clist = NULL;

                size = tgr.revgraph[i].neighbors.size();

                if ( size == 0 ) {
                	continue;
                }

                encodediff(tgr.revgraph[i].neighbors);
                csize = cc->Compress(tgr.revgraph[i].neighbors.data(), uintbuffer, size);

                reverse[i].size = size;
                reverse[i].csize = csize;
                reverse[i].clist = new uint[csize];
                memcpy(reverse[i].clist, uintbuffer, csize * sizeof(uint));


                //force free memory
                vector<uint>().swap(tgr.revgraph[i].neighbors);
                tgr.revgraph[i].neighbors.clear();

        }
        
        delete [] uintbuffer;
        delete [] ucharbuffer;

}




void TGraph::decodetime(uint v, uint *res) {
        //if (tgraph[v].changes == 0) return;
         
        cc->Decompress(tgraph[v].ctime, res, tgraph[v].changes);
        decodediff(res, tgraph[v].changes);
}

void TGraph::decodeneigh(uint v, uint *res) {
        //if (tgraph[v].changes == 0) return;
        etdc_decode(etdctable, etdcsize, tgraph[v].cneighbors, 
                        tgraph[v].csize_neighbors, res, tgraph[v].changes);
        
}

void TGraph::decodereverse(uint v, uint *res) {
        //if (tgraph[v].changes == 0) return;

        cc->Decompress(reverse[v].clist, res, reverse[v].size);
        decodediff(res, reverse[v].size);
}


void TGraph::direct_point(uint v, uint t, uint *res)  {
        if (v>=nodes || tgraph[v].changes == 0) return;
        uint *timep = new uint[BLOCKSIZE*tgraph[v].changes];
        uint *nodep = new uint[BLOCKSIZE*tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        
        uint i=0;
        
        
        unordered_map <uint,uint> thash;
        for(uint j=0; j < tgraph[v].changes; j++) {
          if (timep[j] > t) break;
          thash[nodep[j]]++;
        }
        for ( unordered_map<uint,uint>::iterator it = thash.begin(); it != thash.end(); ++it ) {
          if (it->second %2 == 1) {
//            printf("first %u second %u\n",it->first,it->second);
            res[++i] = it->first;
          }
          
        }
        *res = i;
        
        
        /*
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] > t) break;
                res[++i] = nodep[j];
        }
        
        *res = i;
        qsort(&res[1], *res, sizeof(unsigned int), compare);
        xor_arraysort(res);
        */
        delete [] timep;
        delete [] nodep;
}

void TGraph::direct_interval(uint v, uint tstart, uint tend, uint semantic, uint *res)  {
        if (v>=nodes|| tgraph[v].changes == 0) return;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].changes];
        uint *nodep = new uint[BLOCKSIZE*tgraph[v].changes];
        
        uint *buffer = new uint[BLOCKSIZE*tgraph[v].changes];
        uint *interval = new uint[BLOCKSIZE*tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        
        
        
        unordered_map <uint,uint> thashp; //point 
        unordered_map <uint,uint> thashi; //interval
        for(uint j=0; j < tgraph[v].changes; j++) {
          if (timep[j] <= tstart) {
              thashp[nodep[j]]++;
          }
          
          else if (timep[j] > tstart && timep[j]<= tend) {
              thashi[nodep[j]] = 1;
          }
          
          if (timep[j] > tend) break;
        }
        
        
        uint i=0;
        for ( unordered_map<uint,uint>::iterator it = thashp.begin(); it != thashp.end(); ++it ) {
          if (it->second %2 == 1) {
            buffer[++i] = it->first;
          }
          
        }
        *buffer = i;
        
        uint k=0;
        for ( unordered_map<uint,uint>::iterator it = thashi.begin(); it != thashi.end(); ++it ) {
            interval[++k] = it->first;        
        }
        *interval = k;
        
        qsort(&buffer[1], *buffer, sizeof(unsigned int), compare);    
        qsort(&interval[1], *interval, sizeof(unsigned int), compare);
        
        //this semantic filter is O(d) where d is the out degree of the node
        if (semantic == 0) {
                merge_arraysort(res, buffer, interval);
        }
        else if (semantic == 1) {
                diff_arraysort(buffer, interval);
                memcpy(res, buffer, (*buffer+1)*sizeof(uint));
        }
        
        /*
        uint i=0;
        uint k=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] <= tstart) {
                        buffer[++i] = nodep[j];
                }
                
                else if (timep[j] > tstart && timep[j]<= tend) {
                        interval[++k] = nodep[j];
                }
                
                if (timep[j] > tend) break;
        }
        *buffer = i;
        *interval = k;
        
        qsort(&buffer[1], *buffer, sizeof(unsigned int), compare);
        xor_arraysort(buffer);
        
        qsort(&interval[1], *interval, sizeof(unsigned int), compare);
        remove_duplicates(interval);
        
        if (semantic == 0) {
                merge_arraysort(res, buffer, interval);
        }
        else if (semantic == 1) {
                diff_arraysort(buffer, interval);
                memcpy(res, buffer, (*buffer+1)*sizeof(uint));
        }
        */
        
        delete [] timep;
        delete [] nodep;
        delete [] buffer;
        delete [] interval;
}


void TGraph::direct_weak(uint v, uint tstart, uint tend, uint *res)  {
        direct_interval(v, tstart, tend, 0, res);
}


void TGraph::direct_strong(uint v, uint tstart, uint tend, uint *res)  {
        direct_interval(v, tstart, tend, 1, res);
}


uint TGraph::snapshot(uint t){
        uint *buffer = new uint [BUFFER];
        
        uint edges=0;
        for(uint v=0; v < nodes; v++) {
                direct_point(v, t, buffer);
                edges += *buffer;
        }
        
        delete [] buffer;
        
        return edges;
}

int TGraph::edge_point(uint v, uint u, uint t){
        if (v>=nodes || tgraph[v].changes == 0) return 0;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].changes];
        uint *nodep = new uint[BLOCKSIZE*tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] > t) break;
                if( u == nodep[j]) occ++;
        }
        
        delete [] timep;
        delete [] nodep;

        return (occ%2);
}

int TGraph::edge_interval(uint v, uint u, uint tstart, uint tend, uint semantic){
        if (v>=nodes || tgraph[v].changes == 0) return 0;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].changes];
        uint *nodep = new uint[BLOCKSIZE*tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;
        uint occinterval=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] <= tstart) {
                        if( u == nodep[j]) occ++;
                }
                else if (timep[j] > tstart && timep[j]<= tend) {
                        if( u == nodep[j]) occinterval++;
                }
                
                if (timep[j] > tend) break;
                
        }
        
        delete [] timep;
        delete [] nodep;
        
        if (semantic == 0) {
                if ( (occ%2) || occinterval > 0)
                        return 1;
        }
        else if (semantic == 1) {
                if ( (occ%2) && occinterval == 0 ) {
                        return 1;
                }
        }
        return 0;
        
        
}

int TGraph::edge_weak(uint u, uint v, uint tstart, uint tend){
        return edge_interval(u, v, tstart, tend, 0);
}

int TGraph::edge_strong(uint u, uint v, uint tstart, uint tend){
        return edge_interval(u, v, tstart, tend, 1);
}

int TGraph::edge_next(uint v, uint u, uint t){
        if (v>=nodes || tgraph[v].changes == 0) return -1;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].changes];
        uint *nodep = new uint[BLOCKSIZE*tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;
        uint tnext=-1;
	uint wasOff=0;
	uint k;
        for(uint j=0; j < tgraph[v].changes; j++) {
	
                if( u == nodep[j]) { 
	                if (timep[j] > t) {
				wasOff = 1;
				k = j;
	                        break;
	                }
                        occ++;
                }
        }
        
	if (occ%2==1) {
		tnext = t;
	}
	else if (wasOff && occ%2==0) {
		tnext = timep[k];
	}

        delete [] timep;
        delete [] nodep;
        
        return tnext;
}




void TGraph::reverse_point(uint v, uint t, uint *res) {
	if (v>=nodes || reverse[v].size == 0) return;

	uint *nodep = new uint[reverse[v].size+100];

	decodereverse(v, nodep);

        uint i=0;

        for(uint j=0; j < reverse[v].size; j++) {
        	if (edge_point(nodep[j], v, t)) {
        		res[++i] = nodep[j];
        	}

        }
        *res = i;
        delete [] nodep;
}

void TGraph::reverse_weak(uint v, uint tstart, uint tend, uint *res) {
	if (v>=nodes || reverse[v].size == 0) return;

	uint *nodep = new uint[reverse[v].size+100];

	decodereverse(v, nodep);

        uint i=0;

        for(uint j=0; j < reverse[v].size; j++) {
        	if (edge_weak(nodep[j], v, tstart, tend)) {
        		res[++i] = nodep[j];
        	}

        }
        *res = i;
        delete [] nodep;
}

void TGraph::reverse_strong(uint v, uint tstart, uint tend, uint *res) {
	if (v>=nodes || reverse[v].size == 0) return;

	uint *nodep = new uint[reverse[v].size+100];

	decodereverse(v, nodep);

        uint i=0;

        for(uint j=0; j < reverse[v].size; j++) {
        	if (edge_strong(nodep[j], v, tstart, tend)) {
        		res[++i] = nodep[j];
        	}

        }
        *res = i;
        delete [] nodep;
}
