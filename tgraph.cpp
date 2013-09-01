#include "tgraphreader.h"
#include "tgraph.h"


#include "uthash.h"
#include "etdc.h"

#include "cppUtils.h"
#include "debug.h"

#include "arraysort.h"


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
                cds_utils::saveValue<unsigned char>(out, tgraph[i].cneighbors, tgraph[i].csize_neighbors);
                
                cds_utils::saveValue<uint>(out, tgraph[i].ctime, tgraph[i].csize_time);
        }
        
        cds_utils::saveValue<uint>(out, etdctable, etdcsize);
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
                tg->tgraph[i].cneighbors = cds_utils::loadValue<unsigned char>(in, tg->tgraph[i].csize_neighbors);
                
                tg->tgraph[i].ctime = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_time);
        }
        
        tg->etdctable = cds_utils::loadValue<uint>(in, tg->etdcsize);
        
        tg->loadpolicy();
        
        return tg;
}






void TGraph::create(TGraphReader &tgr) {
        //uint *nodesbuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        //uint *timebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        unsigned char *ccnodesbuffer = new unsigned char[ (BUFFER/BLOCKSIZE+1) * BLOCKSIZE ];
        uint *cctimebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
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
                
                if (node_changes == 0) { continue;}
                 
                csize_neigh = etdc_encode(&table, tgr.tgraph[i].neighbors.data(), node_changes, ccnodesbuffer);
                tgraph[i].cneighbors = new unsigned char [csize_neigh];
                memcpy(tgraph[i].cneighbors, ccnodesbuffer, csize_neigh);
                
                encodediff(tgr.tgraph[i].timepoints);
                csize_time = cc->Compress(tgr.tgraph[i].timepoints.data(), cctimebuffer, node_changes);
                tgraph[i].ctime = new uint [csize_time];
                memcpy(tgraph[i].ctime, cctimebuffer, csize_time * sizeof(uint));
                //printf("Compression time ratio: %f\n", (float)csize_time/node_changes);


                tgraph[i].changes = node_changes;
                tgraph[i].csize_time = csize_time;
                tgraph[i].csize_neighbors = csize_neigh;
        }
        fprintf(stderr, "\n");
        
        
        
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


void TGraph::direct_point(uint v, uint t, uint *res)  {
        if (v>=nodes || tgraph[v].changes == 0) return;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint i=0;
        
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] > t) break;
                res[++i] = nodep[j];
        }
        
        *res = i;
        qsort(&res[1], *res, sizeof(unsigned int), compare);
        xor_arraysort(res);
}

void TGraph::direct_weak(uint v, uint tstart, uint tend, uint *res)  {
        if (v>=nodes|| tgraph[v].changes == 0) return;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        uint *buffer = new uint[tgraph[v].changes];
        uint *interval = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
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
        
        merge_arraysort(res, buffer, interval);
}


void TGraph::direct_strong(uint v, uint tstart, uint tend, uint *res)  {
        if (v>=nodes || tgraph[v].changes == 0) return;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        //uint *buffer = new uint[tgraph[v].changes];
        uint *interval = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint i=0;
        uint k=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] <= tstart) {
                        res[++i] = nodep[j];
                }
                
                else if (timep[j] > tstart && timep[j]<= tend) {
                        interval[++k] = nodep[j];
                }
                
                if (timep[j] > tend) break;
        }
        *res = i;
        *interval = k;
        
        qsort(&res[1], *res, sizeof(unsigned int), compare);
        xor_arraysort(res);
        
        qsort(&interval[1], *interval, sizeof(unsigned int), compare);
        remove_duplicates(interval);
        
        diff_arraysort(res, interval);
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

int TGraph::edge_point(uint u, uint v, uint t){
        if (v>=nodes || tgraph[v].changes == 0) return 0;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] > t) break;
                if( v == nodep[j]) occ++;
        }
        
        return (occ%2);
}

int TGraph::edge_weak(uint u, uint v, uint tstart, uint tend){
        if (v>=nodes || tgraph[v].changes == 0) return 0;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;
        uint occinterval=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] <= tstart) {
                        if( v == nodep[j]) occ++;
                }
                else if (timep[j] > tstart && timep[j]<= tend) {
                        if( v == nodep[j]) occinterval++;
                }
                
                if (timep[j] > tend) break;
                
        }
        
        if ( (occ%2) || occinterval > 0) {
                return 1;
        }
        return 0;
        
        
}

int TGraph::edge_strong(uint u, uint v, uint tstart, uint tend){
        if (v>=nodes || tgraph[v].changes == 0) return 0;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;
        uint occinterval=0;
        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] <= tstart) {
                        if( v == nodep[j]) occ++;
                }
                else if (timep[j] > tstart && timep[j]<= tend) {
                        if( v == nodep[j]) occinterval++;
                }
                
                if (timep[j] > tend) break;
                
        }
        
        if ( (occ%2) && occinterval > 0 ) {
                return 1;
        }
        return 0;
}

int TGraph::edge_next(uint u, uint v, uint t){
        if (v>=nodes || tgraph[v].changes == 0) return 0;
        
        uint *timep = new uint[tgraph[v].changes];
        uint *nodep = new uint[tgraph[v].changes];
        
        decodetime(v, timep);
        decodeneigh(v, nodep);
        
        uint occ=0;

        for(uint j=0; j < tgraph[v].changes; j++) {
                if (timep[j] > t) {
                        if (occ%2 == 1) return t;
                        
                }
                
                if( v == nodep[j]) occ++;
        }
        
        return -1;
}