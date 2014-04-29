#include <cstdio>
#include <fstream>

#include <map>
#include <vector>
#include <set>

#include <assert.h>
#include <getopt.h>
#include "tgraph.h"
#include "tgraphreader.h"

#include "debug.h"

struct opts {
	char *outfile;
	char c[100];
};

int readopts(int argc, char **argv, struct opts *opts) {
	int o;
	
	
	// Default options
	opts->c[0] = 0;

	while ((o = getopt(argc, argv, "c:")) != -1) {
		switch (o) {
			case 'c':
			strcpy(opts->c, optarg);
			break;
			default: /* '?' */
			break;
		}
	}
	
        if (optind >= argc || (argc-optind) < 1 || opts->c[0] == 0) {
		fprintf(stderr, "%s -c \"compression string\" <outputfile> \n", argv[0]);
		fprintf(stderr, "Expected argument for options: \n");
		fprintf(stderr, "Block Coding: pfor, turbo-rice\n");
			fprintf(stderr, "Non-Block Coding: RICE, S9, S16, VBYTE, NULL\n");
				fprintf(stderr, "Parameter -c must be a non-block code, or a blockcode:blocksize:non-blockcode:paddingsize\n");
					fprintf(stderr, "Example non-blockingcode: -c \"rice\"\n");
							fprintf(stderr, "Example block-code: -c \"pfor:128:rice:128\"\n");
//static const char kRiceCoding[] = "rice";
//static const char kTurboRiceCoding[] = "turbo-rice";
//static const char kPForDeltaCoding[] = "pfor";
//static const char kS9Coding[] = "s9";
//static const char kS16Coding[] = "s16";
//static const char kVarByteCoding[] = "vbyte";
//static const char kNullCoding[] = "null";
		
		
		fprintf(stderr, "Expected argument after options\n");
		exit(EXIT_FAILURE);
        }
	
		LOG("Compressing with '%s'",opts->c);
	
	opts->outfile = argv[optind];
	
	return optind;

}

/*
TGraphReader* readcontacts() {
	uint nodes, edges, lifetime, contacts;
	uint u,v,a,b;

	vector < map<uint, vector<uint> > > btable;

	vector < set<uint> > revgraph; //edges are in the form v,u

	scanf("%u %u %u %u", &nodes, &edges, &lifetime, &contacts);

	map<uint, vector<uint> > tm;
	btable.insert(btable.begin(), nodes, tm);
  set<uint> ts;
  revgraph.insert(revgraph.begin(), nodes, ts);

	uint c_read = 0;
	while( EOF != scanf("%u %u %u %u", &u, &v, &a, &b)) {
		c_read++;
		if(c_read%500000==0) fprintf(stderr, "Processing %.1f%%\r", (float)c_read/contacts*100);

		btable[u][a].push_back(v);


		//reverse node
		revgraph[v].insert(u);

		if (b == lifetime-1) continue;

		btable[u][b].push_back(v);
	}
	fprintf(stderr, "Processing %.1f%%\r", (float)c_read/contacts*100);
	assert(c_read == contacts);


	TGraphReader *tgraphreader = new TGraphReader(nodes,edges,2*contacts,lifetime);

	//reverse neighbors
	set<uint>::iterator its;
	for(uint i = 0; i < nodes; i++) {
		if(i%10000==0) fprintf(stderr, "Copying reverse %.1f%%\r", (float)i/nodes*100);
    
		for( its = revgraph[i].begin(); its != revgraph[i].end(); ++its) {
			tgraphreader->addReverseEdge(i, *its);
		}
    revgraph[i].clear();
	}
  vector< set<uint> >().swap(revgraph);
  revgraph.clear();





	map<uint, vector<uint> >::iterator it;

	//temporal graph
	for(uint i = 0; i < nodes; i++) {
		if(i%10000==0) fprintf(stderr, "Copying temporal %.1f%%\r", (float)i/nodes*100);
    
		for( it = btable[i].begin(); it != btable[i].end(); ++it) {
			for(uint j = 0; j < (it->second).size(); j++ ) {
				tgraphreader->addChange(i, (it->second).at(j), it->first);
			}
      vector<uint>().swap(it->second);
      (it->second).clear();
      
		}
    btable[i].clear(); 
	}



	return tgraphreader;
}

*/

TGraphReader* readcontacts() {
	uint nodes, edges, lifetime, contacts;
	uint u,v,a,b;

	vector < set<uint> > revgraph; //edges are in the form v,u

	scanf("%u %u %u %u", &nodes, &edges, &lifetime, &contacts);

  set<uint> ts;
  revgraph.insert(revgraph.begin(), nodes, ts);

	TGraphReader *tgraphreader = new TGraphReader(nodes,edges,2*contacts,lifetime);

	uint c_read = 0;
	while( EOF != scanf("%u %u %u %u", &u, &v, &a, &b)) {
		c_read++;
		if(c_read%500000==0) fprintf(stderr, "Processing %.1f%%\r", (float)c_read/contacts*100);

		tgraphreader->addChange(u, v, a);

		//reverse node
		revgraph[v].insert(u);

		if (b == lifetime-1) continue;

		tgraphreader->addChange(u, v, b);
	}
	fprintf(stderr, "Processing %.1f%%\r", (float)c_read/contacts*100);
	assert(c_read == contacts);

	//reverse neighbors
	set<uint>::iterator its;
	for(uint i = 0; i < nodes; i++) {
		if(i%10000==0) fprintf(stderr, "Copying reverse %.1f%%\r", (float)i/nodes*100);
    
		for( its = revgraph[i].begin(); its != revgraph[i].end(); ++its) {
			tgraphreader->addReverseEdge(i, *its);
		}
    revgraph[i].clear();
	}
  vector< set<uint> >().swap(revgraph);
  revgraph.clear();

	return tgraphreader;
}

int main(int argc, char *argv[]) {

        int optind;
        TGraph tg;
        struct opts opts;

        TGraphReader *tgraphreader;

	optind = readopts(argc, argv, &opts);
        
        /*        uint nodes, edges, changes, maxtime;
        uint u,v,t,o;
        uint p;
        //scanf("%d %d %d %d", &nodes, &edges, &changes, &maxtime);
        scanf("%u %u %u", &nodes,&changes, &maxtime);
        
        TGraphReader tgraphreader(nodes,edges,changes,maxtime);
        
        p = 0;
        while ( EOF != scanf("%u %u %u %u", &u, &v, &t, &o) ) {
                if (p%10000==0) fprintf(stderr, "Loading: %0.2f%%\r", (float)p/changes*100);
                p++;
                
                tgraphreader.addChange(u,v,t);
        }
        
        fprintf(stderr, "\n");
        */

	tg.set_policy(opts.c);

	tgraphreader = readcontacts();

        tg.create(*tgraphreader);
        
        INFO("Saving structure...");
	ofstream f;
	f.open(opts.outfile, ios::binary);
	tg.save(f);
	f.close();
        
        
        
        return 0;
}
