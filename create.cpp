#include <cstdio>
#include <fstream>
#include "tgraphreader.h"
#include "tgraph.h"
#include "debug.h"

int main(int argc, char *argv[]) {
        uint nodes, edges, changes, maxtime;
        uint u,v,t,o;
        uint p;
        
        TGraph tg;
        
        if (argc < 2) {
                fprintf(stderr, "Usage: %s <outputfile>\n", argv[0]);
                return 1;
        }
        
        
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
        
        tg.create(tgraphreader);
        
        INFO("Saving structure...");
	ofstream f;
	f.open(argv[1], ios::binary);
	tg.save(f);
	f.close();
        
        
        
        return 0;
}
