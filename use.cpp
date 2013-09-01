/*
 * use.c
 *
 *  Created on: Mar 20, 2013
 *      Author: diegocaro
 */

#include <cstdio>
#include "arraysort.h"

#include "tgraph.h"

int main(int argc, char *argv[]) {
	TGraph *index;
	ifstream f;
	uint *buffer;
	uint i, j;
	
        if (argc < 2) {
                fprintf(stderr, "Usage: %s <inputgraph>\n", argv[0]);
                return 1;
        }

	buffer = (uint *)malloc(sizeof(uint)*BUFFER);
        
	printf("Reading file: '%s'\n", argv[1]);

	f.open(argv[1], ios::binary);
	index = TGraph::load(f);
	f.close();
        
        
      
	uint a,b;
	printf("Write your query: ");
	while(EOF != scanf("%u %u", &a, &b)) {
		index->direct_point(a, b, buffer);
		printf("t=%u -> ", 0);
		print_arraysort(buffer);
                printf("Write your query: ");
	}
	return 0;
       
/*
	uint a,b,c;
	printf("Write your query: ");
	while(EOF != scanf("%u %u %u", &a, &b, &c)) {
		index->direct_weak(a, b, c, buffer);
		printf("t=%u -> ", 0);
		print_arraysort(buffer);
                printf("Write your query: ");
	}
	return 0;
        
        
        
	uint a,b;
	printf("Write your query: ");
	while(EOF != scanf("%u %u", &a, &b)) {
		index->direct_point(a, b, buffer);
		printf("t=%u -> ", 0);
		print_arraysort(buffer);
                printf("Write your query: ");
	}
	return 0;
*/

/*
	for ( i = 0; i < 5; i++) {
		printf("start(%u): %u\n", i, start(index.map, i));
		for( j = 0; j <= index.maxtime; j++) {
			get_neighbors_point(buffer, &index, i, j);
			printf("t=%u -> ", j);print_arraysort(buffer);
		}
	}
 
	for ( i = 0; i < 5; i++) {
		printf("direct weak %u\n", i);
		for( j = 0; j <= index.maxtime; j++) {
			printf("time interval [0, %u)\n", j);
			get_neighbors_weak(buffer, &index, i, 0, j);
			printf("t=%u -> ", j);print_arraysort(buffer);
		}
	}

	for ( i = 0; i < 5; i++) {
		printf("direct strong %u\n", i);
		for( j = 0; j <= index.maxtime; j++) {
			printf("time interval [0, %u)\n", j);
			get_neighbors_strong(buffer, &index, i, 0, j);
			printf("t=%u -> ", j);print_arraysort(buffer);
		}
	}

 
	for ( i = 0; i < 5; i++) {
		printf("Reverses of node %u\n", i);
		for( j = 0; j <= 5; j++) {
			get_reverse_point(buffer, &index, i, j);
			printf("t=%u <- ", j);print_arraysort(buffer);
		}
	}
  
	for ( i = 0; i < 5; i++) {
		printf("Reverses of node %u (slow)\n", i);
		for( j = 0; j <= 5; j++) {
			get_reverse_point_slow(buffer, &index, i, j);
			printf("t=%u <- ", j);print_arraysort(buffer);
		}
	}

	for ( i = 0; i < 5; i++) {
			printf("Reverses weak %u\n", i);
			for( j = 0; j <= index.maxtime; j++) {
				printf("time interval [0, %u)\n", j);
				get_reverse_weak(buffer, &index, i, 0, j);
				printf("t=%u <- ", j);print_arraysort(buffer);
			}
		}


	for ( i = 0; i < 5; i++) {
			printf("Reverses strong %u\n", i);
			for( j = 0; j <= index.maxtime; j++) {
				printf("time interval [0, %u)\n", j);
				get_reverse_strong(buffer, &index, i, 0, j);
				printf("t=%u <- ", j);print_arraysort(buffer);
			}
		}
*/
	free(buffer);
	return 0;
}

