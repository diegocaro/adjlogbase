#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <map>

#include "arraysort.h"
#include "timing.h"

using namespace std;

void generateRandom(uint *seq, uint n, uint sigma) {
  for(uint i =0; i < n; i++) {
    seq[i] = rand()%sigma;
  }
}

#define BUFFER 10000000

void copyBuffer(uint *buffer, uint n, uint *seq) {
  for(uint k=0; k < n; k++) {
    buffer[k+1] = seq[k];
  }
  *buffer = n;
}


// O( n log n )
void countSort(const uint *seq, uint *res) {
  for(uint j=1; j <= *seq; j++) {
   res[j] = seq[j];
  }
  *res = *seq;
  
  qsort(&res[1], *res, sizeof(unsigned int), compare);
  xor_arraysort(res);
}

// O( n )
void countHash(const uint *seq, uint *res) {
  unordered_map <uint,uint> thash;
  for(uint j=1; j <= *seq; j++) {
    thash[seq[j]]++;
  }
        
  uint i=0;
  for ( unordered_map<uint,uint>::iterator it = thash.begin(); it != thash.end(); ++it ) {
    if (it->second %2 == 1) {
      res[++i] = it->first;
    }
  }
  *res = i;
}

// n log (sigma)
void countMap(const uint *seq, uint *res) {
  map <uint,uint> tmap;
  for(uint j=1; j <= *seq; j++) {
    tmap[seq[j]]++;
  }
  uint i=0;
  for ( map<uint,uint>::iterator it = tmap.begin(); it != tmap.end(); ++it ) {
    if (it->second %2 == 1) {
      res[++i] = it->first;
    }
  }
  *res = i;
}

// n log (sigma)
void countSet(const uint *seq, uint *res) {
  set <uint> tset;
  set<uint>::iterator it;

  for(uint j=1; j <= *seq; j++) {
    it = tset.find(seq[j]);
    
    if (it == tset.end()) {
      tset.insert(seq[j]);
    }
    else {
      tset.erase(it);
    }
      
  }
  
  uint i=0;
  for ( set<uint>::iterator it = tset.begin(); it != tset.end(); ++it ) {
    res[++i] = *it;
  }
  *res = i;
  
}


// O(n)
void countUnsrtSet(const uint *seq, uint *res) {
  unordered_set <uint> tset;
  unordered_set<uint>::iterator it;

  for(uint j=1; j <= *seq; j++) {
    it = tset.find(seq[j]);
    
    if (it == tset.end()) {
      tset.insert(seq[j]);
    }
    else {
      tset.erase(it);
    }
      
  }
  
  uint i=0;
  for ( unordered_set<uint>::iterator it = tset.begin(); it != tset.end(); ++it ) {
    res[++i] = *it;
  }
  *res = i;
  
}


int main() {
  uint nsigmas=16;
  uint sigmas[16] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
  
  uint N = 21;
  uint Ns[21] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152};
  
  uint samples = 4;
  
  srand (time(NULL));
  uint *seq = new uint[BUFFER];
  uint *buffer = new uint[BUFFER];

  unsigned long elapsed;
  uint n, sigma;
  for(uint i=0; i < N; i++) {

    for(uint j=0; j < nsigmas; j++) {
      n = Ns[i];
      sigma = sigmas[j];
      
      if (sigma > n) continue;

      generateRandom(&seq[1], n, sigma);
      *seq = n;
      
      for(uint k=0; k < samples; k++) {
      
    	startClockTime();
      countSort(seq, buffer);
      elapsed = endClockTime();
      printf("%s\t%u\t%u\t%ld\t%ld\n","countSort", n, sigma, elapsed, elapsed / n);

    	startClockTime();
      countHash(seq, buffer);
      elapsed = endClockTime();
      printf("%s\t%u\t%u\t%ld\t%ld\n","countHash", n, sigma, elapsed, elapsed / n);
      
    	startClockTime();
      countMap(seq, buffer);
      elapsed = endClockTime();
      printf("%s\t%u\t%u\t%ld\t%ld\n","countMap", n, sigma, elapsed, elapsed / n);

    	startClockTime();      
      countSet(seq, buffer);
      elapsed = endClockTime();
      printf("%s\t%u\t%u\t%ld\t%ld\n","countSet", n, sigma, elapsed, elapsed / n);
      
    	startClockTime();      
      countUnsrtSet(seq, buffer);
      elapsed = endClockTime();
      printf("%s\t%u\t%u\t%ld\t%ld\n","countUnsrtSet", n, sigma, elapsed, elapsed / n);
      
    }
      
    }
    
    
  }
  delete [] seq;
  delete [] buffer;
  return 0;
}