DEPEN_COMPRESS=../compresslists/
DEPEN_COMPRESS2=../compresslists/compression_toolkit/

DEPEN_ETDC=../etdc/

CC=g++
DEBUG = -g -DDEBUG

EXPERIMENTS = -DEXPERIMENTS

CFLAGS=-c -Wall -I$(DEPEN_COMPRESS) -I$(DEPEN_COMPRESS2) -I$(DEPEN_ETDC) -O3 
LDFLAGS=-Wall -lm 

EXT=

SOURCES=tgraph.cpp arraysort.cpp
OBJECTS=$(SOURCES:.cpp=.o)

MAINSRC=create.cpp use.cpp benchmark.cpp exp-query.cpp #use.cpp getsize.cpp
MAINOBJ=$(MAINSRC:.cpp=.o)

COMPRESSOBJS=compression_toolkit/coding_factory.o compression_toolkit/pfor_coding.o compression_toolkit/rice_coding2.o compression_toolkit/s9_coding.o compression_toolkit/vbyte_coding.o coding_policy.o compression_toolkit/null_coding.o compression_toolkit/rice_coding.o compression_toolkit/s16_coding.o compression_toolkit/unpack.o
DEPENOBJS+=$(addprefix $(DEPEN_COMPRESS), $(COMPRESSOBJS))

ETDCOBJS=etdc.o
DEPENOBJS+=$(addprefix $(DEPEN_ETDC), $(ETDCOBJS))

EXECUTABLE=create use benchmark exp-query #getsize text2adjacencybin

#all: CFLAGS += -O9
all: compresslists etdc $(SOURCES) $(MAINSRC) $(MAINOBJ) $(OBJECTS) $(EXECUTABLE) 

compresslists: 
	$(MAKE) -C $(DEPEN_COMPRESS)

etdc: 
	$(MAKE) -C $(DEPEN_ETDC)

debug: CFLAGS += $(DEBUG)
debug: clean all

experiments: CFLAGS += $(EXPERIMENTS)
experiments: clean all

#$(EXECUTABLE): $(MAINOBJ) $(OBJECTS)
#	$(CC) $(LDFLAGS) $(OBJECTS) $(DEPENOBJS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: cleanobj cleanexe

cleanobj:
	rm -f $(OBJECTS) $(MAINOBJ)

cleanexe:
	rm -f $(EXECUTABLE)

create: create.o
	$(CC) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

use: use.o tgraph.o
	$(CC) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

benchmark: benchmark.o
	$(CC) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

getsize: getsize.o
	$(CC) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

text2adjacencybin: text2adjacencybin.o
	$(CC) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)
	
exp-query: exp-query.o
	$(CC) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)