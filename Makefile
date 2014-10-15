DEPEN_COMPRESS=../compresslists/
DEPEN_COMPRESS2=../compresslists/compression_toolkit/

DEPEN_ETDC=../etdc/

DEBUG = -g -DDEBUG

EXPERIMENTS = -DEXPERIMENTS

UNAME := $(shell uname)
LDFLAGSRT = 
ifeq ($(UNAME), Linux)
LDFLAGSRT = -lrt
endif

CFLAGS=-c -Wall -I$(DEPEN_COMPRESS) -I$(DEPEN_COMPRESS2) -I$(DEPEN_ETDC) -O3 -std=c++0x
LDFLAGS=-Wall -lm $(LDFLAGSRT)

EXT=

SOURCES=tgraph.cpp arraysort.cpp
OBJECTS=$(SOURCES:.cpp=.o)

MAINSRC=create.cpp use.cpp benchmark.cpp exp-query.cpp benchmark2.cpp #use.cpp getsize.cpp
MAINOBJ=$(MAINSRC:.cpp=.o)

COMPRESSOBJS=compression_toolkit/coding_factory.o compression_toolkit/pfor_coding.o compression_toolkit/rice_coding2.o compression_toolkit/s9_coding.o compression_toolkit/vbyte_coding.o coding_policy.o compression_toolkit/null_coding.o compression_toolkit/rice_coding.o compression_toolkit/s16_coding.o compression_toolkit/unpack.o
DEPENOBJS+=$(addprefix $(DEPEN_COMPRESS), $(COMPRESSOBJS))

ETDCOBJS=etdc.o
DEPENOBJS+=$(addprefix $(DEPEN_ETDC), $(ETDCOBJS))

EXECUTABLE=create use benchmark exp-query benchmark2 #getsize text2adjacencybin

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
	$(CXX) $(CFLAGS) $< -o $@

clean: cleanobj cleanexe

cleanobj:
	rm -f $(OBJECTS) $(MAINOBJ)

cleanexe:
	rm -f $(EXECUTABLE)

create: create.o
	$(CXX) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

use: use.o tgraph.o
	$(CXX) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

benchmark: benchmark.o
	$(CXX) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

benchmark2: benchmark2.o
	$(CXX) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

getsize: getsize.o
	$(CXX) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)

exp-query: exp-query.o
	$(CXX) $(OBJECTS) $(DEPENOBJS) $< -o $@$(EXT) $(LDFLAGS)
