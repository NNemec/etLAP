all: exp

CPPFLAGS =
CPPFLAGS += -I.
#CPPFLAGS += -DNDEBUG
#CPPFLAGS += -DFLATLOOP
#CPPFLAGS += -DTWODIMARRAY
#CPPFLAGS += -DDEBUG

CXXFLAGS = $(CPPFLAGS)
CXXFLAGS += -Wall
#CXXFLAGS += -save-temps
CXXFLAGS += -O3
CXXFLAGS += -finline-limit=2000
#CXXFLAGS += -finline-limit=5000
CXXFLAGS += -funroll-loops
CXXFLAGS += -march=pentium
#CXXFLAGS += -march=athlon
#CXXFLAGS += -g

#CXXFLAGS += -DASSIGN_POLICY=AP_norefs
CXXFLAGS += -DASSIGN_POLICY=AP_manual
#CXXFLAGS += -DASSIGN_POLICY=AP_automatic
#CXXFLAGS += -DASSIGN_POLICY=AP_conservative

#CXX = g++-3.0
CXX = g++-3.2
#CXX = g++
#CXX = mpiCC

exp: exp.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) exp.cpp -o exp 

benchmark: exp
	@uname -a >> benchmarks
	@date >> benchmarks
	@echo $(CXX) $(CXXFLAGS) >> benchmarks
	@exp | tee -a benchmarks
	@echo >> benchmarks

fsimple: fsimple.cpp etLAP/*.h Makefile
	$(TIME) $(CXX) $(CXXFLAGS) fsimple.cpp -o fsimple 

tryout: tryout.cpp etLAP/*.h Makefile
	$(TIME) $(CXX) $(CXXFLAGS) tryout.cpp -o tryout

clean:
	rm -f *.ii *.o *.s
	rm -f fexp vexp fsimple tryout
