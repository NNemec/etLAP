all: try_gsl

CPPFLAGS =
CPPFLAGS += -I.
CPPFLAGS += -DNDEBUG
#CPPFLAGS += -DNO_FLATLOOP
#CPPFLAGS += -DTWODIMARRAY
#CPPFLAGS += -DDEBUG

CXXFLAGS =
# $(CPPFLAGS)
CXXFLAGS += -Wall
#CXXFLAGS += -save-temps
CXXFLAGS += -O3
CXXFLAGS += -finline-limit=2000
#CXXFLAGS += -finline-limit=5000
CXXFLAGS += -funroll-loops
#CXXFLAGS += -march=pentium
CXXFLAGS += -march=athlon
CXXFLAGS += -mfpmath=sse
CXXFLAGS += -mmmx
CXXFLAGS += -msse
#CXXFLAGS += -g
#CXXFLAGS += -fno-inline

#CXXFLAGS += -DASSIGN_POLICY=AP_norefs
CXXFLAGS += -DASSIGN_POLICY=AP_manual
#CXXFLAGS += -DASSIGN_POLICY=AP_automatic
#CXXFLAGS += -DASSIGN_POLICY=AP_conservative

#CXX = g++-3.0
CXX = g++-3.2
#CXX = g++
#CXX = mpiCC

try_gsl: try_gsl.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) try_gsl.cpp -lgsl -lcblas -latlas -o try_gsl

exp: exp.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) exp.cpp -o exp 

benchmark: exp
	@uname -a >> benchmarks
	@date >> benchmarks
	@echo $(CXX) $(CXXFLAGS) $(CPPFLAGS) >> benchmarks
	@exp | tee -a benchmarks
	@echo >> benchmarks

fsimple: fsimple.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) fsimple.cpp -o fsimple 

tryout: tryout.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) tryout.cpp -o tryout

clean:
	rm -f *.ii *.o *.s
	rm -f exp fsimple tryout try_gsl
