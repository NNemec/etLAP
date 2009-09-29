default: sun_gen

CPPFLAGS =
CPPFLAGS += -I.
CPPFLAGS += -DNDEBUG
#CPPFLAGS += -DNO_FLATLOOP
#CPPFLAGS += -DTWODIMARRAY
#CPPFLAGS += -DDEBUG

CXXFLAGS =
# $(CPPFLAGS)
CXXFLAGS += -Wall
CXXFLAGS += -std=c++0x -pedantic
#CXXFLAGS += -save-temps
CXXFLAGS += -O3
#CXXFLAGS += -finline-limit=2000
#CXXFLAGS += -finline-limit=5000
#CXXFLAGS += -funroll-loops
#CXXFLAGS += -march=pentium
#CXXFLAGS += -march=athlon
#CXXFLAGS += -mfpmath=sse
#CXXFLAGS += -mmmx
#CXXFLAGS += -msse
#CXXFLAGS += -g
#CXXFLAGS += -fno-inline

#CXXFLAGS += -DASSIGN_POLICY=AP_norefs
CXXFLAGS += -DASSIGN_POLICY=AP_manual
#CXXFLAGS += -DASSIGN_POLICY=AP_automatic
#CXXFLAGS += -DASSIGN_POLICY=AP_conservative

#CXX = icc
CXX = g++
#CXX = mpiCC

sun_gen: sun_gen.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) sun_gen.cpp -o sun_gen

try_octave: try_octave.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) try_octave.cpp -L/usr/lib/octave-2.1.35 -loctave -lcruft -lblas -llapack -ldl -lreadline -lfftw -lkpathsea -lg2c -o try_octave

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
