all: sun_gen exp fsimple benchmark # try_octave  try_gsl

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

test_vector: test_vector.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

sun_gen: sun_gen.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

try_octave: try_octave.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ -L/usr/lib/octave-2.1.35 -loctave -lcruft -lblas -llapack -ldl -lreadline -lfftw -lkpathsea -lg2c

try_gsl: try_gsl.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ -lgsl -lcblas -latlas

exp: exp.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

benchmark: exp
	@uname -a >> benchmarks
	@date >> benchmarks
	@echo $(CXX) $(CXXFLAGS) $(CPPFLAGS) >> benchmarks
	@exp | tee -a benchmarks
	@echo >> benchmarks

fsimple: fsimple.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

tryout: tryout.cpp etLAP/*.h Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.ii *.o *.s
	rm -f exp fsimple tryout try_gsl
