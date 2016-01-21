CXXFLAGS=-O3 -DDSFMT_MEXP=216091 -DUSE_SYS_RAND=1

all: opt4fcm7

opt4fcm7: acor.o de.o pso.o rcga.o fcm.o main.o main_opt.o opt.o dSFMT.o rand_gen.o
	g++ $(CXXFLAGS) acor.o de.o pso.o rcga.o opt.o fcm.o main.o main_opt.o dSFMT.o rand_gen.o -o opt4fcm7

acor.o: acor.cpp
	g++ -c $(CXXFLAGS) acor.cpp

de.o: de.cpp
	g++ -c $(CXXFLAGS) de.cpp

pso.o: pso.cpp
	g++ -c $(CXXFLAGS) pso.cpp

rcga.o: rcga.cpp
	g++ -c $(CXXFLAGS) rcga.cpp

fcm.o: fcm.cpp
	g++ -c $(CXXFLAGS) fcm.cpp

main.o: main.cpp
	g++ -c $(CXXFLAGS) main.cpp

main_opt.o: main_opt.cpp
	g++ -c $(CXXFLAGS) main_opt.cpp

opt.o: opt.cpp
	g++ -c $(CXXFLAGS) opt.cpp

dSFMT.o: twister/dSFMT.c
	g++ -c $(CXXFLAGS) twister/dSFMT.c

rand_gen.o: rand_gen.cpp
	g++ -c $(CXXFLAGS) rand_gen.cpp

clean:
	rm -rf *.o opt4fcm7

