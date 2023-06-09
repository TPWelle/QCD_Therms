CC= g++ -g -std=c++11
OMPFLAGS = -fopenmp
GSLFLAGS = -lgsl -lgslcblas -lm
DIR=$(shell pwd)
#-L/usr/local/lib -L/usr/src/debug

IFLAGS =   -I include
LFLAGS =  $(OMPFLAGS) $(GSLFLAGS) -I/usr/local/include
CFLAGS = -g -O3 -fPIC $(OMPFLAGS) $(IFLAGS) -fpermissive -fvisibility=hidden
EXECUTABLE = program

IFLAGS_PYBIND = $(shell python3 -m pybind11 --includes)
OBJECTS := $(patsubst src/%.cpp, obj/%.o,$(wildcard src/*.cpp))


compute : therm.so
	python3 -W ignore compute.py

opti : therm.so
	python3 -W ignore opti.py

obj/%.o : src/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS_PYBIND) -c -o $@ $^

therm.so : $(OBJECTS)
	$(CC) $(CFLAGS) -shared -Wl,-soname,therm $(IFLAGS_PYBIND) $^ -o therm.so $(LFLAGS)

clean:
	rm -f *.o *.so *.a
	rm -f -r obj/*.o


