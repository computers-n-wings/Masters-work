# Declare variables
CCP = mpicxx
CXXFLAGS = -std=c++11 -Wall -O2
HDRS = myfunctions.h solvers.h
OBJS = coursework.o myfunctions.o solvers.o
LDLIBSP = -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas -lboost_program_options

# Default target
all: Task5 clean

# This target will compile the files
compile: $(OBJS) $(HDRS)
	$(CCP) $(CXXFLAGS) $(OBJS) $(HDRS) $(LDLIBSP) -o run

# This target will run the executable with parameters for task 1
Task1: compile
	mpirun -np 1 ./run 
	python $@.py &

# This target will run the executable with parameters for task 2
Task2: compile
	mpirun -np 1 ./run --time_dependence 1 --scheme 0
	# python $@.py &

# This target will run the executable with parameters for task 3
Task3: compile
	mpirun -np 1 ./run --time_dependence 1 --scheme 1
	# python Task1.py &

Task4: compile
	mpirun -np 2 ./run --time_dependence 1 --scheme 0 --parallel 1

Task5: compile
	mpirun -np 4 ./run --Nx 24 --time_dependence 1 --scheme 1 --parallel 1 --Nt 10000

# This target will remove .o files and compile executable
clean:
	-rm -f *.o run

%.o: %.cpp $(HDRS)
	$(CCP) $(CXXFLAGS) -o $@ -c $<

