# Declare variables
CC = g++
CXXFLAGS = -std=c++11 -Wall
HDRS = myfunctions.h
OBJS = coursework.o myfunctions.o
LDLIBS = -llapack -lblas -lboost_program_options

# All the targets
all: compile Task2 clean

# This target will compile the files
compile: $(OBJS) $(HDRS)
	$(CC) $(CXXFLAGS) $(OBJS) $(HDRS) $(LDLIBS) -O2 -o compile

# This target will run the executable with parameters for task 1
Task1: compile
	./compile  
	python $@.py &

# This target will run the executable with parameters for task 2
Task2: compile
	./compile --Nx 24 --time_dependence 1 --scheme 0 --T 1.0 --Nt 10000 
	python Task1.py &
	# python $@.py &

# This target will run the executable with parameters for task 3
Task3: compile
	./compile --Nx 24 --time_dependence 1 --scheme 1 --T 1.0 --Nt 100
	# python Task1.py &

# This target will remove .o files and compile executable
clean:
	-rm -f *.o  compile 

%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $<

# # Object containing the main function
# coursework.o: coursework.cpp
# 	$(CC) $(CXXFLAGS) $(HDRS) $(OBJS) $(LDLIBS)

# # Object containing all of my functions
# myfunctions.o: myfunctions.cpp 
# 	$(CC) $(CXXFLAGS) $(HDRS) myfunctions.cpp $(LDLIBS)

