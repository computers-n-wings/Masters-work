# Declare variables
CC = g++
CXXFLAGS = -std=c++11 -Wall
HDRS = myfunctions.h
OBJS = coursework.o myfunctions.o
LDLIBS = -llapack -lblas

# All the targets
all: compile task1 clean

# This target will compile the files
compile: $(OBJS) $(HDRS)
	$(CC) $(CXXFLAGS) $(OBJS) $(HDRS) $(LDLIBS) -o compile

# This target will run the out file
task1: compile
	./compile

# This target will remove .o files and compile executable
clean :
	-rm -f *.o compile 

%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $<

# # Object containing the main function
# coursework.o: coursework.cpp
# 	$(CC) $(CXXFLAGS) $(HDRS) $(OBJS) $(LDLIBS)

# # Object containing all of my functions
# myfunctions.o: myfunctions.cpp 
# 	$(CC) $(CXXFLAGS) $(HDRS) myfunctions.cpp $(LDLIBS)

