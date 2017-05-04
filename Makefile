NAME := factor
CC := g++
#CPPFLAGS := -g -pg file-arcs -ftest-coverage	# Test Flags
CPPFLAGS := -O3
SRC_par := factor.cpp
SRC_seq := factor_seq.cpp
SRC_broken := factor_broken.cpp
FLAGS_par := -fopenmp

all: seq par broken

par:
	$(CC) $(CPPFLAGS) $(FLAGS_par) $(SRC_par) -o $(NAME)_par

broken:
	$(CC) $(CPPFLAGS) $(FLAGS_par) $(SRC_broken) -o $(NAME)_broken

seq:
	$(CC) $(CPPFLAGS) $(SRC_seq) -o $(NAME)_seq

clean:
	rm -f $(NAME)_par $(NAME)_seq $(NAME)_broken
