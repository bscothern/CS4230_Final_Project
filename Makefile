NAME := factor
CC := g++
<<<<<<< HEAD
CPPFLAGS := -O3
=======
CPPFLAGS := -g -pg #-O3
>>>>>>> c346a9c24f421a118017d4616a658fef56dc2100
SRC_par := factor.cpp
SRC_seq := factor_seq.cpp
SRC_mike := mike_factor.cpp
FLAGS_par := -fopenmp

all: seq par mike

par:
	$(CC) $(CPPFLAGS) $(FLAGS_par) $(SRC_par) -o $(NAME)_par

mike:
	$(CC) $(CPPFLAGS) $(FLAGS_par) $(SRC_mike) -o $(NAME)_mike

seq:
	$(CC) $(CPPFLAGS) $(SRC_seq) -o $(NAME)_seq

clean:
	rm -f $(NAME)_par $(NAME)_seq
