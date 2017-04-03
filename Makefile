NAME := factor
CC := g++
CPPFLAGS := -O3
SRC_par := factor.cpp
SRC_seq := factor_seq.cpp

all: seq par

par:
	$(CC) $(CPPFLAGS) $(SRC_par) -o $(NAME)_par

seq:
	$(CC) $(CPPFLAGS) $(SRC_seq) -o $(NAME)_seq

clean:
	rm -f $(NAME)_par $(NAME)_seq
