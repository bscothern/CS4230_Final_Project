NAME := factor
CC := g++
CPPFLAGS := -O3
SRC_seq := factor_seq.cpp

seq:
	$(CC) $(CPPFLAGS) $(SRC_seq) -o $(NAME)_seq

clean:
	rm -f $(NAME) $(NAME)_seq
