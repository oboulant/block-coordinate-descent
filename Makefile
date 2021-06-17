CC=gcc
CFLAGS=-I. -Iread_data -Iebcd -Wall
OBJ = ./examples/main.o ./read_data/read_data.o ./ebcd/ebcd.o
ODIR=obj

$(ODIR)/%.o: %.c %.h
	$(CC) -c -o $@ $< 

all: main

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	find . -type f -name '*.o' -exec rm {} +
	rm -f main