# main: ./examples/main.c ./read_data/read_data.c ./ebcd/ebcd.c
# 	gcc -I. -Iread_data -Iebcd -o main ./examples/main.c ./read_data/read_data.c ./ebcd/ebcd.c
# clean:
# 	rm -f main

CC=gcc
CFLAGS=-I. -Iread_data -Iebcd
OBJ = ./examples/main.o ./read_data/read_data.o ./ebcd/ebcd.o
ODIR=obj

$(ODIR)/%.o: %.c %.h
	$(CC) -c -o $@ $< 

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	find . -type f -name '*.o' -exec rm {} +
	rm -f main