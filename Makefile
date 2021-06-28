CC=gcc
CFLAGS=-Wall
INC_EBCD_FLAGS := -Iebcd
INC_READ_FLAGS := -Iread_data
INC_UNITY_FLAGS := -Itests/unity/src
OBJ = ./read_data/read_data.o ./ebcd/ebcd.o
LIBRARIES = -lm


%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

all: main

main: $(OBJ) ./examples/main.o 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBRARIES) 
	
./examples/main.o:
	$(CC) -c -o ./examples/main.o ./examples/main.c $(INC_EBCD_FLAGS) $(INC_READ_FLAGS) $(CFLAGS) 

test: $(OBJ_TEST)
	$(CC) -c -o ./tests/TestEbcd.o ./tests/TestEbcd.c $(INC_UNITY_FLAGS) $(INC_EBCD_FLAGS) $(CFLAGS)
	$(CC) -c -o ./tests/unity/src/unity.o ./tests/unity/src/unity.c $(INC_UNITY_FLAGS) $(CFLAGS)
	$(CC) -o ./test ./tests/TestEbcd.o ./tests/unity/src/unity.o $(CFLAGS) $(LIBRARIES) 

.PHONY: clean

clean:
	find . -type f -name '*.o' -exec rm {} +
	rm -f main
	rm -f test