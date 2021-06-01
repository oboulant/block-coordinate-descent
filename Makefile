main: ./ebcd/main.c ./ebcd/read_data.c ./ebcd/ebcd.c
	gcc -o main ./ebcd/main.c ./ebcd/read_data.c ./ebcd/ebcd.c -I.

clean:
	rm -f main