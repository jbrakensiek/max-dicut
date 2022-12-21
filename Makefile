CPPFLAGS=-g -Wall -Werror -O2
LDFLAGS=-g
LDLIBS=-larb -lflint

main: main.o common.o lower.o
	gcc $(LDFLAGS) -o main main.o common.o lower.o $(LDLIBS)

main.o: main.c
	gcc $(CPPFLAGS) -c main.c

common.o: common.c common.h
	gcc $(CPPFLAGS) -c common.c

lower.o: lower.c lower.h
	gcc $(CPPFLAGS) -c lower.c

.PHONY clean:
	rm -f *.o
	rm -f main

#https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile
