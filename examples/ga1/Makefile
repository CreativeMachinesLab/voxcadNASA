CXX=g++
CC=g++

OPTFLAGS = -g -Wall


NaiveGA_src = \
Individual.cpp

OBJECT_FILES = \
Individual.o

all: naiveGA

naiveGA:	main.o naiveGA.o
		$(CC) $(CFLAGS) main.o -o naiveGA $(OBJECT_FILES)

main.o:		main.cpp
		$(CC) $(CFLAGS) -c main.cpp

naiveGA.o: $(NaiveGA_src)
	$(CC) $(CFLAGS) $(NaiveGA_src) -shared -fPIC -o $(OBJECT_FILES)

clean:
	rm -rf *.o naiveGA */*.o
