CC = gcc

BASEFLAGS = -Wall -Wextra -fopenmp
NODEBUG_FLAGS = -dNDEBUG 
DEBUG_FLAGS = -g

LDLIBS = -lm

OBJS = gene2pic.o lodepng.o NearestNeighbourUpscale.o

EXE = gene2pic

debug: CFLAGS = $(BASEFLAGS) $(DEBUG_FLAGS)
debug: $(EXE)

release: CFLAGS = $(BASEFLAGS) $(NODEBUG_FLAGS) 
release: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $(EXE) $(LDLIBS)

gene2pic.o: gene2pic.c gene2pic.h
	$(CC) $(CFLAGS) -c gene2pic.c

NearestNeighbourUpscale.o: NearestNeighbourUpscale.c NearestNeighbourUpscale.h
	$(CC) $(CFLAGS) -c NearestNeighbourUpscale.c

lodepng.o: LODEPNG/lodepng.c LODEPNG/lodepng.h
	$(CC) $(CFLAGS) -c LODEPNG/lodepng.c

clean:
	-rm -f $(OBJS)
	-rm -f *~
	-rm -f $(EXE)
	-rm -f $(EXE)_d