

BIN = ../bin
EXE = $(BIN)/micro-ph

default : $(EXE)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

$(EXE) : micro-ph.o solvers.o
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) -o $@ $^

clean :
	rm -f *.o
