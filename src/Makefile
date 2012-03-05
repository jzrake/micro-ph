

BIN = ../bin
EXE = $(BIN)/micro-ph

default : $(EXE)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

$(BIN) : 
	mkdir $(BIN)

$(EXE) : micro-ph.o
	$(CC) $(CFLAGS) -o $@ $^

clean :
	rm -f *.o
