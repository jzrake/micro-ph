

BIN = bin
EXE = $(BIN)/micro-ph


export CFLAGS = -Wall -O2


default : $(EXE)

$(EXE) : FORCE
	@make -C src


FORCE :
	

clean :
	make -C src clean
	rm -f $(EXE)
