CC		= gcc
CFLAGS= -I. -lgsl -lgslcblas -lm
VPATH	= aux
DEPS	= cubature.h
ODIR	= obj
_OBJ	= main.o hcubature.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

run: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

