CC=		cc
CFLAGS=		-g -Wall -O2 -Wno-unused-function
CPPFLAGS=
INCLUDES=
PROG=		amplicon-tk
BINDIR=		/usr/local/bin
LIBS=		-lz
OBJS=		uniques.o patch.o bin.o collapse.o mapping.o level.o lca.o voting.o kstring.o utils.o

.SUFFIXES:.c .o
.PHONY:all clean

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

amplicon-tk:$(OBJS) amplicon-tk.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

install:all
		install amplicon-tk $(BINDIR)

clean:
		rm -fr *.o amplicon-tk