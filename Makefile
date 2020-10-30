CC=gcc
LDFLAGS_MACOS=-lomp
LDFLAGS_LINUX=-fopenmp
CFLAGS=-Xpreprocessor -fopenmp

### Tries to set appropriate LDFLAGS accorting to the result of `uname -s`
### If this fails set manually as appropriate
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	LDFLAGS=$(LDFLAGS_LINUX)
endif
ifeq ($(UNAME), Darwin)
	LDFLAGS=$(LDFLAGS_MACOS)
endif


all: runme

%.o: %.c
	$(CC) $(CFLAGS) -c $<

runme: runme.o
	$(CC) -o $@ $< $(LDFLAGS)

clean:
	$(RM) runme.o
