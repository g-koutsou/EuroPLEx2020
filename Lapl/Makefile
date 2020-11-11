CC=cc
LDFLAGS_MACOS=-lomp
LDFLAGS_LINUX=-fopenmp
CFLAGS=-O3 -Xpreprocessor -fopenmp -std=gnu99
LIBS=-lm

### Tries to set appropriate LDFLAGS accorting to the result of `uname -s`
### If this fails set manually as appropriate
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	LDFLAGS=$(LDFLAGS_LINUX)
endif
ifeq ($(UNAME), Darwin)
	LDFLAGS=$(LDFLAGS_MACOS)
endif

ifneq ($(VL),)
CFLAGS += -DVL=$(VL)
endif

all:
	@echo make what?

.PHONY: A 

A: lapl
B: laplv

lapl: lapl.c 
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

laplv: laplv.c 
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	$(RM) lapl laplv
