# =============================================================================
# Makefile — Instant Q Engine
# =============================================================================
CC      = g++
CFLAGS  = -O3 -march=native -fopenmp -std=c++17 -Wall
LDFLAGS = -lz -fopenmp

PREFIX  ?= $(HOME)/.local

SRC     = instant_q.cpp
BIN     = instant_q

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

install: $(BIN)
	mkdir -p $(PREFIX)/bin
	cp $(BIN) $(PREFIX)/bin/

clean:
	rm -f $(BIN)

test:
	@echo "Compile test..."
	$(CC) $(CFLAGS) -o $(BIN) $(SRC) $(LDFLAGS)
	./$(BIN) --help
	@echo "OK"

.PHONY: all install clean test
