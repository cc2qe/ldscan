# define our object and binary directories
export BIN_DIR  = bin
export SRC_DIR  = src
export C = gcc
export CFLAGS = -O3

all: 
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)
	$(C) -o $(BIN_DIR)/ldscan $(SRC_DIR)/ldscan.c $(CFLAGS)
