#  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick
#  This code is licensed under the MIT License.

# https://stackoverflow.com/questions/40621451/makefile-automatically-compile-all-c-files-keeping-o-files-in-separate-folde

SRC=src
LIB=src
BIN=bin
OBJ=bin
TEST=test
EXP=experiments
CPPFLAGS=-g -Wall -O2
LDFLAGS=-g
LDLIBS=-lflint
CC=g++

LIB_SRC=$(wildcard $(SRC)/*.cpp)
LIB_OBJ=$(patsubst $(SRC)/%.cpp, $(OBJ)/lib_%.o, $(LIB_SRC))
LIB_HPP=$(wildcard $(LIB)/*.hpp)

TEST_SRC=$(wildcard $(TEST)/*.cpp)
TEST_OBJ=$(patsubst $(TEST)/%.cpp, $(OBJ)/test_%.o, $(TEST_SRC))
TEST_BIN=$(patsubst $(TEST)/%.cpp, $(BIN)/test_%.bin, $(TEST_SRC))

EXP_SRC=$(wildcard $(EXP)/*.cpp)
EXP_HPP=$(wildcard $(EXP)/*.hpp)
EXP_OBJ=$(patsubst $(EXP)/%.cpp, $(OBJ)/exp_%.o, $(EXP_SRC))
EXP_BIN=$(patsubst $(EXP)/%.cpp, $(BIN)/exp_%.bin, $(EXP_SRC))

ALL_OBJ=$(LIB_OBJ) $(TEST_OBJ) $(EXP_OBJ)
ALL_BIN=$(TEST_BIN) $(EXP_BIN)

all: $(ALL_BIN)

$(BIN)/test_%.bin: $(OBJ)/test_%.o $(LIB_OBJ)
	$(CC) $(LDFLAGS) -o $@ $< $(LIB_OBJ) $(LDLIBS)

$(BIN)/exp_%.bin: $(OBJ)/exp_%.o $(LIB_OBJ)
	$(CC) $(LDFLAGS) -o $@ $< $(LIB_OBJ) $(LDLIBS)

$(OBJ)/lib_%.o: $(SRC)/%.cpp $(LIB_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) -c $< -o $@

$(OBJ)/test_%.o: $(TEST)/%.cpp $(LIB_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) -c $< -o $@

$(OBJ)/exp_%.o: $(EXP)/%.cpp $(LIB_HPP) $(EXP_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) -c $< -o $@

.PHONY clean:
	rm -f $(OBJ)/*.o
	rm -f $(BIN)/*.bin
