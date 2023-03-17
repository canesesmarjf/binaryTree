#!/bin/bash

# Variables:
# =====================================================================================================================
COMPILER= g++
SYS=$(shell uname)
LIB_ROOT = $(PWD)/
ARMA_INCL = $(LIB_ROOT)arma_libs/include/
ARMA_LIBS = $(LIB_ROOT)arma_libs/lib/
INCL = -I $(ARMA_INCL) -I include/
# LIBS = -L $(ARMA_LIBS) -larmadillo
ifeq ($(shell uname),Darwin)
    LIBS = -L $(ARMA_LIBS) -larmadillo -headerpad_max_install_names
else
    LIBS = -L $(ARMA_LIBS) -larmadillo
endif
OBJ = main.o BinaryTree.o

# =====================================================================================================================
all: bin/BinarySearch.exe

bin/BinarySearch.exe: obj/main.o obj/BinaryTree.o

	$(COMPILER) -o $@ $^ $(LIBS)

	if [ $(SYS) = "Darwin" ]; then \
		echo "Additional steps for compilation on OS... DONE!"; \
		install_name_tool -change @rpath/libarmadillo.9.dylib $(ARMA_LIBS)libarmadillo.9.dylib $@; \
	fi

obj/main.o: src/main.cpp
	$(COMPILER) -c $^ -o $@ $(INCL) -std=c++11

obj/BinaryTree.o: src/BinaryTree.cpp include/BinaryTree.h
	$(COMPILER) -c $< -o $@ $(INCL) -std=c++11

clean: cleanBin cleanObj

cleanBin:
	-rm -r bin/*

cleanObj:
	-rm -r obj/*
