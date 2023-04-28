#!/bin/bash

# # Variables:
# # =====================================================================================================================
# COMPILER= g++
# SYS=$(shell uname)
# LIB_ROOT = $(PWD)/
# ARMA_INCL = $(LIB_ROOT)arma_libs/include/
# ARMA_LIBS = $(LIB_ROOT)arma_libs/lib/
# INCL = -I $(ARMA_INCL) -I include/
# # LIBS = -L $(ARMA_LIBS) -larmadillo
# ifeq ($(shell uname),Darwin)
#     LIBS = -L $(ARMA_LIBS) -larmadillo -headerpad_max_install_names
# else
#     LIBS = -L $(ARMA_LIBS) -larmadillo
# endif
# OBJ = main_1.o BinaryTree.o
# OPT = -g
#
# # =====================================================================================================================
# all: bin/main_1.exe
#
# bin/main_1.exe: obj/main_1.o obj/BinaryTree.o
#
# 	$(COMPILER) -o $@ $^ $(LIBS)
#
# 	if [ $(SYS) = "Darwin" ]; then \
# 		echo "Additional steps for compilation on OS... DONE!"; \
# 		install_name_tool -change @rpath/libarmadillo.9.dylib $(ARMA_LIBS)libarmadillo.9.dylib $@; \
# 	fi
#
# obj/main_1.o: src/main_1.cpp
# 	$(COMPILER) $(OPT) -c $^ -o $@ $(INCL) -std=c++11
#
# obj/BinaryTree.o: src/BinaryTree.cpp include/BinaryTree.h
# 	$(COMPILER) $(OPT) -c $< -o $@ $(INCL) -std=c++11
#
# clean: cleanBin cleanObj
#
# cleanBin:
# 	-rm -r bin/*
#
# cleanObj:
# 	-rm -r obj/*

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
SRCS = main_1.cpp main_2.cpp main_3.cpp main_4.cpp
OBJS = $(SRCS:%.cpp=obj/%.o)
EXES = $(SRCS:%.cpp=bin/%.exe)
OPT = -g

# =====================================================================================================================
all: $(EXES)

bin/%.exe: obj/%.o obj/BinaryTree.o
	$(COMPILER) -o $@ $^ $(LIBS)

	if [ $(SYS) = "Darwin" ]; then \
		echo "Additional steps for compilation on OS... DONE!"; \
		install_name_tool -change @rpath/libarmadillo.9.dylib $(ARMA_LIBS)libarmadillo.9.dylib $@; \
	fi

obj/%.o: src/%.cpp
	$(COMPILER) $(OPT) -c $< -o $@ $(INCL) -std=c++11

obj/BinaryTree.o: src/BinaryTree.cpp include/BinaryTree.h
	$(COMPILER) $(OPT) -c $< -o $@ $(INCL) -std=c++11

clean:
	rm -rf obj/* bin/*

.SUFFIXES:
