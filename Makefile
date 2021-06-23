CXX    := g++
LIB    := ar cr
RFLAGS := -std=c++11 -g
FLAGS  := 
LDFLAGS :=
#LDLIBS := -lNumTools
LDLIBS := 

INCLUDE = -I/usr/local/include
LIBS    = -L/usr/local/lib

INCLUDE += -I$(CURDIR)/include

SRC_FILES = $(wildcard src/*.cpp)
HEADER_FILES = $(wildcard include/*.hpp)
OBJ_DIR := obj/
LIB_DIR := $(CURDIR)/lib/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)
TEST_DIR := tests/
INSTALL_DIR := /usr/local

LIBRARY := $(LIB_DIR)libPrimitiveSolver.a

.PHONY: all dirs clean

all : dirs $(LIBRARY)

objs : dirs $(OBJ_FILES)

dirs : $(LIB_DIR) $(OBJ_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(LIBRARY) : $(OBJ_FILES)
	$(LIB) $(LIBRARY) $(OBJ_FILES) $(LDFLAGS) $(LDLIBS)

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean
clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(LIBRARY)
	cd $(TEST_DIR) && $(MAKE) clean

.PHONY: install
install:
	mkdir -p $(INSTALL_DIR)/include/PrimitiveSolver
	cp $(HEADER_FILES) $(INSTALL_DIR)/include/PrimitiveSolver
	cp $(LIBRARY) $(INSTALL_DIR)/lib

# Variables to export for making tests.
export CXX
export RFLAGS
export FLAGS
export LDFLAGS
export INCLUDE
export LIBS
export LIB_DIR
export LIBRARY

.PHONY: build_tests
build_tests : $(LIBRARY)
	$(MAKE) -C $(TEST_DIR)

.PHONY: test
test : $(LIBRARY)
	cd $(TEST_DIR) && $(MAKE) run
