include config.mk

CXX      := g++
CPPFLAGS += $(HDF5_INCLUDE)
LIB      := ar cr
RFLAGS   := -std=c++11 -O0 -g -Wall -Wpedantic
FLAGS    := 
LDFLAGS  += $(HDF5_LIBS)

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
BENCHMARK_DIR := benchmark/
POINT_DIR := point-debugger/
FUZZER_DIR := fuzzer/
INSTALL_DIR := /usr/local

LIBRARY_NAME := libPrimitiveSolver.a

LIBRARY := $(LIB_DIR)$(LIBRARY_NAME)

.PHONY: all dirs clean

all : dirs $(LIBRARY)

objs : dirs $(OBJ_FILES)

dirs : $(LIB_DIR) $(OBJ_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(LIBRARY) : $(OBJ_FILES)
	$(LIB) $(LIBRARY) $(OBJ_FILES)

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(RFLAGS) $(FLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean
clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(LIBRARY)
	cd $(TEST_DIR) && $(MAKE) clean
	cd $(BENCHMARK_DIR) && $(MAKE) clean
	cd $(POINT_DIR) && $(MAKE) clean

.PHONY: install
install:
	mkdir -p $(INSTALL_DIR)/include/PrimitiveSolver
	mkdir -p $(INSTALL_DIR)/lib
	cp $(HEADER_FILES) $(INSTALL_DIR)/include/PrimitiveSolver
	cp $(LIBRARY) $(INSTALL_DIR)/lib

.PHONY: uninstall
uninstall:
	rm -rf $(INSTALL_DIR)/include/PrimitiveSolver
	rm -rf $(INSTALL_DIR)/lib/$(LIBRARY_NAME)

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

.PHONY: build_benchmark
build_benchmark : $(LIBRARY)
	$(MAKE) -C $(BENCHMARK_DIR)

.PHONY: benchmark
benchmark : $(LIBRARY)
	cd $(BENCHMARK_DIR) && $(MAKE) run

.PHONY: benchmark_stress
benchmark_stress : $(LIBRARY)
	cd $(BENCHMARK_DIR) && $(MAKE) stress

.PHONY: build_debugger
build_debugger : $(LIBRARY)
	$(MAKE) -C $(POINT_DIR)

.PHONY: build_fuzzer
build_fuzzer: $(LIBRARY)
	$(MAKE) -C $(FUZZER_DIR)
