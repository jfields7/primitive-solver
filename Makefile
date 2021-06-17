CXX    := g++
LIB    := ar
RFLAGS := -std=c++11 -g
FLAGS  := 
LDFLAGS :=
LDLIBS := -lNumTools

INCLUDE = -I/usr/local/include
LIBS    = -L/usr/local/lib

INCLUDE += -I./include

SRC_FILES = $(wildcard src/*.cpp)
OBJ_DIR := obj/
LIB_DIR := lib/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

LIBRARY := $(LIB_DIR)PrimitiveSolver

.PHONY: all dirs clean

all : dirs $(LIBRARY)

objs : dirs $(OBJ_FILES)

dirs : $(LIB_DIR) $(OBJ_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(LIBRARY) : $(OBJ_FILES)
	$(LIB) $(RFLAGS) $(FLAGS) -o $@ $(OBJ_FILES) $(LDFLAGS) $(LDLIBS)

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) -c $< -o $@

clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(LIBRARY)
