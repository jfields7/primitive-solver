INCLUDE += -I./include

SRC_FILES = $(wildcard src/*.cpp)
OBJ_DIR := obj/
BIN_DIR := bin/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

.PHONY: all dirs clean fuzzer

all : dirs fuzzer

objs : dirs $(OBJ_FILES)

dirs : $(BIN_DIR) $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) -c  $< -o  $@

fuzzer: $(OBJ_FILES)
	$(CXX) $(RFLAGS) $(FLAGS) $(OBJ_FILES) $(LIBRARY) $(LDFLAGS) $(LDLIBS) -o $(BIN_DIR)/fuzzer

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(BIN_DIR)/fuzzer
