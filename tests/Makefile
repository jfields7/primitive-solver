INCLUDE += -I./include

SRC_FILES = $(wildcard src/*.cpp)
OBJ_DIR := obj/
BIN_DIR := bin/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

.PHONY: all dirs clean tests run

all : dirs tests

objs : dirs $(OBJ_FILES)

dirs : $(BIN_DIR) $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) -c $< -o $@

eos_idealgas : $(OBJ_FILES)
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) eos/test_idealgas.cpp -c -o $(OBJ_DIR)/test_idealgas.o
	$(CXX) $(RFLAGS) $(FLAGS) $(OBJ_FILES) $(OBJ_DIR)/test_idealgas.o $(LIBRARY) -o $(BIN_DIR)/test_idealgas

error_resetfloor : $(OBJ_FILE)
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) error/test_resetfloor.cpp -c -o $(OBJ_DIR)/test_resetfloor.o
	$(CXX) $(RFLAGS) $(FLAGS) $(OBJ_FILES) $(OBJ_DIR)/test_resetfloor.o $(LIBRARY) -o $(BIN_DIR)/test_resetfloor

tests: dirs eos_idealgas error_resetfloor

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(BIN_DIR)*

run: dirs tests
	cd bin; \
	./test_idealgas; \
	./test_resetfloor
