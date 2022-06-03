INCLUDE += -I./include

SRC_FILES = $(wildcard src/*.cpp)
OBJ_DIR := obj/
BIN_DIR := bin/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

.PHONY: all dirs clean benchmark run

all : dirs benchmark

objs : dirs $(OBJ_FILES)

dirs : $(BIN_DIR) $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) -c  $< -o  $@

benchmark: $(OBJ_FILES)
	$(CXX) $(RFLAGS) $(FLAGS) $(OBJ_FILES) $(LIBRARY) $(LDFLAGS) $(LDLIBS) -o $(BIN_DIR)/benchmark

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(BIN_DIR)*

run: dirs benchmark
	cd bin; \
	./benchmark -s -i ../pars/basic.par

stress: dirs benchmark
	cd bin; \
	time ./benchmark -i ../pars/stress_test.par
