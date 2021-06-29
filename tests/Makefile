INCLUDE += -I./include

SRC_FILES = $(wildcard src/*.cpp)
OBJ_DIR := obj/
BIN_DIR := bin/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

# A function for making tests
define make_test
	$(CXX) $(RFLAGS) $(FLAGS) $(INCLUDE) $(1)/$(2).cpp -c -o $(OBJ_DIR)/$(2).o
	$(CXX) $(RFLAGS) $(FLAGS) $(OBJ_FILES) $(OBJ_DIR)/$(2).o $(LIBRARY) $(LDFLAGS) $(LDLIBS) -o $(BIN_DIR)/$(2)
endef

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
	$(call make_test,eos,test_idealgas)

error_resetfloor : $(OBJ_FILE)
	$(call make_test,error,test_resetfloor)

primitive_idealgas : $(OBJ_FILE)
	$(call make_test,primitive,test_primitive_idealgas)

tests: dirs eos_idealgas error_resetfloor primitive_idealgas

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(BIN_DIR)*

run: dirs tests
	cd bin; \
	./test_idealgas; \
	./test_resetfloor; \
  ./test_primitive_idealgas
