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

eos_eoscompose : $(OBJ_FILES)
	$(call make_test,eos,test_eoscompose)

eos_idealgas : $(OBJ_FILES)
	$(call make_test,eos,test_idealgas)

eos_piecewise : $(OBJ_FILES)
	$(call make_test,eos,test_piecewise)

coldeos_eoscompose : $(OBJ_FILES)
	$(call make_test,coldeos,test_coldeoscompose)

coldeos_polytrope : $(OBJ_FILES)
	$(call make_test,coldeos,test_polytrope)

coldeos_pwp : $(OBJ_FILES)
	$(call make_test,coldeos,test_pwp)

error_resetfloor : $(OBJ_FILE)
	$(call make_test,error,test_resetfloor)

primitive_idealgas : $(OBJ_FILE)
	$(call make_test,primitive,test_primitive_idealgas)

primitive_piecewise : $(OBJ_FILE)
	$(call make_test,primitive,test_primitive_piecewise)

primitive_eoscompose : $(OBJ_FILE)
	$(call make_test,primitive,test_primitive_eoscompose)

primitive_unphysical : $(OBJ_FILE)
	$(call make_test,primitive,test_primitive_unphysical)

unit_systems : $(OBJ_FILE)
	$(call make_test,units,test_units)

tests: dirs coldeos_pwp coldeos_eoscompose coldeos_polytrope eos_eoscompose eos_idealgas eos_piecewise error_resetfloor primitive_idealgas primitive_piecewise primitive_eoscompose primitive_unphysical unit_systems

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(BIN_DIR)*

run: dirs tests
	cd bin; \
	./test_eoscompose; \
	./test_idealgas; \
	./test_piecewise; \
	./test_resetfloor; \
	./test_primitive_idealgas; \
	./test_primitive_piecewise; \
	./test_primitive_eoscompose; \
	./test_primitive_unphysical; \
	./test_units
