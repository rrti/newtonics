CC = g++
CFLAGS = -Wall -Wextra -g -std=c++11 -O3 -fno-strict-aliasing
LFLAGS = -lGL -lGLEW -lGLU -lSDL -lboost_chrono -lboost_system -lboost_thread

MKDIR = mkdir
TARGET = test

## BASE_SRC_DIR = ./src
BASE_SRC_DIR = ./
BASE_OBJ_DIR = ./obj
BASE_BIN_DIR = ./bin

MATH_SRC_DIR = $(BASE_SRC_DIR)/math
MATH_OBJ_DIR = $(BASE_OBJ_DIR)/math
PHYS_SRC_DIR = $(BASE_SRC_DIR)/physics
PHYS_OBJ_DIR = $(BASE_OBJ_DIR)/physics
MISC_SRC_DIR = $(BASE_SRC_DIR)/misc
MISC_OBJ_DIR = $(BASE_OBJ_DIR)/misc

BASE_OBS = \
	$(BASE_OBJ_DIR)/newtonics.o \
	$(BASE_OBJ_DIR)/world_state.o \
	$(BASE_OBJ_DIR)/program_state.o \
	$(BASE_OBJ_DIR)/program_clock.o
MATH_OBS = \
	$(MATH_OBJ_DIR)/matrix.o \
	$(MATH_OBJ_DIR)/plane.o \
	$(MATH_OBJ_DIR)/point.o \
	$(MATH_OBJ_DIR)/raw_triangle.o \
	$(MATH_OBJ_DIR)/ray.o \
	$(MATH_OBJ_DIR)/tuple.o \
	$(MATH_OBJ_DIR)/vector.o
PHYS_OBS = \
	$(PHYS_OBJ_DIR)/collision_mesh.o \
	$(PHYS_OBJ_DIR)/planar_body.o \
	$(PHYS_OBJ_DIR)/rigid_body.o
MISC_OBS = \
	$(MISC_OBJ_DIR)/opengl_wrapper.o \
	$(MISC_OBJ_DIR)/scoped_timer.o

OBJECTS = \
	$(BASE_OBS) \
	$(MATH_OBS) \
	$(PHYS_OBS) \
	$(MISC_OBS)


objects: $(OBJECTS)


$(BASE_OBJ_DIR)/%.o: $(BASE_SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c   -o $@    $<

$(MATH_OBJ_DIR)/%.o: $(MATH_SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c   -o $@    $<

$(PHYS_OBJ_DIR)/%.o: $(PHYS_SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c   -o $@    $<

$(MISC_OBJ_DIR)/%.o: $(MISC_SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c   -o $@    $<



dirs:
	if [ ! -d $(BASE_BIN_DIR) ]; then $(MKDIR) $(BASE_BIN_DIR); fi
	if [ ! -d $(BASE_OBJ_DIR) ]; then $(MKDIR) $(BASE_OBJ_DIR); fi
	if [ ! -d $(MATH_OBJ_DIR) ]; then $(MKDIR) $(MATH_OBJ_DIR); fi
	if [ ! -d $(PHYS_OBJ_DIR) ]; then $(MKDIR) $(PHYS_OBJ_DIR); fi
	if [ ! -d $(MISC_OBJ_DIR) ]; then $(MKDIR) $(MISC_OBJ_DIR); fi

target:
	## LFLAGS after OBJECTS, other way around causes linking errors
	$(CC) -o $(BASE_BIN_DIR)/$(TARGET)  $(OBJECTS)  $(LFLAGS)

all:
	make dirs
	make -j 4 objects
	make      target

clean:
	rm -f -r $(BASE_BIN_DIR)
	rm -f -r $(BASE_OBJ_DIR)

