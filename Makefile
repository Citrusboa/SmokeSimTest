# Compiler and target
CC = g++
TARGET = main

# Release / Debug flags
RELEASE := 1
ifeq ($(RELEASE),1)
    CXXFLAGS += -O3
else
    CXXFLAGS += -g
endif

# Include paths
INCLUDE_PATH = -I/mingw64/include

# Linker libraries
LIBS = -lglew32 -lglfw3 -lopengl32 -lglu32

# Eigen is header-only — nothing to link
CXXFLAGS += -I/mingw64/include/eigen3

# Source files
SRC := $(shell find . -name "*.cpp")

# Object rule
%.o : %.cpp
	@echo "Compiling $< ..."
	@$(CC) $(INCLUDE_PATH) $(CXXFLAGS) -o $@ -c $<

# Final target
$(TARGET): $(SRC:.cpp=.o)
	@echo "Linking $@..."
	@$(CC) -o $@ $(SRC:.cpp=.o) $(LIBS)

.PHONY: all clean

all: $(TARGET)

clean:
	\rm -f *.o
	\rm -f $(TARGET)
