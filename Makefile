############################################################
#     SFG Simulator — Full Pipeline Makefile
############################################################

CXX      = g++
CXXFLAGS = -std=c++17 -O3 -g -fopenmp

# ----------------------------------------------------------
# Include / Library paths
# ----------------------------------------------------------
INCLUDES = -Iinclude \
           -I$(OPENBLAS_ROOT)/include \
           -I$(HDF5_ROOT)/include

LIBS = -L$(OPENBLAS_ROOT)/lib \
       -L$(HDF5_ROOT)/lib \
       -lopenblas \
       -lhdf5_cpp -lhdf5

# ----------------------------------------------------------
# Automatically collect all .cpp files
# ----------------------------------------------------------
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

TARGET = sfg_simulator

# ----------------------------------------------------------
# Build rules
# ----------------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LIBS) -fopenmp

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# ----------------------------------------------------------
# Clean up
# ----------------------------------------------------------
clean:
	rm -f src/*.o $(TARGET)
