CCPP = clang++

SSE_FLAG = 

# Agregamos la ruta de libomp a CFLAGS
CFLAGS = -std=c++17 -w -O3 -Xpreprocessor -fopenmp -Iinclude -I/opt/homebrew/opt/libomp/include
CFLAGS += $(shell pkg-config --cflags opencv4)

# Agregamos libomp al LDFLAGS
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp $(shell pkg-config --libs opencv4)

SOURCES_CPP := $(shell find . -name '*.cpp')
OBJ := $(SOURCES_CPP:%.cpp=%.o)
HEADERS := $(shell find . -name '.h' -o -name '.hpp')

all: CPM

.cpp.o:
	$(CCPP) -c $(CFLAGS) -o $@ $<

CPM: $(OBJ)
	$(CCPP) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ) CPM
