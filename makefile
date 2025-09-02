#CC          = /opt/homebrew/bin/g++-15 -std=c++11  -O3 -fopenmp
#CC          = /opt/homebrew/bin/gcc-15 -std=c++11  -O3 -fopenmp
CC          = /opt/homebrew/opt/llvm/bin/clang++ -std=c++11  -O3 -fopenmp


PRGM = DG
DEPS = $(shell find ./src -name "*.hpp")
SRC = $(shell find ./src -name "*.cpp")
OBJ = $(SRC:%.c=%.o) 

all: $(PRGM)

%.o: %.c $(DEPS)
	$(CC) -c $< -o $@

$(PRGM): $(OBJ)
	$(CC) $^ -o $@


