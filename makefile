CC          = g++ -std=c++11  -O3

PRGM = DG
DEPS = $(shell find ./src -name "*.hpp")
SRC = $(shell find ./src -name "*.cpp")
OBJ = $(SRC:%.c=%.o) 

all: $(PRGM)

%.o: %.c $(DEPS)
	$(CC) -c $< -o $@

$(PRGM): $(OBJ)
	$(CC) $^ -o $@


