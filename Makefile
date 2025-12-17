NAME := surf
CXX := g++
CC := gcc
LD := g++

DEFINES := -DUSE_BARE_C_BESSEL

LDFLAGS = -lopenblas -ldivERGe -fopenmp
CXXFLAGS = $(DEFINES) $(INCLUDES) $(FLAGS) -std=c++17 -Wall -Wextra -pedantic -fopenmp -MMD
CFLAGS = $(DEFINES) $(INCLUDES) $(FLAGS) -std=c11 -Wall -Wextra -pedantic -fopenmp -MMD

-include Makefile.local

SRC_C := $(wildcard *.c interactpolate/*.c)
OBJ_C := $(patsubst %.c,%.c.o,$(SRC_C))
SRC_CPP := $(wildcard *.cpp)
OBJ_CPP := $(patsubst %.cpp,%.cpp.o,$(SRC_CPP))

.PHONY: all clean

all: $(NAME)

-include *.d

$(NAME): $(OBJ_CPP) $(OBJ_C)
	$(LD) $^ -o $@ $(LDFLAGS)

%.cpp.o: %.cpp %.hpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.cpp.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.c.o: %.c %.h
	$(CC) -c $< -o $@ $(CFLAGS)

%.c.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	-@$(RM) -f *.cpp.o *.cpp.d *.c.o *.c.d interactpolate/*.c.o interactpolate/*.c.d $(NAME)
