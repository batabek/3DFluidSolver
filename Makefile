CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

TARGET := FluidSolver
SRCS   := Main.cpp Fluid3D.cpp
OBJS   := $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp Fluid3D.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all run clean
