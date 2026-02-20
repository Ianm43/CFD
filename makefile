CXX := g++
MOC := moc
CXXFLAGS := -O3 -g -std=c++23 -Wall -fopenmp


SRC_DIR := src
BUILD_DIR := build
BIN_DIR := bin
TARGET := $(BIN_DIR)/INCO
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
HEADERS := $(wildcard *.h)


.PHONY : all, clean

all : $(TARGET) 

$(TARGET) : $(SOURCES) $(HEADERS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SOURCES) $(HEADERS) -o $(TARGET)

run : $(TARGET)
	./$(TARGET)
clean :
	rm -f INCO.exe

