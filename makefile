CXX := g++
MOC := moc

SRC_DIR := src
BUILD_DIR := build
BIN_DIR := bin
TARGET := $(BIN_DIR)/INCO
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
HEADERS := $(wildcard $(SRC_DIR)/*.h)

CXXFLAGS := -O3 -g -std=c++23 -Wall -fopenmp

.PHONY : all, clean

all : $(TARGET) $(HEADERS) $(SOURCES)

$(TARGET) : $(SOURCES) $(HEADERS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SOURCES) $(HEADERS) -o $(TARGET)

run : $(TARGET)
	./$(TARGET)
clean :
	rm -f INCO.exe

