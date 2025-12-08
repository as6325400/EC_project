CXX := g++-14
CXXFLAGS := -std=c++17 -O2 -pipe -Wall -Wextra -Wshadow -Wconversion

SRC_DIR := src
VALIDATE_DIR := validate
TESTCASE_DIR := testcases
BIN_DIR := bin

MAIN_SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
VALIDATE_SOURCES := $(wildcard $(VALIDATE_DIR)/*.cpp)
VALIDATE_BINARIES := $(patsubst $(VALIDATE_DIR)/%.cpp,$(BIN_DIR)/%,$(VALIDATE_SOURCES))

.PHONY: all clean judge

all: $(BIN_DIR)/main $(VALIDATE_BINARIES)

$(BIN_DIR):
	@mkdir -p $@

$(BIN_DIR)/main: $(MAIN_SOURCES) | $(BIN_DIR)
	@if [ -z "$(MAIN_SOURCES)" ]; then \
		echo "No source files found in $(SRC_DIR). Add your solution source files first."; \
		exit 1; \
	fi
	$(CXX) $(CXXFLAGS) $(MAIN_SOURCES) -o $@

$(BIN_DIR)/%: $(VALIDATE_DIR)/%.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

judge: $(BIN_DIR)/checker
	./scripts/judge.sh $(BIN_DIR)/main

clean:
	rm -rf $(BIN_DIR)
