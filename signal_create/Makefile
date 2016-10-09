#######################################################
#
#    FFT Test
#    ===========================================
#    Elouan Autret
#    -------------------------------------------
#
#######################################################


#######################################################
# Extra Build Options
#######################################################

USE_GNUPLOT= 0

#######################################################
# TOOLCHAINS
#######################################################


export TOOLCHAIN = linux-local


#######################################################
# FILES
#######################################################

# Directories
export BUILD_DIR = build
SRC_DIR = ./src
OUTPUT_DIR = ./
INC = -I./src/

# Sources

CORE = Signal/signal.cpp


NODES =

SYSTEM_SERVICES =

SRC_MAIN = main.cpp

SRC = $(CORE)




# Object files
OBJECTS = $(addprefix $(BUILD_DIR)/, $(SRC:.cpp=.o))
OBJECT_MAIN = $(addprefix $(BUILD_DIR)/, $(SRC_MAIN:.cpp=.o))

# Target Output
EXECUTABLE = create_signal.run
OBJECT_FILE = $(BUILD_DIR)/objects.tmp


#######################################################
# TOOLS
#######################################################
MKDIR_P = mkdir -p

CPPFLAGS = -g -Wall -pedantic -Werror -std=c++14


C_TOOLCHAIN = 1
CC = gcc
CXX = g++
SIZE = size


DEFINES = -DGNUPLOT=$(USE_GNUPLOT)

#######################################################
# Rules
#######################################################

.PHONY: clean

all: $(EXECUTABLE) stats

	#  Create the directories needed
$(BUILD_DIR):
	@$(MKDIR_P) $(BUILD_DIR)

# Link and build
$(EXECUTABLE) : $(BUILD_DIR) $(OBJECTS) $(OBJECT_MAIN)
	rm -f $(OBJECT_FILE)
	@echo Linking object files
	@echo -n " " $(OBJECTS) >> $(OBJECT_FILE)
	$(CXX) $(LDFLAGS) $(INC) @$(OBJECT_FILE) $(OBJECT_MAIN) -Wl,-rpath=./ -o $@
	@echo Built using toolchain: $(TOOLCHAIN)

# Compile CPP files into the build folder
$(BUILD_DIR)/%.o:$(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo Compiling CPP File: $@

	@$(CXX) -c $(CPPFLAGS) $(INC) -o ./$@ $< $(DEFINES)



#####################################################################
# Tool Rules

stats:$(EXECUTABLE)
	@echo Final executable size:
	$(SIZE) $(EXECUTABLE)

clean:
	@echo Removing existing object files and executable
	@rm -f -r $(BUILD_DIR)
	@rm -f $(EXECUTABLE)
	@echo DONE
