# Only the Cuba library is required for the numerical integration


# Recola path
RECOLA_PATH = /Users/rhorry/Research/BUILDS
# your CUBA 4.0 path (testing: /usr/local/Cuba-4.0)
# CUBA_PATH = /Users/rhorry/Research/BUILDS/CUBA/

CXX=g++-14
# compiler flags: -Wall -O2
CXXFLAGS = -Wall -O2 -std=c++11

# compiler name: g++ is default
FF = gfortran
#
INCLUDE_PATH =
LIBS_PATH = 

# For external |M|^2 from Recola
ifneq "$(RECOLA_PATH)" ""
INCLUDE_PATH += -I$(RECOLA_PATH)/include
LIBS_PATH += -L$(RECOLA_PATH)/lib -lrecola
endif

# Define code compilation structure
SRCDIR = ./src
BINDIR = ./bin
OBJDIR = ./obj
INCDIR = ./inc

# Define various source files (i.e. the core files)
SOURCES = tools.cc var.cc ME2_Analytic.cc dsigma.cc 
# Append to this those for the various processes

# Combine these
SOURCES += $(SIGMA_SOURCES)

# Object list obtained from the source list
_OBJ = $(SOURCES:.cc=.o)
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

# Add local source/header files 
INCLUDE_PATH += -I$(SRCDIR) -I$(INCDIR)

CXXFLAGS += $(INCLUDE_PATH)
CXXFLAGS += "-Wl,-no_warn_duplicate_libraries"

# Build all binaries
all: $(BINDIR)/dSigma.exe

# Test program for differential cross-section computation
$(BINDIR)/dSigma.exe: $(SRCDIR)/wrapper_dsigma.cc $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $< $(LIBS_PATH)	

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f $(BINDIR)/*.exe
	rm -vf $(OBJDIR)/*.o

.PHONY: clean
