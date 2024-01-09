CXX      = g++

CXXFLAGS= $(shell root-config --cflags)
LIBS    = $(shell root-config --libs ) 

SOURCES = NtupleVarsTProxy.cc AnalyzeTProxytBSM.cc
HEADERS = NtupleVarsTProxy.h AnalyzeTProxytBSM.h 
#LinkDef.h
OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = analyzeTProxytBSM

all: $(SOURCES) $(EXECUTABLE)

%.o: %.cc $(HEADERS)
	@echo Compiling $<...
	$(CXX) $(CXXFLAGS) -c -o  $@  $< 

$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $(EXECUTABLE) ..."
	@echo "@$(CXX) -g $(LIBS) $(OBJECTS)+ -o $@"
	@$(CXX) -o  $@ $^ $(LIBS) -lTMVA
	@echo "done"

# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:  $(OBJECTS)

# The default target, which gives instructions, can be called regardless of whether or not files need to be updated.
.PHONY : clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

###
Ntuplevariables.o: NtupleVarsTProxy.h 
AnalyzeLightBSM.o: NtupleVarsTProxy.h AnalyzeTProxytBSM.h 
#Ntuplevariables.o: NtupleVarsTProxy.h LinkDef.h
#AnalyzeLightBSM.o: NtupleVarsTProxy.h AnalyzeTProxytBSM.h LinkDef.h
