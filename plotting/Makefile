PROGRAM       = DrawingClasses
version       = JetBackgroundStudy
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -D$(version) # Add -g for debugger
LDFLAGS       = -O
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)

CXXFLAGS  += $(INCS)

# Put to hdrsdict header all files the inherit TObject
#HDRSDICT = SplitCanvas.h
# HDRS += $(HDRSDICT) ... nanoDict.h

HDRS += AlgorithmLibrary.h JetBackgroundCard.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM) :    $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -shared -lEG -L$(PWD) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM).so
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) $(PROGRAM).o $(PROGRAM).so

cl:  clean $(PROGRAM)

# Dictionary is needed for all classes inheriting TObject from root
#nanoDict.cc: $(HDRSDICT)
#		@echo "Generating dictionary ..."
#		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
#		@rootcling nanoDict.cc -c -D$(version) $(HDRSDICT)
