# subdirs
SRCDIR = src
INCDIR = inc
OBJDIR = obj
BINDIR = bin

# Final executable
TARGET = $(BINDIR)/dodri

# Using != runs the command immediately, so the actual flags appear
# in the compile and link commands
GM_CPP	!= GraphicsMagick++-config --cppflags
GM_LD	!= GraphicsMagick++-config --ldflags --libs

# compiler rules
# -pthread appears to be unused
# Allow user to override default compiler and optional flags using make
# arguments or environment variables
CXX		?= g++
CXXFLAGS	?= -O3 -Wall -fno-strict-aliasing

# Add required flags
CXXFLAGS	+= -I$(INCDIR) $(GM_CPP)
LDFLAGS		+= $(GM_LD)

DEBUG		= -g

SOURCES		:= $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/io/*.cpp)
INCLUDES	:= $(wildcard $(INCDIR)/*.h)
OBJECTS		:= $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Default install tools, can be overridden by user
MKDIR		?= mkdir
INSTALL		?= install

# Support staged installs.  Most package managers will provide
# DESTDIR and PREFIX.
DESTDIR		?= stagedir
PREFIX		?= /usr/local

all: $(TARGET)

$(BINDIR):
	$(MKDIR) -p $(BINDIR)

$(OBJDIR)/io:
	$(MKDIR) -p $(OBJDIR)/io

$(TARGET): $(BINDIR) $(OBJDIR)/io $(OBJECTS) $(INCLUDES)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY:
clean:
	rm -f $(OBJECTS) $(TARGET)

debug: CXXFLAGS = -I$(INCDIR) $(IMGMGK) $(DEBUG)

debug: all

install: all
	$(MKDIR) -p $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -c $(TARGET) $(DESTDIR)$(PREFIX)/bin

