# Final executable
TARGET=dodri

# subdirs
SRCDIR = src
INCDIR = inc
OBJDIR = obj
BINDIR = bin

# compiler rules
CC = g++ #-std=c++11
DEBUG = -g
OPT = -O3 -Wall
IMGMGK= -fno-strict-aliasing -pthread  `GraphicsMagick++-config --cppflags --ldflags --libs`

CFLAGS=  -I$(INCDIR) $(OPT) $(IMGMGK) 

SOURCES := $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/io/*.cpp)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o) 

$(BINDIR)/$(TARGET): $(OBJECTS) $(INCLUDES)
	$(CC) $(OBJECTS) $(CFLAGS) -o $@

all: $(BINDIR)/$(TARGET) 
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCLUDES) 
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: 
clean:
	rm -f $(OBJECTS) $(BINDIR)/$(TARGET)
debug:	CFLAGS = -I$(INCDIR) $(IMGMGK) $(DEBUG) 

debug: all 
