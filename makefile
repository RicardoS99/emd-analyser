ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

CC := g++
INC := -I src/
CPPFLAGS := -g -m64 -std=c++11 -pedantic -Wall $(INC) $(ROOTINC)
LIBFLAGS := -shared
VPATH = src/
LIB := lib/

SRCS := $(notdir $(wildcard src/*.C)) \
	$(notdir $(wildcard src/*.cpp))
SRCS1 := $(filter %.cpp,  $(SRCS))
SRCS2 := $(filter %.C,  $(SRCS))
OBJSM := $(SRCS1:.cpp=.o)
OBJSM += $(SRCS2:.C=.o)
OBJS = $(addprefix bin/, $(OBJSM))
OLIB = $(addprefix bin/, $(OBJSM))


all:

SRC: $(OBJS)

bin/%.exe: bin/%.o $(OBJS)
	@echo "making executable $<...[$@]"
	@$(CC) $(CPPFLAGS) -o $@ $^ $(INC) $(ROOTINC) $(ROOTLIB)
	@echo "[$@] created"
	@echo

bin/%.o: %.C
	@echo "compiling $<...[$@]"
	@$(CC) $(CPPFLAGS) -c -o $@ $^
	@echo "[$@] compiled."
	@echo


bin/%.o: %.cpp
	@echo "compiling $<...[$@]"
	@$(CC) $(CPPFLAGS) -c -o $@ $^
	@echo "[$@] compiled."
	@echo

clean:
	@ rm -fv bin/*.o bin/*.exe log/*.txt bin/*.eps

cleandata:
	@ rm -fv *.pdf
