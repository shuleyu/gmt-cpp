# Related directories.
# SACDIR    := /opt/sac-101.6
# SACDIR2   := /opt/sac

# Compile parameters & dirs, some could be overwritten in Run.sh
# Notice: the order of library names in LIBS could matter.
#   If libA.a depends on libB.a, then -lA should appears before -lB.

COMP      := c++ -std=c++14 -Wall # -fPIC
OUTDIR    := .
INCDIR    := -I.
LIBDIR    := -L/usr/local/Cellar/gmt@5/5.4.5_4/lib
LIBS      := -lgmt -lm

# all *cpp files
SRCFILES  := $(wildcard *.cpp)
DEPFILES  := $(patsubst %.cpp, $(OUTDIR)/%.d, $(SRCFILES))
OBJS      := $(patsubst %.d, %.o, $(DEPFILES))

# main files
MAINS     := $(filter-out %.fun.cpp, $(SRCFILES))
EXEFILES  := $(patsubst %.cpp, $(OUTDIR)/%.out, $(MAINS))

# function files
FUNFILES  := $(wildcard *fun.cpp)
FUNOBJS   := $(patsubst %.cpp, $(OUTDIR)/%.o, $(FUNFILES))

all: $(EXEFILES) $(OBJS)
	@echo > /dev/null

# Resolve dependencies automatically.
-include $(DEPFILES)

%.out: %.o $(FUNOBJS)
	@echo "Updating: $@ ..."
	@$(COMP) -o $@ $^ $(INCDIR) $(LIBDIR) $(LIBS)

$(OUTDIR)/%.o: %.cpp
	@echo "Updating: $@ ..."
	@$(COMP) -MD -MP -c $< -o $@ $(INCDIR)

clean:
	rm -f $(OUTDIR)/*.out $(OUTDIR)/*.o $(OUTDIR)/*.d
