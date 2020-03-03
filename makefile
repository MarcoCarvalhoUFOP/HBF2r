Gpp = gcc -O3 -Wno-unused-result
SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.d)
OBJSDIR = obj
CC = gcc

OBJS = $(patsubst %, $(OBJSDIR)/%, $(SRCS:.c=.o) )
DEPS = $(patsubst %, $(OBJSDIR)/%, $(SRCS:.c=.d) )


all: hbf2r

hbf2r: $(OBJS)
		$(Gpp) $^ -o $@

$(OBJSDIR)/%.o: %.c
		$(Gpp) -MMD -MP -c $< -o $@

.PHONY: clean

clean:
		$(RM) $(OBJS) $(DEPS)

-include $(DEPS)
