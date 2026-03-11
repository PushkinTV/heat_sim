CC = gcc
CFLAGS = -Wall -Wextra -O2 -lm

TARGET = heat_sim
SRCS = main.c heat.c
OBJS = $(SRCS:.c=.o)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c heat.h settings.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean
