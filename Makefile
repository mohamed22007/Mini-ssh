# Compilateur et options
CC = gcc
CFLAGS = -Wall -Wextra -g
TARGET = test

SRC = test.c bigNmb.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
