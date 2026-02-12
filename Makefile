
CC = gcc
CFLAGS = -Wall -Wextra -O3 -g

# Nom de l'exécutable
TARGET = test_prog

SRC_C = test.c bigNmb.c
SRC_ASM = operation.s 

OBJ = $(SRC_C:.c=.o) $(SRC_ASM:.s=.o)


all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.s
	$(CC) $(CFLAGS) -c $< -o $@

# Lancer le programme
run: $(TARGET)
	./$(TARGET)

# Nettoyer les fichiers générés
clean:
	rm -f $(OBJ) $(TARGET)
