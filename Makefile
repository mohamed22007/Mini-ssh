CC = gcc
CFLAGS = -Wall -Wextra -O3 -g

# Nom de l'exécutable final
TARGET = Mssh

# Fichiers sources
SRC_C = main.c rsa.c bigNmb.c
SRC_ASM = operation.s 

# Fichiers objets générés
OBJ = $(SRC_C:.c=.o) $(SRC_ASM:.s=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

# Règle de compilation pour les fichiers C
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Règle de compilation pour l'Assembleur
%.o: %.s
	$(CC) $(CFLAGS) -c $< -o $@

# Raccourci pour tester le programme avec keygen (exemple)
run: $(TARGET)
	./$(TARGET) keygen id_rsa.pub id_rsa

# Nettoyer les fichiers générés
clean:
	rm -f $(OBJ) $(TARGET)