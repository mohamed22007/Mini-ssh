#ifndef BIGNMB_H
#define BIGNMB_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// On définit la capacité max (en mots de 64 bits)
#define T_MAX 32 
// On definie capacite limiter pour que le resultat ne depasse 
// pas T_MAX
#define T_LIMIT 16

// Type pour definie les grand nombre 
// Premier indice est le taille de nombre (blocs considerer)
typedef uint64_t* bignmb;

// --- Fonctions C ---
// Initialise un Nombre 
bignmb new_big(uint64_t val);
// Supprimer 
void free_big(bignmb n);
// Afficher un Gran nombre
void print_big(bignmb n);
// Genere un nombre aleatoire 
bignmb gen_aleatoire();


// --- Fonctions Assenmbleur ---
// Addiction deux nombres
extern uint64_t asm_add(uint64_t* dest, uint64_t* src, uint64_t len);

// Soustrait b de a. avec a > b
extern uint64_t asm_sub(uint64_t* a, uint64_t* b, uint64_t len);

// multiplie le Grand Nombre 'src' par un mot de 64 bits 'scalaire'.
extern void asm_mult(uint64_t* dest, uint64_t* src, uint64_t scalaire, uint64_t len);

// --- Operation Arthimetiquer ---
// Addition a = a + b
bignmb Add_big(bignmb a, bignmb b);
// Addition a = a - b
bignmb Sub_big(bignmb a, bignmb b);
// Addition a = a + b*n
bignmb Mult_big(bignmb a, bignmb b);


#endif