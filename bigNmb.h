#ifndef BIGNMB_H
#define BIGNMB_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h> 

// On définit la capacité max (en mots de 64 bits)
#define T_MAX 128 
// On definie capacite limiter pour que le resultat ne depasse pas T_MAX
#define T_LIMIT 32

// Type pour definie les grand nombre 
// Premier indice est le taille de nombre (blocs considerer)
typedef uint64_t* bignmb;

typedef struct 
{
    bignmb n;
    bignmb e;
    bignmb d;
} cles;


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
extern uint64_t asm_mult(uint64_t* dest, uint64_t* src, uint64_t scalaire, uint64_t len);

// --- Fonctions Utilitaires internes ---
void copy_big(bignmb dest, bignmb src);
void Div2_big(bignmb a);
bignmb inverse(bignmb a, bignmb m);
uint64_t big_mod_small(bignmb n, uint64_t d);
bignmb Sub_mod_big(bignmb a, bignmb b, bignmb n);
// Modulo a % b
bignmb Mod_big(bignmb a, bignmb b);

// --- Operation Arthimetiquer ---
// Addition a = a + b
bignmb Add_big(bignmb a, bignmb b);
// Addition a = a - b
bignmb Sub_big(bignmb a, bignmb b);
// Addition a = a + b*n
bignmb Mult_big(bignmb a, bignmb b);
// Modulo de res = a*b (mod n)
bignmb Montgomery_Mult(bignmb a, bignmb b, bignmb n);
// Fonction de comparaion si a>b reurn 1 sino return 0 et si a==b return 2
int Comp_big(bignmb a, bignmb b);
// Fonction pour pré-calculer la constante de Montgomery
bignmb Calculer_R2_Mod_N(bignmb n);
// Puissance modulaire (prend R2 en paramètre maintenant)
bignmb Puiss_big(bignmb a, bignmb e, bignmb n, bignmb precalc_r2);

// --- Fonction de Tests ---
// Décompose n-1 en 2^s * d
// Retourne s et modifie d sur place
int decomposition_n_minus_1(bignmb n, bignmb d);
// Test de Miller-Rabin (k itérations)
int Miller_Rabin(bignmb n, int k);
// Vérifier est ce qu'il est divisible par un petit nombre premier (< 256)
int est_divisible_par_petits_premiers(bignmb n);


// Fonction principale pour genere un nombre premier 
bignmb Gen_premier();
// Fonctio trouber les cles pubilque et 
cles Gen_cles();

#endif