#ifndef BIG_NUM
#define BIG_NUM

#define T_MAX 70

#include<stdint.h>

// Type pour definir les grands nombre de 1024 bits et moins
typedef struct
{
    int Length; //Longeur actuel de nombre
    uint32_t num[T_MAX]; // N = sommes des (2^32)^(i)*num[i]
}bignmb;

// From int to bignmb
bignmb Int_big(uint32_t a);

// Additon entre deux grands nombres
bignmb Add_big(bignmb a, bignmb b);

// Multiplication entre deux grands nombres 
bignmb Mult_big(bignmb a, bignmb b);

// Soustraction entre deux grand Nombres
bignmb Sous_big(bignmb a, bignmb b);

// Modulo a = b mod n
bignmb Modul_big(bignmb a, bignmb n);

// Puissance modulo b = a^(n) mod m
bignmb Puis_big(bignmb a, bignmb n, bignmb m);

// Generer un big nombre aleatoire 
bignmb Ale_big();

// Afficher un nombre big
void Aff_big(bignmb a);

// Comparer deux grand nombres
int Comp_big(bignmb a, bignmb b);


#endif