#define _POSIX_C_SOURCE 199309L // <--- INDISPENSABLE pour CLOCK_MONOTONIC
#include <stdio.h>
#include <stdlib.h> // Nécessaire pour free()
#include <time.h>
#include "bigNmb.h"

// Fonction utilitaire pour mesurer le temps
double get_time_sec() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts); // Maintenant c'est défini !
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

int main() {
    // ... le reste de votre code est parfait ...
    // (J'ai recopié la suite pour être sûr que vous avez le contexte)
    
    printf("=== Test de Performance BigInt (T_MAX = %d mots / %d bits) ===\n", T_MAX, T_MAX*64);

    // 1. Génération
    printf("\n[1] Génération de nombres aléatoires...\n");
    bignmb a = gen_aleatoire(T_MAX); 
    bignmb b = gen_aleatoire(T_MAX);

    print_big(a); // Petite correction: print_big n'a pas besoin de "A=" avant si vous l'avez mis dedans, sinon gardez votre printf
    print_big(b);

    double start, end;
    bignmb res;

    // 2. Test Addition
    printf("\n[2] Test Addition (1 million d'itérations)...\n");
    start = get_time_sec();
    for(int i = 0; i < 1000000; i++) {
        res = Add_big(a, b);
        free_big(res); 
    }
    end = get_time_sec();
    printf("-> Temps Total : %.4f sec\n", end - start);
    printf("-> Vitesse     : %.2f ns/op\n", ((end - start) / 1000000) * 1e9);


    // 3. Test Soustraction
    if (a[0] < b[0]) { bignmb tmp = a; a = b; b = tmp; } 
    
    printf("\n[3] Test Soustraction (1 million d'itérations)...\n");
    start = get_time_sec();
    for(int i = 0; i < 1000000; i++) {
        res = Sub_big(a, b);
        free_big(res);
    }
    end = get_time_sec();
    printf("-> Temps Total : %.4f sec\n", end - start);
    printf("-> Vitesse     : %.2f ns/op\n", ((end - start) / 1000000) * 1e9);


    // 4. Test Multiplication
    printf("\n[4] Test Multiplication (10 000 itérations)...\n");
    start = get_time_sec();
    for(int i = 0; i < 10000; i++) {
        res = Mult_big(a, b);
        free_big(res);
    }
    end = get_time_sec();
    printf("-> Temps Total : %.4f sec\n", end - start);
    printf("-> Vitesse     : %.2f us/op\n", ((end - start) / 10000) * 1e6);

    free_big(a);
    free_big(b);
    
    printf("\n=== Fin des Tests ===\n");
    return 0;
}