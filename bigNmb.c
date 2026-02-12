#include "bigNmb.h"
#include <string.h> 

// Initialisation d'un BigInt
bignmb new_big(uint64_t val) {
    // Allocation de T_MAX + 1 (Index 0 pour la taille)
    bignmb res = (uint64_t*)calloc(T_MAX + 1, sizeof(uint64_t));
    if (!res) { perror("Alloc failed"); exit(1); }
    
    res[0] = 1;      
    res[1] = val;    
    return res;
}

void free_big(bignmb n){
    if (n) free(n);
}

void print_big(bignmb n) {
    if (!n) return;
    printf("0x");
    // On affiche du plus poids fort au plus faible
    for (int i = n[0]; i >= 1; i--) {
        if (i == n[0]) printf("%lx", n[i]); 
        else printf("%016lx", n[i]);        
    }
    printf("\n");
}

// Génération aléatoire depuis /dev/urandom
bignmb gen_aleatoire() {
    // On utilise calloc pour que tout le reste soit à 0
    bignmb res = (uint64_t*)calloc(T_MAX + 1, sizeof(uint64_t));
    
    if (!res) {
        fprintf(stderr, "Erreur d'allocation mémoire\n");
        exit(1);
    }

    // On fixe la taille utilisée
    res[0] = T_LIMIT;

    // Lecture depuis /dev/urandom
    FILE *fp = fopen("/dev/urandom", "rb");
    if (!fp) {
        perror("Erreur ouverture /dev/urandom");
        free(res);
        exit(1);
    }
    
    // On lit exactement LIMITE_BLOCS 
    size_t lus = fread(&res[1], sizeof(uint64_t), T_LIMIT, fp);
    fclose(fp);

    if (lus != T_LIMIT) {
        fprintf(stderr, "Erreur de lecture aléatoire\n");
        free(res);
        exit(1);
    }

    // OPTIMISATIONS POUR MILLER-RABIN 
    res[T_LIMIT] |= (1ULL << 63);
    res[1] |= 1ULL;

    return res;
}


bignmb Add_big(bignmb a, bignmb b) {
    // Identifier Grand et Petit pour la taille
    bignmb grand = (a[0] >= b[0]) ? a : b;
    bignmb petit = (a[0] < b[0]) ? a : b;

    // Créer le résultat 
    bignmb res = new_big(0); 
    memcpy(res, grand, (grand[0] + 1) * sizeof(uint64_t));

    // Appel ASM : res += petit
    uint64_t carry = asm_add(&res[1], &petit[1], petit[0]);
    uint64_t i = petit[0] + 1;
    
    while (carry && i <= res[0]) {
        res[i]++;
        carry = (res[i] == 0) ? 1 : 0;
        i++;
    }

    if (carry) {
        if (res[0] < T_MAX) {
            res[0]++;
            res[res[0]] = 1;
        }
    }

    return res;
}

bignmb Sub_big(bignmb a, bignmb b) {
    
    // Copier a dans resultat
    bignmb res = new_big(0);
    memcpy(res, a, (a[0] + 1) * sizeof(uint64_t));

    // Appel ASM : res -= b
    uint64_t borrow = asm_sub(&res[1], &b[1], b[0]);

    uint64_t i = b[0] + 1;
    while (borrow && i <= res[0]) {
        uint64_t temp = res[i];
        res[i]--; 
        borrow = (temp == 0) ? 1 : 0;
        i++;
    }

    while (res[0] > 1 && res[res[0]] == 0) {
        res[0]--;
    }

    return res;
}

bignmb Mult_big(bignmb a, bignmb b) {

    bignmb res = (uint64_t*)calloc(T_MAX * 2 + 1, sizeof(uint64_t));
    res[0] = a[0] + b[0] - 1; // Estimation minimale
    if (res[0] > T_MAX) res[0] = T_MAX; // Sécurité

    for (int i = 0; i < b[0]; i++) {
        uint64_t scalaire = b[i+1]; 
        if (scalaire == 0) continue;

        // On appelle asm_mul
        uint64_t carry = asm_mul(&res[1 + i], &a[1], scalaire, a[0]);
        
        // Ajouter la retenue finale de la multiplication au mot suivant
        res[1 + i + a[0]] += carry; 
    }

    // Recalculer la taille réelle
    int max_len = a[0] + b[0] + 1; 
    if (max_len > T_MAX) max_len = T_MAX;
    res[0] = max_len;

    // Supprimer les zéros de tête
    while (res[0] > 1 && res[res[0]] == 0) {
        res[0]--;
    }

    return res;
}