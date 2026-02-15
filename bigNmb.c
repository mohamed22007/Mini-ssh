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
    for (uint64_t i = n[0]; i >= 1; i--) {
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
    if (!res) exit(1);
    
    res[0] = a[0] + b[0]; 
    if (res[0] > T_MAX) res[0] = T_MAX; 

    for (uint64_t i = 0; i < b[0]; i++) {
        uint64_t scalaire = b[i+1]; 
        if (scalaire == 0) continue;

        // On appelle asm_mul
        asm_mult(&res[1 + i], &a[1], scalaire, a[0]);
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

void shift_left_R(bignmb res, bignmb a, int shift_amount) {
    int old_len = a[0];
    for (int i = old_len; i >= 1; i--) {
        res[i + shift_amount] = a[i];
    }

    // On remplit le bas avec des zéros
    for (int i = 1; i <= shift_amount; i++) {
        res[i] = 0;
    }

    // Mise à jour de la taille
    res[0] = old_len + shift_amount;
}


void shift_right_R(bignmb res, bignmb a, int shift_amount) {
    int old_len = a[0];

    if (shift_amount >= old_len) {
        res[0] = 1;
        res[1] = 0;
        return;
    }

    // On copie la partie haute vers le bas
    int new_len = old_len - shift_amount;
    
    for (int i = 1; i <= new_len; i++) {
        res[i] = a[i + shift_amount];
    }

    // Mise à jour taille
    res[0] = new_len;
}

// Calcule l'inverse de n0 modulo 2^64
uint64_t calcul_n_prime(uint64_t n0) {
    uint64_t x = n0;
    x = x * (2 - n0 * x);
    x = x * (2 - n0 * x);
    x = x * (2 - n0 * x);
    x = x * (2 - n0 * x);
    x = x * (2 - n0 * x);
    return -x; 
}

int Comp_big(bignmb a, bignmb b) {
    if (a[0] > b[0]) return 1;  // a est plus long
    if (a[0] < b[0]) return -1; // b est plus long

    // on compare de haut en bas
    for (int i = a[0]; i >= 1; i--) {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0; // Ils sont strictement égaux
}

bignmb Montgomery_Mult(bignmb a, bignmb b, bignmb n) {
    uint64_t k = n[0]; 
    uint64_t n_inv = calcul_n_prime(n[1]); 

    // Allocation sur la pile 
    uint64_t T[T_MAX * 2 + 2];
    memset(T, 0, sizeof(T));

    // Multiplication a * b stockée dans T (Integrated Multiplication)
    for (int i = 0; i < (int)a[0]; i++) {
        uint64_t u = a[i+1];
        uint64_t carry = 0;
        for (int j = 0; j < (int)b[0]; j++) {
            // Utilisation de __int128 pour multiplication 64x64 -> 128
            unsigned __int128 prod = (unsigned __int128)u * b[j+1] + T[i+j] + carry;
            T[i+j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        T[i + b[0]] += carry;
    }

    // Boucle de réduction de Montgomery
    for (int i = 0; i < (int)k; i++) {
        // m = T[i] * n_inv (mod 2^64)
        uint64_t m_val = T[i] * n_inv;
        
        // T = T + m * n
        uint64_t carry = 0;
        for (int j = 0; j < (int)k; j++) {
            unsigned __int128 prod = (unsigned __int128)m_val * n[j+1] + T[i+j] + carry;
            T[i+j] = (uint64_t)prod; 
            carry = (uint64_t)(prod >> 64);
        }
        
        // Propagation retenue
        unsigned __int128 sum = (unsigned __int128)T[i+k] + carry;
        T[i+k] = (uint64_t)sum;
        T[i+k+1] += (uint64_t)(sum >> 64);
    }

    // Réduction finale et copie
    bignmb res = new_big(0);
    res[0] = k;
    
    for(int i=0; i < (int)k; i++) {
        if (k+i < T_MAX*2 + 2) res[i+1] = T[k+i];
    }
    
    // Nettoyage taille
    while (res[0] > 1 && res[res[0]] == 0) res[0]--;

    // Si Res >= N, on soustrait N
    if (Comp_big(res, n) >= 0) {
        bignmb tmp = Sub_big(res, n);
        free_big(res);
        return tmp;
    } 
    return res;
}

bignmb Calculer_R2_Mod_N(bignmb n) {
    uint64_t k = n[0];
    
    // Allocation sur la pile 
    uint64_t R[T_MAX + 2]; 
    memset(R, 0, sizeof(R));
    
    // Initialisation à 1
    R[0] = 1;
    R[1] = 1;

    // R^2 = 2^(128 * k)
    uint64_t total_bits = 128 * k;

    for (uint64_t i = 0; i < total_bits; i++) {
        // Doublement (Shift Left de 1 bit)
        uint64_t carry = 0;
        for (int j = 1; j <= (int)R[0]; j++) {
            uint64_t next_carry = (R[j] >> 63);
            R[j] = (R[j] << 1) | carry;
            carry = next_carry;
        }
        if (carry) {
            R[0]++;
            R[R[0]] = carry;
        }

        // Modulo : Si R >= n, R = R - n
        // Comparaison
        int plus_grand = 0;
        if (R[0] > n[0]) plus_grand = 1;
        else if (R[0] == n[0]) {
            for (int j = n[0]; j >= 1; j--) {
                if (R[j] > n[j]) { plus_grand = 1; break; }
                if (R[j] < n[j]) { plus_grand = 0; break; }
            }
        }

        // Soustraction in-place
        if (plus_grand) {
            uint64_t borrow = 0;
            for (int j = 1; j <= (int)n[0]; j++) {
                uint64_t sub = n[j];
                uint64_t r_val = R[j];
                uint64_t diff = r_val - sub - borrow;
                
                if (r_val < sub + borrow) borrow = 1;
                else borrow = 0;
                // Cas special overflow sub+borrow
                if (sub == 0xFFFFFFFFFFFFFFFF && borrow == 1) borrow = 1;

                R[j] = diff;
            }
            // Propager borrow
            int idx = n[0] + 1;
            while (borrow && idx <= (int)R[0]) {
                if (R[idx] == 0) { R[idx] = 0xFFFFFFFFFFFFFFFF; borrow=1; }
                else { R[idx]--; borrow=0; }
                idx++;
            }
            while (R[0] > 1 && R[R[0]] == 0) R[0]--;
        }
    }

    bignmb res_final = new_big(0);
    res_final[0] = R[0];
    for(int i=1; i<=(int)R[0]; i++) res_final[i] = R[i];
    return res_final;
}

// a^e mod n (Optimisé: prend R2 en argument)
bignmb Puiss_big(bignmb a, bignmb e, bignmb n, bignmb R2){
    
    // Passage dans le domaine de Montgomery
    // a_barre = a * R mod n
    bignmb a_barre = Montgomery_Mult(a, R2, n);
    
    // x_barre = 1 * R mod n
    bignmb un = new_big(1);
    bignmb x_barre = Montgomery_Mult(un, R2, n);
    free_big(un);
    
    // Exponentiation binaire
    for (int i = e[0]; i >= 1; i--) {
        uint64_t w = e[i];
        int start_bit = (i == e[0]) ? 63 : 63;
        
        // Skip les zeros du debut
        if (i == e[0]) {
             while (start_bit >= 0 && !((w >> start_bit) & 1)) start_bit--;
        }

        for (int j = start_bit; j >= 0; j--) {
            // Square
            bignmb sq = Montgomery_Mult(x_barre, x_barre, n);
            free_big(x_barre);
            x_barre = sq;
            
            // Multiply
            if ((w >> j) & 1) {
                bignmb mul = Montgomery_Mult(x_barre, a_barre, n);
                free_big(x_barre);
                x_barre = mul;
            }
        }
    }
    
    // Sortie du domaine
    bignmb un_final = new_big(1);
    bignmb res = Montgomery_Mult(x_barre, un_final, n);
    
    free_big(un_final);
    free_big(a_barre);
    free_big(x_barre);

    return res;
}


// Copie profonde de bignmb (utile pour Miller Rabin)
void copy_big(bignmb dest, bignmb src) {
    memcpy(dest, src, (src[0] + 1) * sizeof(uint64_t));
}

// Fonction interne pour décomposer n-1 = 2^s * d
int decomposition_n_minus_1(bignmb n, bignmb d) {
    // d = n - 1
    copy_big(d, n);
    
    // Soustraction de 1 
    if (d[1] >= 1) {
        d[1] -= 1;
    } else {
        int i = 1;
        while (d[i] == 0) {
            d[i] = 0xFFFFFFFFFFFFFFFF;
            i++;
        }
        d[i]--;
    }

    int s = 0;
    
    // On vérifie le bit de poids faible du premier mot
    while ((d[1] & 1) == 0) {
        s++;
        
        uint64_t carry = 0;
        for (int i = d[0]; i >= 1; i--) {
            uint64_t next_carry = (d[i] & 1) << 63; 
            d[i] = (d[i] >> 1) | carry;
            carry = next_carry;
        }
        
        if (d[0] > 1 && d[d[0]] == 0) d[0]--;
        
        if (d[0] == 1 && d[1] == 0) break; 
    }
    return s;
}

int Miller_Rabin(bignmb n, int k) {

    // Cas triviaux
    if (n[0] == 1 && n[1] <= 3) return 1;
    if ((n[1] & 1) == 0) return 0; 

    // Décomposition n-1 = 2^s * d
    bignmb d = new_big(0);
    int s = decomposition_n_minus_1(n, d);
    
    // Pré-calculs constants
    bignmb R2 = Calculer_R2_Mod_N(n);
    bignmb n_minus_1 = Sub_big(n, new_big(1)); // n-1 pour comparaisons
    bignmb one = new_big(1);
    bignmb two = new_big(2); 

    int est_premier = 1; 

    for (int i = 0; i < k; i++) {
        // a = random entre [2, n-2]
        bignmb a = gen_aleatoire();

        // Sécurite
        if (Comp_big(a, n) >= 0) {
            a[0] = n[0]; // On force la taille
            if (a[a[0]] >= n[n[0]]) a[a[0]] = n[n[0]] - 1;
        }
        if (a[0] == 1 && a[1] < 2) a[1] = 2;

        // x = a^d mod n
        bignmb x = Puiss_big(a, d, n, R2);
        free_big(a);

        // Si x == 1 ou x == n-1, ce témoin ne prouve rien, on continue
        if (Comp_big(x, one) == 0 || Comp_big(x, n_minus_1) == 0) {
            free_big(x);
            continue;
        }

       
        int temoin_compose = 1; 
        for (int r = 1; r < s; r++) {
            // x = x^2 mod n
            bignmb x_prev = x;
            x = Puiss_big(x, two, n, R2);
            free_big(x_prev);

            if (Comp_big(x, n_minus_1) == 0) {
                temoin_compose = 0; 
                break;
            }
            
            if (Comp_big(x, one) == 0) {
                free_big(x);
                est_premier = 0;
                goto fin;
            }
        }

        free_big(x);
        if (temoin_compose) {
            est_premier = 0;
            goto fin;
        }
    }

fin:
    free_big(d);
    free_big(R2);
    free_big(n_minus_1);
    free_big(one);
    free_big(two);
    
    return est_premier;
}


// Liste des nombres premiers < 256 pour le crible rapide
static const uint16_t petits_premiers[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 
    211, 223, 227, 229, 233, 239, 241, 251
};
#define NB_PETITS_PREMIERS 53 

// Fonction utilitaire : (BigInt % uint64_t)
uint64_t big_mod_small(bignmb n, uint64_t d) {
    uint64_t remainder = 0;
    for (int i = n[0]; i >= 1; i--) {
        unsigned __int128 temp = ((unsigned __int128)remainder << 64) | n[i];
        remainder = (uint64_t)(temp % d);
    }
    return remainder;
}

// Vérification avec petit nombre
int est_divisible_par_petits_premiers(bignmb n) {
    for (int i = 0; i < NB_PETITS_PREMIERS; i++) {
        if (big_mod_small(n, petits_premiers[i]) == 0) return 1;
    }
    return 0;
}

// Fonction principale : Génère un grand nombre premier
bignmb Gen_premier() {
    bignmb n = NULL;
    
    while (1) {
        if (n) free_big(n);
        
        // Génération aléatoire
        n = gen_aleatoire();
        
        n[1] |= 1ULL; // Impair
        
        // Verification avec petit nombre
        if (est_divisible_par_petits_premiers(n)) {
            continue; 
        }

        // Miller-Rabin (Test probabiliste)
        if (Miller_Rabin(n, 5)) {
            // Pour La sécurite
            if (Miller_Rabin(n, 20)) {
                return n; // Trouvé 
            }
        }
    }
}