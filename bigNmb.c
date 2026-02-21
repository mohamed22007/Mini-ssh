#include "bigNmb.h"
#include <string.h> 
#define NMB_TEST 35

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

        // asm_mult retourne la retenue finale du dernier mot
        uint64_t carry = asm_mult(&res[1 + i], &a[1], scalaire, a[0]);

        // Propager la retenue dans les mots suivants
        uint64_t j = 1 + i + a[0];
        while (carry && j <= T_MAX) {
            res[j] += carry;
            carry = (res[j] < carry) ? 1 : 0;
            j++;
        }
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
            if (Miller_Rabin(n, NMB_TEST)) {
                return n; // Trouvé 
            }
        }
    }
}

// Soustraction modulaire : calcule (a - b) mod n
bignmb Sub_mod_big(bignmb a, bignmb b, bignmb n) {
    bignmb res;
    if (Comp_big(a, b) >= 0) {
        res = Sub_big(a, b);
    } else {
        bignmb tmp = Add_big(a, n);
        res = Sub_big(tmp, b);
        free_big(tmp); // On libère le temporaire
    }
    return res;
}

// Divise un BigInt par 2 (décalage bit à bit vers la droite)
void Div2_big(bignmb a) {
    uint64_t carry = 0;
    for (int i = a[0]; i >= 1; i--) {
        uint64_t next_carry = (a[i] & 1) << 63; // On récupère le bit sortant
        a[i] = (a[i] >> 1) | carry;
        carry = next_carry;
    }
    // Nettoyage de la taille si le mot de poids fort devient 0
    if (a[0] > 1 && a[a[0]] == 0) {
        a[0]--;
    }
}

// Modulo robuste : retourne a % b
bignmb Mod_big(bignmb a, bignmb b) {
    bignmb R = new_big(0);
    copy_big(R, a);
    if (Comp_big(R, b) < 0) return R;

    bignmb D = new_big(0);
    copy_big(D, b);

    int a_bits = 0, b_bits = 0;
    for(int i = a[0]; i >= 1; i--) { if(a[i]) { a_bits = (i-1)*64 + 63; uint64_t t=a[i]; while(!(t&0x8000000000000000ULL)) { t<<=1; a_bits--; } break; } }
    for(int i = b[0]; i >= 1; i--) { if(b[i]) { b_bits = (i-1)*64 + 63; uint64_t t=b[i]; while(!(t&0x8000000000000000ULL)) { t<<=1; b_bits--; } break; } }

    int shift = a_bits - b_bits;
    if (shift < 0) { free_big(D); return R; }

    for (int k = 0; k < shift; k++) {
        uint64_t carry = 0;
        for (int i = 1; i <= (int)D[0]; i++) {
            uint64_t nc = (D[i] >> 63);
            D[i] = (D[i] << 1) | carry;
            carry = nc;
        }
        if (carry && D[0] < T_MAX) { D[0]++; D[D[0]] = carry; }
    }

    for (int k = 0; k <= shift; k++) {
        if (Comp_big(R, D) >= 0) {
            bignmb tmp = Sub_big(R, D);
            free_big(R);
            R = tmp;
        }
        Div2_big(D);
    }

    free_big(D);
    return R;
}

// Euclide étendu : retourne x tel que a*x ≡ 1 (mod m)
// Représentation signée : on maintient (old_r, r) et (old_s, s)
// avec le signe stocké séparément (0 = positif, 1 = négatif)
bignmb inverse(bignmb a, bignmb m) {
    // old_r = a, r = m
    bignmb old_r = new_big(0); copy_big(old_r, a);
    bignmb r     = new_big(0); copy_big(r, m);

    // old_s = 1, s = 0  (coefficients de Bézout pour a)
    bignmb old_s = new_big(1);
    bignmb s     = new_big(0);
    int old_s_neg = 0;  // signe de old_s (0=positif, 1=négatif)
    int s_neg     = 0;

    bignmb zero = new_big(0);

    while (Comp_big(r, zero) != 0) {
        // quotient = old_r / r  (on calcule via Mod_big : old_r = q*r + rem)
        // On a besoin de q = (old_r - rem) / r
        // On utilise : rem = Mod_big(old_r, r), puis q*r = old_r - rem

        bignmb rem = Mod_big(old_r, r);

        // Calculer quotient q = (old_r - rem) / r par soustraction répétée
        // Mais c'est trop lent. On calcule q différemment :
        // On fait la division bit à bit
        bignmb numerator = Sub_big(old_r, rem);
        
        // Division exacte : q = numerator / r  (numerator est divisible par r)
        // On fait ça avec un shift + soustraction
        bignmb q = new_big(0);
        if (Comp_big(numerator, zero) != 0) {
            // Calculer nb de bits de numerator et r
            int n_bits = 0, d_bits = 0;
            for (int i = numerator[0]; i >= 1; i--) {
                if (numerator[i]) {
                    n_bits = (i-1)*64 + 63;
                    uint64_t t = numerator[i];
                    while (!(t & 0x8000000000000000ULL)) { t <<= 1; n_bits--; }
                    break;
                }
            }
            for (int i = r[0]; i >= 1; i--) {
                if (r[i]) {
                    d_bits = (i-1)*64 + 63;
                    uint64_t t = r[i];
                    while (!(t & 0x8000000000000000ULL)) { t <<= 1; d_bits--; }
                    break;
                }
            }
            int shift = n_bits - d_bits;

            // D = r << shift
            bignmb D = new_big(0); copy_big(D, r);
            for (int k = 0; k < shift; k++) {
                uint64_t carry = 0;
                for (int i = 1; i <= (int)D[0]; i++) {
                    uint64_t nc = (D[i] >> 63);
                    D[i] = (D[i] << 1) | carry;
                    carry = nc;
                }
                if (carry && D[0] < T_MAX) { D[0]++; D[D[0]] = carry; }
            }

            bignmb N2 = new_big(0); copy_big(N2, numerator);
            for (int k = shift; k >= 0; k--) {
                if (Comp_big(N2, D) >= 0) {
                    bignmb t2 = Sub_big(N2, D); free_big(N2); N2 = t2;
                    // q += 1 << k
                    // On ajoute 2^k à q : trouver le mot et le bit
                    int word = k / 64 + 1;
                    int bit  = k % 64;
                    if (word <= T_MAX) {
                        if (q[0] < (uint64_t)word) q[0] = word;
                        q[word] += (1ULL << bit);
                        // propager retenue
                        uint64_t carry2 = 0;
                        if (q[word] < (1ULL << bit)) carry2 = 1; // overflow
                        // simple carry propagation
                        for (int wi = word; carry2 && wi <= T_MAX; wi++) {
                            q[wi] += carry2;
                            carry2 = (q[wi] == 0) ? 1 : 0;
                        }
                    }
                }
                // D >>= 1
                uint64_t carry = 0;
                for (int i = (int)D[0]; i >= 1; i--) {
                    uint64_t nc = (D[i] & 1) << 63;
                    D[i] = (D[i] >> 1) | carry;
                    carry = nc;
                }
                if (D[0] > 1 && D[D[0]] == 0) D[0]--;
            }
            free_big(D);
            free_big(N2);
        }
        free_big(numerator);

        // Rotation : old_r, r = r, rem
        free_big(old_r);
        old_r = r;
        r = rem;

        // new_s = old_s - q * s
        // avec gestion du signe
        bignmb qs = Mult_big(q, s);
        free_big(q);

        // new_s_neg et new_s = old_s - q*s (signé)
        bignmb new_s;
        int new_s_neg;

        if (old_s_neg == s_neg) {
            // Même signe : soustraction
            if (Comp_big(old_s, qs) >= 0) {
                new_s = Sub_big(old_s, qs);
                new_s_neg = old_s_neg;
            } else {
                new_s = Sub_big(qs, old_s);
                new_s_neg = !old_s_neg;
            }
        } else {
            // Signes opposés : addition
            new_s = Add_big(old_s, qs);
            new_s_neg = old_s_neg;
        }
        free_big(qs);
        free_big(old_s);
        old_s     = s;
        old_s_neg = s_neg;
        s         = new_s;
        s_neg     = new_s_neg;
    }

    // Résultat dans old_s, de signe old_s_neg
    // Si négatif, on ramène dans [0, m) : résultat = m - old_s
    bignmb result;
    if (old_s_neg) {
        result = Sub_big(m, old_s);
    } else {
        result = new_big(0);
        copy_big(result, old_s);
    }

    // Réduction finale mod m au cas où
    while (Comp_big(result, m) >= 0) {
        bignmb t = Sub_big(result, m);
        free_big(result);
        result = t;
    }

    free_big(old_r); free_big(r);
    free_big(old_s); free_big(s);
    free_big(zero);
    return result;
}

// Géneration de cles public et prive 
cles Gen_cles(){
    // container de Resulats
    cles resultat;
    // Cles pubilc e (Valeur Standard)
    bignmb e = new_big(65537);
    resultat.e= e;

    bignmb un = new_big(1);
    // Indentifient d'euler
    bignmb phi_n;
    bignmb n;

    // Pour assurer que e et phi_n sont premier entre eux
    do {
        // Les nombres premiers necessaire pour le codege
        bignmb p = Gen_premier();
        bignmb q = Gen_premier();

        // n et euler de n 
        bignmb n = Mult_big(p,q);
        bignmb phi_n = Mult_big(Sub_big(p,un),Sub_big(q,un));

        //nettoyer
        free(p);
        free(q);

    } while (big_mod_small(phi_n, e) != 0);
    
    resultat.n = n ;
    // Trouver le inverse de e 
    bignmb d = inverse(e,phi_n);
    resultat.d = d;

    // Nettoyage 
    free(phi_n);
    free(un);

    return resultat ;
}