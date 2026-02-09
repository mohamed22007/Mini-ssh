#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "bigNmb.h"

bignmb Int_big(uint32_t a){
    bignmb resultat;
    for(int i=0; i<T_MAX; i++) resultat.num[i] = 0;
    
    resultat.num[0] = a;
    resultat.Length = 1; // Longueur minimum de 1
    return resultat;
}

bignmb Add_big(bignmb a, bignmb b){
    bignmb resultat;
    // Initialisation propre
    for(int i=0; i<T_MAX; i++) resultat.num[i] = 0;

    // Pointeurs pour identifier petit et grand 
    bignmb *grand = (a.Length >= b.Length) ? &a : &b;
    bignmb *petit = (a.Length >= b.Length) ? &b : &a;

    uint32_t reste = 0;
    int i;

    for(i = 0; i < petit->Length; i++){
        uint64_t res = (uint64_t)grand->num[i] + petit->num[i] + reste;
        
        resultat.num[i] = (uint32_t)res;     // Stocke les 32 bits du bas
        reste = (uint32_t)(res >> 32);       // Récupère la retenue 
    }

    // Continuer avec le reste du grand nombre
    for(; i < grand->Length; i++){
        uint64_t res = (uint64_t)grand->num[i] + reste;
        
        resultat.num[i] = (uint32_t)res;
        reste = (uint32_t)(res >> 32); 
    }

    // ÉTAPE 3 : S'il reste une retenue à la toute fin
    if (reste > 0){
        resultat.num[i] = reste;
        resultat.Length = grand->Length + 1;
    } else {
        resultat.Length = grand->Length;
    }


    return resultat;
}

// Verifier que a >= b
int Comp_big (bignmb a, bignmb b){
    if (a.Length > b.Length) return 1;
    if (a.Length < b.Length) return 0;

    // Si les tailles sont égales, on compare bloc par bloc en descendant
    for (int i = a.Length - 1; i >= 0; i--) {
        if (a.num[i] > b.num[i]) return 1;
        if (a.num[i] < b.num[i]) return 0;
    }
    
    return 2; // Donc a = b est vrai
}

// Soustraction 
bignmb Sous_big(bignmb a , bignmb b){

    int comp = Comp_big(a,b);
    if (comp == 0){
        printf("Soustraction de a plus petit que b");
        exit(1) ;
    }

    // Initialiser le resultat
    bignmb resultat ;
    for(int j = 0 ; j < T_MAX; j++) resultat.num[j] = 0 ;
    resultat.Length = a.Length ;

    // Si il sont egaux 
    if (comp == 2){
        resultat.Length = 0 ;
        return resultat ;
    }

    int32_t borrow = 0; // Emprunt (0 ou 1)
    int i;

    // Soustraction sur la longueur de B
    for (i = 0; i < b.Length; i++) {
        // On prend A, on enlève B, on enlève l'emprunt précédent
        int64_t diff = (int64_t)a.num[i] - b.num[i] - borrow;

        if (diff < 0) {
            resultat.num[i] = (uint32_t)(diff + (1LL << 32));
            borrow = 1;
        } else {
            resultat.num[i] = (uint32_t)diff;
            borrow = 0;
        }
    }

    for (; i < a.Length; i++) {
        int64_t diff = (int64_t)a.num[i] - borrow;

        if (diff < 0) {
            resultat.num[i] = (uint32_t)(diff + (1LL << 32));
            borrow = 1;
        } else {
            resultat.num[i] = (uint32_t)diff;
            borrow = 0;
        }
    }

    // Pour determiner le taille de resultat 
    while (resultat.Length > 1 && resultat.num[resultat.Length - 1] == 0) {
        resultat.Length--;
    }

    return resultat ;

}

// Afficher 
void Aff_big(bignmb a) {
    printf("0x");
    // On commence par la fin 
    if (a.Length == 0) { printf("0"); return; }
    
    printf("%X", a.num[a.Length - 1]); // Le premier sans les zéros devant
    
    for (int i = a.Length - 2; i >= 0; i--) {
        printf("%08X", a.num[i]); // Les suivants avec padding
    }
    printf("\n");
}

// Fonction auxiliere pour multiplier entre 32 bits nombre te big nombre 
bignmb mult(int32_t n, bignmb a){
    bignmb resultat ;
    for(int i = 0; i < T_MAX; i++) resultat.num[i] = 0 ;
    uint32_t reste = 0;

    // Si on multiplie par 0, on retourne 0 tout de suite
    if (n == 0) {
        resultat.Length = 1;
        resultat.num[0] = 0;
        return resultat;
    }

    int i;
    for (i = 0; i< a.Length ; i ++){
        uint64_t res = (uint64_t)n * a.num[i] + reste;

        resultat.num[i] = (uint32_t)res ;
        reste = (uint32_t)(res >> 32);
    }

    resultat.Length = a.Length;
    if ( reste == 0){
        return resultat;
    } else {
        resultat.num[resultat.Length] = reste ;
        resultat.Length ++ ;
        return resultat ;
    }
    
}

// Fonction auxauliere pour dcalage :
void deccaler (bignmb* a, int pad){
    if (pad + a->Length >= T_MAX){
        printf("Decallage hors zone");
        exit(1);
    }

    for (int i = a->Length -1; i >= 0; i--){
        a->num[i + pad ] = a->num[i];
    }
    a->Length += pad ;
    for (int i = 0; i < pad ; i++){
        a->num[i] = 0 ;
    }
}

// Multilication 
bignmb Mult_big(bignmb a , bignmb b){
    bignmb *grand = (a.Length >= b.Length) ? &a : &b;
    bignmb *petit = (a.Length >= b.Length) ? &b : &a;

    bignmb resultat ;
    for(int i = 0; i < T_MAX; i++) resultat.num[i] = 0 ;
    resultat.Length = 0 ;

    for(int i = 0; i < petit->Length ; i++){
        bignmb res = mult(petit->num[i], (*grand));
        deccaler(&res, i);
        resultat = Add_big(resultat,res);
    }
    return resultat ;

}

// Pour verifier si un big nombre n est nulle ou
int nulle(bignmb a){
    if (a.Length == 0) return 1;
    if (a.Length == 1 && a.num[0] == 0) return 1;
    return 0;
}


// Fonction auxiliere poru calculer le bit efficace 
int bit_mvc(bignmb a){
    if (nulle(a)) return 0;
    
    // On regarde le dernier bloc
    uint32_t top = a.num[a.Length - 1];
    int bits_in_top = 0;
    
    // On cherche la position du bit le plus fort dans ce bloc
    while (top > 0) {
        top >>= 1; // Décalage à droite
        bits_in_top++;
    }
    
    // Total = (Blocs complets * 32) + bits du dernier bloc
    return (a.Length - 1) * 32 + bits_in_top;
}

// Décale 'a' vers la gauche de 'shift' bits (équivalent à a * 2^shift)
bignmb Decalage_Gauche(bignmb a, int shift) {

    if (shift == 0) return a;
    bignmb res = Int_big(0);
    
    int dec_blocs = shift / 32;      // Nombre de blocs complets de décalage
    int dec_bits = shift % 32;       // Nombre de bits restants
    
    res.Length = a.Length + dec_blocs;
    if (dec_bits > 0) res.Length++; // Une retenue possible peut augmenter la taille
    
    // Attention aux dépassements de tableau
    if (res.Length > T_MAX) {
        printf("Erreur: Depassement capacite max\n");
        exit(1);
    }

    uint32_t carry = 0;
    
    // On remplit les blocs de décalage avec des zéros
    for (int i = 0; i < dec_blocs; i++) res.num[i] = 0;
    
    // On traite les anciens blocs
    for (int i = 0; i < a.Length; i++) {
        uint64_t val = (uint64_t)a.num[i] << dec_bits; // On décale les bits
        res.num[i + dec_blocs] = (uint32_t)val | carry; // On met la partie basse + retenue précédente
        carry = (uint32_t)(val >> 32); // On garde la partie haute comme retenue
    }
    
    // S'il reste une retenue à la fin
    if (carry > 0) {
        res.num[a.Length + dec_blocs] = carry;
    } else {
        // Ajustement si la longueur a été surestimée
        if (res.num[res.Length - 1] == 0 && res.Length > 1) res.Length--;
    }
    
    return res;
}

  
bignmb Modul_big(bignmb a, bignmb n) {

    if (nulle(n)) {
        printf("Erreur: Division par 0\n");
        exit(1);
    }
    
    // Si a < n, le reste est a
    if (Comp_big(a, n) == 0) return a;
    if (Comp_big(a, n) == 2) return Int_big(0);

    bignmb reste = a;
    
    while (Comp_big(reste, n) >= 1) {
        
        //Calculer la différence de taille en bits
        int bits_r = bit_mvc(reste);
        int bits_n = bit_mvc(n);
        int diff = bits_r - bits_n;
        
        // Sécurité
        if (diff < 0) diff = 0;

        // Créer une version de n décalée 
        bignmb temp = Decalage_Gauche(n, diff);
        
        // Si temp est trop grand  on recule d'un cran
        if (Comp_big(temp, reste) == 1) {
            temp = Decalage_Gauche(n, diff - 1);
        }

        reste = Sous_big(reste, temp);
    }

    return reste;
} 


// Foction auxiliere :
int Est_pair(bignmb n){
    if (nulle(n)) return 1; 
    return (n.num[0] % 2) == 0; // Retourne 1 si le reste est 0
}

// Diviser par deux 
bignmb Div_2(bignmb a) {
    bignmb res;
    for(int i=0; i<T_MAX; i++) res.num[i]=0;
    if (nulle(a)) { res.Length=1; return res; }
    
    res.Length = a.Length;
    uint32_t carry = 0;
    for (int i = a.Length - 1; i >= 0; i--) {
        uint32_t val = a.num[i];
        res.num[i] = (val >> 1) | (carry << 31);
        carry = val & 1;
    }
    while (res.Length > 1 && res.num[res.Length - 1] == 0) res.Length--;
    return res;
}

// Puissance modulaire  a^(n) mod m
bignmb Puis_big(bignmb a, bignmb n, bignmb m){
    bignmb res = Int_big(1);
    
    // Sécurité : on commence par réduire la base
    bignmb base = Modul_big(a, m); 
    bignmb exp = n;

    while (!nulle(exp)){
        // Si exposant impair
        if (!Est_pair(exp)) {
            res = Mult_big(res, base);
            res = Modul_big(res, m); 
        }
        
        base = Mult_big(base, base);
        base = Modul_big(base, m);
        
        exp = Div_2(exp);
    }
    return res;
}

// Generer un nombre Aleatoire en C 
bignmb Ale_big(){
    FILE *fp;

    bignmb resultat ; 
    for (int k = 0; k < T_MAX; k++) resultat.num[k] = 0;

    // Ouverture en lecture binaire 
    fp = fopen("/dev/urandom", "rb");
    if (fp == NULL) {
        printf("Erreur d'ouverture de /dev/urandom\n");
        exit(1);
    }

    // Lire les valeur aleatoire 
    size_t lus = fread(resultat.num, sizeof(uint32_t), T_MAX, fp);

    if (lus != T_MAX) {
        fprintf(stderr, "Erreur : lecture partielle depuis /dev/urandom\n");
        exit(1);
    }

    fclose(fp);

    int i = 32 ;
    while ((resultat.num[i] == 0)&&(i != 0)){
        i -- ;
    }
    resultat.Length = i + 1 ;
    return resultat ;
}