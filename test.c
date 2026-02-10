#include <stdio.h>
#include <assert.h>
#include "bigNmb.h"

void test_addition() {
    bignmb a = Int_big(100);
    bignmb b = Int_big(50);
    bignmb res = Add_big(a, b);
    
    // 100 + 50 = 150
    // On vérifie le premier bloc
    assert(res.num[0] == 150);
    printf("[OK] Addition\n");
}

void test_multiplication() {
    bignmb a = Int_big(10);
    bignmb b = Int_big(10);
    bignmb res = Mult_big(a, b);
    
    assert(res.num[0] == 100);
    printf("[OK] Multiplication\n");
}

void test_modulo() {
    bignmb a = Int_big(100);
    bignmb n = Int_big(30);
    // 100 mod 30 = 10
    bignmb res = Modul_big(a, n);
    
    assert(res.num[0] == 10);
    printf("[OK] Modulo\n");
}

void test_puissance() {
    // 2^10 mod 1000 = 1024 mod 1000 = 24
    bignmb mod = Ale_big();
    bignmb exp = Int_big(10);
    bignmb base = Add_big(mod, Int_big(9));
    
    bignmb res = Puis_big(base, exp, mod);
    assert(res.num[0] == 3486784401);
    printf("[OK] Puissance\n");
}

void test_aleatoire() {
    bignmb rnd = Ale_big();
    // Juste vérifier qu'il ne crash pas et Length > 0
    assert(rnd.Length > 0);
    printf("[OK] Aleatoire\n");
}

int main() {
    printf("--- Lancement des tests ---\n");
    test_addition();
    test_multiplication();
    test_modulo();
    test_puissance();
    test_aleatoire();
    printf("--- Tous les tests sont valides ---\n");
    return 0; // Retourner 0 signifie SUCCÈS pour CI/CD
}