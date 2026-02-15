#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "bigNmb.h"

double get_time_sec() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

void test_add(bignmb a, bignmb b) {
    int iter = 1000000;
    double start = get_time_sec();
    for(int i = 0; i < iter; i++) {
        bignmb r = Add_big(a, b);
        free_big(r);
    }
    double end = get_time_sec();
    printf("ADD  : %.4f sec | %.2f ns/op\n",
           end - start,
           ((end - start)/iter)*1e9);
}

void test_sub(bignmb a, bignmb b) {
    if (Comp_big(a,b) == -1) {
        bignmb tmp = a; a = b; b = tmp;
    }

    int iter = 1000000;
    double start = get_time_sec();
    for(int i = 0; i < iter; i++) {
        bignmb r = Sub_big(a, b);
        free_big(r);
    }
    double end = get_time_sec();
    printf("SUB  : %.4f sec | %.2f ns/op\n",
           end - start,
           ((end - start)/iter)*1e9);
}

void test_mult(bignmb a, bignmb b) {
    int iter = 100000;
    double start = get_time_sec();
    for(int i = 0; i < iter; i++) {
        bignmb r = Mult_big(a, b);
        free_big(r);
    }
    double end = get_time_sec();
    printf("MULT : %.4f sec | %.2f us/op\n",
           end - start,
           ((end - start)/iter)*1e6);
}

void test_pow(bignmb a, bignmb b) {
    bignmb n = gen_aleatoire();
    n[1] |= 1ULL;

    double start = get_time_sec();
    bignmb R2 = Calculer_R2_Mod_N(n);
    double end = get_time_sec();

    printf("R2   : %.4f sec (precalc)\n", end - start);

    int iter = 100;
    start = get_time_sec();
    for(int i = 0; i < iter; i++) {
        bignmb r = Puiss_big(a, b, n, R2);
        free_big(r);
    }
    end = get_time_sec();

    printf("POW  : %.4f sec | %.2f ms/op\n",
           end - start,
           ((end - start)/iter)*1000);

    free_big(n);
    free_big(R2);
}

void test_miller() {
    bignmb p = new_big(0);
    p[0] = 1;
    p[1] = 0x1FFFFFFFFFFFFFFF;

    bignmb c = new_big(0);
    c[0] = 1;
    c[1] = 0x1FFFFFFFFFFFFFFE;

    double start = get_time_sec();
    Miller_Rabin(p,10);
    double end = get_time_sec();
    printf("MR(P): %.6f sec\n", end - start);

    start = get_time_sec();
    Miller_Rabin(c,10);
    end = get_time_sec();
    printf("MR(C): %.6f sec\n", end - start);

    free_big(p);
    free_big(c);
}

void test_prime_gen() {
    double start, end;

    start = get_time_sec();
    bignmb p = Gen_premier();
    end = get_time_sec();
    printf("P = "); print_big(p);
    printf("Time P : %.4f sec\n", end - start);

    start = get_time_sec();
    bignmb q = Gen_premier();
    end = get_time_sec();
    printf("Q = "); print_big(q);
    printf("Time Q : %.4f sec\n", end - start);

    free_big(p);
    free_big(q);
}

int main() {
    srand(time(NULL));

    bignmb a = gen_aleatoire();
    bignmb b = gen_aleatoire();

    test_add(a,b);
    test_sub(a,b);
    test_mult(a,b);
    test_pow(a,b);
    test_miller();
    test_prime_gen();

    free_big(a);
    free_big(b);

    return 0;
}
