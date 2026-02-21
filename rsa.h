#ifndef RSA_H
#define RSA_H

#include "bigNmb.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// Structure contenant tous les paramètres d'une clé RSA complète (Norme PKCS#1)
typedef struct {
    bignmb n;       // Modulus (p * q)
    bignmb e;    
    bignmb d;      
    bignmb p;    
    bignmb q;    
    bignmb dp;      // d mod (p-1)
    bignmb dq;      // d mod (q-1)
    bignmb qinv;    // inverse de q mod p
} rsa_keys;

// Génération de l'ensemble des clés et paramètres CRT
rsa_keys generate_rsa_keys();

// Libération de la mémoire de la structure de clés
void free_rsa_keys(rsa_keys* keys);

// Exportation de la clé publique au format OpenSSH (id_rsa.pub)
int export_public_key_ssh(rsa_keys* keys, const char* filename, const char* comment);

// Exportation de la clé privée au format PEM / PKCS#1 (id_rsa)
int export_private_key_pem(rsa_keys* keys, const char* filename);

// Convertit un BigInt en tableau d'octets Big-Endian
size_t big_to_bytes(bignmb a, uint8_t** buffer);
// Encode un tableau d'octets en Base64
void base64_encode(const uint8_t* data, size_t input_length, char* encoded_data);

#endif 