#ifndef RSA_H
#define RSA_H
#define RSA_KEY_SIZE 512        // n = 4096 bits = 512 octets (T_LIMIT=32 => p,q de 2048 bits)
#define BLOC_SIZE (RSA_KEY_SIZE - 11)  // Max PKCS#1 v1.5 : RSA_KEY_SIZE - 11 = 501 octets

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

// Fonction pour chiffrer le message par bloque 
bignmb chiffrer_bloc (bignmb M, bignmb e, bignmb n, bignmb R2);

// Fonction Pour dechiffrement 
bignmb dechiffrement_bloc (bignmb C, rsa_keys key);

// Déchiffrement d'un fichier complet par blocs
int dechiffrer_fichier(const char* fichier_in, const char* fichier_out, rsa_keys cles_privees);

// Chiffrement d'un fichier complet par blocs
int chiffrer_fichier(const char* fichier_in, const char* fichier_out, bignmb e, bignmb n);

// Importation de la clé publique depuis un fichier OpenSSH
int import_public_key_ssh(const char* filename, bignmb* e_out, bignmb* n_out);

// Importation de la clé privée depuis un fichier PEM (retourne la struct)
rsa_keys import_private_key_pem(const char* filename);

#endif