#include "rsa.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// ============================================================================
// PARTIE 1 : GENERATION DES CLES ET PARAMETRES CRT
// ============================================================================

rsa_keys generate_rsa_keys() {
    rsa_keys keys;
    bignmb un = new_big(1);
    bignmb p_minus_1, q_minus_1, phi_n;

    // Regénérer p et q jusqu'à ce que gcd(e=65537, phi_n) == 1
    // 65537 est premier donc gcd(e, phi_n)=1  <=>  phi_n % 65537 != 0
    while (1) {
        keys.p = Gen_premier();
        keys.q = Gen_premier();

        // S'assurer que p > q pour le calcul de qinv
        if (Comp_big(keys.p, keys.q) < 0) {
            bignmb tmp = keys.p; keys.p = keys.q; keys.q = tmp;
        }

        p_minus_1 = Sub_big(keys.p, un);
        q_minus_1 = Sub_big(keys.q, un);
        phi_n = Mult_big(p_minus_1, q_minus_1);

        if (big_mod_small(phi_n, 65537) != 0) break;

        free_big(keys.p); free_big(keys.q);
        free_big(p_minus_1); free_big(q_minus_1); free_big(phi_n);
    }

    // n = p * q
    keys.n = Mult_big(keys.p, keys.q);

    // e = 65537
    keys.e = new_big(65537);

    // d = inverse(e, phi(n))
    keys.d = inverse(keys.e, phi_n);

    // Paramètres CRT
    keys.dp   = Mod_big(keys.d, p_minus_1);
    keys.dq   = Mod_big(keys.d, q_minus_1);
    keys.qinv = inverse(keys.q, keys.p);

    free_big(un); free_big(p_minus_1); free_big(q_minus_1); free_big(phi_n);
    return keys;
}

void free_rsa_keys(rsa_keys* keys) {
    if (keys->n) free_big(keys->n);
    if (keys->e) free_big(keys->e);
    if (keys->d) free_big(keys->d);
    if (keys->p) free_big(keys->p);
    if (keys->q) free_big(keys->q);
    if (keys->dp) free_big(keys->dp);
    if (keys->dq) free_big(keys->dq);
    if (keys->qinv) free_big(keys->qinv);
}

// ============================================================================
// PARTIE 3 : UTILITAIRES DE FORMATAGE BINAIRE ET BASE64
// ============================================================================

// Encodeur Base64 complet et standardisé
static const char base64_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

void base64_encode(const uint8_t* data, size_t input_length, char* encoded_data) {
    int i = 0, j = 0;
    uint32_t octet_a, octet_b, octet_c, triple;

    for (i = 0; i < (int)input_length; i += 3) {
        octet_a = i < input_length ? data[i] : 0;
        octet_b = i + 1 < input_length ? data[i + 1] : 0;
        octet_c = i + 2 < input_length ? data[i + 2] : 0;

        triple = (octet_a << 16) + (octet_b << 8) + octet_c;

        encoded_data[j++] = base64_chars[(triple >> 3 * 6) & 0x3F];
        encoded_data[j++] = base64_chars[(triple >> 2 * 6) & 0x3F];
        encoded_data[j++] = (i + 1 < input_length) ? base64_chars[(triple >> 1 * 6) & 0x3F] : '=';
        encoded_data[j++] = (i + 2 < input_length) ? base64_chars[(triple >> 0 * 6) & 0x3F] : '=';
    }
    encoded_data[j] = '\0';
}
// Convertit un bignmb en tableau d'octets (Big-Endian)
size_t big_to_bytes(bignmb a, uint8_t** buffer) {
    size_t size = a[0] * 8; 
    
    // On alloue avec +2 pour garantir de la place pour le 0x00 de signe
    *buffer = (uint8_t*)calloc(size + 2, sizeof(uint8_t));
    
    int index = 1; // ON COMMENCE À 1 (La sécurité qui évite le Segfault !)
    for (int i = a[0]; i >= 1; i--) {
        for (int j = 7; j >= 0; j--) {
            (*buffer)[index++] = (a[i] >> (j * 8)) & 0xFF;
        }
    }

    // Supprimer les zéros de tête
    int start = 1;
    while (start < index - 1 && (*buffer)[start] == 0) start++;

    // Si le bit de poids fort est 1, on ajoute un 0x00 devant en reculant de 1 case
    if ((*buffer)[start] & 0x80) {
        start--; // Grâce à notre marge, start devient 0 (et pas -1 !)
        (*buffer)[start] = 0x00;
    }

    size_t final_len = index - start;
    memmove(*buffer, *buffer + start, final_len);
    return final_len;
}

// ============================================================================
// PARTIE 4 : EXPORT DES CLES (SSH et PEM)
// ============================================================================

// Écrit un uint32 en Big-Endian dans un buffer
void write_u32_be(uint8_t* buf, uint32_t val) {
    buf[0] = (val >> 24) & 0xFF; buf[1] = (val >> 16) & 0xFF;
    buf[2] = (val >> 8) & 0xFF;  buf[3] = val & 0xFF;
}

// Exporte la clé publique au format id_rsa.pub (OpenSSH)
int export_public_key_ssh(rsa_keys* keys, const char* filename, const char* comment) {
    FILE* f = fopen(filename, "w");
    if (!f) return -1;

    uint8_t* e_bytes; size_t e_len = big_to_bytes(keys->e, &e_bytes);
    uint8_t* n_bytes; size_t n_len = big_to_bytes(keys->n, &n_bytes);

    // Structure SSH : [Len(7) | "ssh-rsa" | Len(e) | e | Len(n) | n]
    size_t total_len = 4 + 7 + 4 + e_len + 4 + n_len;
    uint8_t* blob = (uint8_t*)malloc(total_len);
    size_t offset = 0;

    write_u32_be(blob + offset, 7); offset += 4;
    memcpy(blob + offset, "ssh-rsa", 7); offset += 7;
    
    write_u32_be(blob + offset, e_len); offset += 4;
    memcpy(blob + offset, e_bytes, e_len); offset += e_len;

    write_u32_be(blob + offset, n_len); offset += 4;
    memcpy(blob + offset, n_bytes, n_len); offset += n_len;

    // Encodage Base64
    char* b64 = (char*)malloc(total_len * 2);
    base64_encode(blob, total_len, b64);

    fprintf(f, "ssh-rsa %s %s\n", b64, comment);

    free(e_bytes); free(n_bytes); free(blob); free(b64);
    fclose(f);
    return 0;
}

// Exporte la clé privée au format id_rsa (PEM / ASN.1 DER)
int export_private_key_pem(rsa_keys* keys, const char* filename) {
    // La construction d'une structure ASN.1 DER est complexe. 
    // Il faut concaténer [0x02, Len, Value] pour Version(0), n, e, d, p, q, dp, dq, qinv
    // Puis englober le tout dans une Sequence [0x30, TotalLen, Datas...]
    
    // À ce stade de ton projet, tu devras construire un buffer et l'encoder en Base64.
    // L'architecture de la fonction ressemblera à l'export SSH, mais avec les tags ASN.1.
    
    FILE* f = fopen(filename, "w");
    if (!f) return -1;
    fprintf(f, "-----BEGIN RSA PRIVATE KEY-----\n");
    // fprintf(f, "%s\n", base64_der_blob);
    fprintf(f, "-----END RSA PRIVATE KEY-----\n");
    fclose(f);
    return 0;
}