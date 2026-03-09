#include "rsa.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// GENERATION DES CLES ET PARAMETRES CRT
rsa_keys generate_rsa_keys() {
    rsa_keys keys;
    bignmb un = new_big(1);
    bignmb p_minus_1, q_minus_1, phi_n;

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

// UTILITAIRES DE FORMATAGE BINAIRE ET BASE64

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

// EXPORT DES CLES (SSH et PEM)

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

// ENCODAGE ASN.1 DER (pour clé privée PKCS#1)


static size_t der_len_size(size_t len) {
    if (len < 128)   return 1;
    if (len < 256)   return 2;
    if (len < 65536) return 3;
    return 4;
}


static size_t der_write_len(uint8_t* buf, size_t len) {
    if (len < 128) {
        buf[0] = (uint8_t)len;
        return 1;
    } else if (len < 256) {
        buf[0] = 0x81;
        buf[1] = (uint8_t)len;
        return 2;
    } else if (len < 65536) {
        buf[0] = 0x82;
        buf[1] = (uint8_t)(len >> 8);
        buf[2] = (uint8_t)(len & 0xFF);
        return 3;
    } else {
        buf[0] = 0x83;
        buf[1] = (uint8_t)(len >> 16);
        buf[2] = (uint8_t)(len >> 8);
        buf[3] = (uint8_t)(len & 0xFF);
        return 4;
    }
}


static size_t der_write_integer(uint8_t* buf, const uint8_t* data, size_t data_len) {
    size_t offset = 0;
    buf[offset++] = 0x02;                          
    offset += der_write_len(buf + offset, data_len);
    memcpy(buf + offset, data, data_len);
    offset += data_len;
    return offset;
}

static size_t der_integer_size(size_t data_len) {
    return 1 + der_len_size(data_len) + data_len;
}


int export_private_key_pem(rsa_keys* keys, const char* filename) {

    /*  Convertir chaque composante en octets Big-Endian  */
    uint8_t* ver_bytes   = NULL; size_t ver_len   = 0;
    uint8_t* n_bytes     = NULL; size_t n_len     = 0;
    uint8_t* e_bytes     = NULL; size_t e_len     = 0;
    uint8_t* d_bytes     = NULL; size_t d_len     = 0;
    uint8_t* p_bytes     = NULL; size_t p_len     = 0;
    uint8_t* q_bytes     = NULL; size_t q_len     = 0;
    uint8_t* dp_bytes    = NULL; size_t dp_len    = 0;
    uint8_t* dq_bytes    = NULL; size_t dq_len    = 0;
    uint8_t* qinv_bytes  = NULL; size_t qinv_len  = 0;

    /* version = 0 */
    ver_bytes  = (uint8_t*)malloc(1);
    ver_bytes[0] = 0x00;
    ver_len    = 1;

    n_len    = big_to_bytes(keys->n,    &n_bytes);
    e_len    = big_to_bytes(keys->e,    &e_bytes);
    d_len    = big_to_bytes(keys->d,    &d_bytes);
    p_len    = big_to_bytes(keys->p,    &p_bytes);
    q_len    = big_to_bytes(keys->q,    &q_bytes);
    dp_len   = big_to_bytes(keys->dp,   &dp_bytes);
    dq_len   = big_to_bytes(keys->dq,   &dq_bytes);
    qinv_len = big_to_bytes(keys->qinv, &qinv_bytes);

    /* Calculer la taille du contenu de la SEQUENCE  */
    size_t seq_content_len =
        der_integer_size(ver_len)  +
        der_integer_size(n_len)    +
        der_integer_size(e_len)    +
        der_integer_size(d_len)    +
        der_integer_size(p_len)    +
        der_integer_size(q_len)    +
        der_integer_size(dp_len)   +
        der_integer_size(dq_len)   +
        der_integer_size(qinv_len);

    size_t total_der_len = 1 + der_len_size(seq_content_len) + seq_content_len;

    /* ----- 3. Construire le buffer DER ----- */
    uint8_t* der = (uint8_t*)malloc(total_der_len);
    if (!der) {
        free(ver_bytes); free(n_bytes); free(e_bytes); free(d_bytes);
        free(p_bytes); free(q_bytes); free(dp_bytes); free(dq_bytes);
        free(qinv_bytes);
        return -1;
    }

    size_t offset = 0;

    /* Tag SEQUENCE */
    der[offset++] = 0x30;

    /* Longueur de la séquence */
    offset += der_write_len(der + offset, seq_content_len);

    /* Écriture de chaque INTEGER */
    offset += der_write_integer(der + offset, ver_bytes,  ver_len);
    offset += der_write_integer(der + offset, n_bytes,    n_len);
    offset += der_write_integer(der + offset, e_bytes,    e_len);
    offset += der_write_integer(der + offset, d_bytes,    d_len);
    offset += der_write_integer(der + offset, p_bytes,    p_len);
    offset += der_write_integer(der + offset, q_bytes,    q_len);
    offset += der_write_integer(der + offset, dp_bytes,   dp_len);
    offset += der_write_integer(der + offset, dq_bytes,   dq_len);
    offset += der_write_integer(der + offset, qinv_bytes, qinv_len);

    /*  Encoder en Base64 avec retours à la ligne tous les 64 caractères */
    size_t b64_max = (total_der_len / 3 + 1) * 4 + 4;
    char*  b64_raw = (char*)malloc(b64_max);
    base64_encode(der, total_der_len, b64_raw);

    /* Écrire le fichier PEM  */
    FILE* f = fopen(filename, "w");
    if (!f) {
        free(ver_bytes); free(n_bytes); free(e_bytes); free(d_bytes);
        free(p_bytes); free(q_bytes); free(dp_bytes); free(dq_bytes);
        free(qinv_bytes); free(der); free(b64_raw);
        return -1;
    }

    fprintf(f, "-----BEGIN RSA PRIVATE KEY-----\n");

    /* Écriture du Base64 par blocs de 64 caractères (norme PEM) */
    size_t b64_len = strlen(b64_raw);
    for (size_t i = 0; i < b64_len; i += 64) {
        size_t chunk = (b64_len - i < 64) ? (b64_len - i) : 64;
        fwrite(b64_raw + i, 1, chunk, f);
        fputc('\n', f);
    }

    fprintf(f, "-----END RSA PRIVATE KEY-----\n");

    /* Nettoyage  */
    free(ver_bytes); free(n_bytes); free(e_bytes); free(d_bytes);
    free(p_bytes); free(q_bytes); free(dp_bytes); free(dq_bytes);
    free(qinv_bytes); free(der); free(b64_raw);
    fclose(f);
    return 0;
}