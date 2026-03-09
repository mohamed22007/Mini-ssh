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

bignmb chiffrer_bloc (bignmb M, bignmb e, bignmb n, bignmb R2){
    return Puiss_big(M,e,n,R2);
}

bignmb dechiffrement_bloc (bignmb C, rsa_keys key){
    bignmb R2p = Calculer_R2_Mod_N(key.p);
    bignmb R2q = Calculer_R2_Mod_N(key.q);
    bignmb m1 = Puiss_big(C,key.dp, key.p,R2p);
    bignmb m2 = Puiss_big(C,key.dq,key.q,R2q);
    bignmb h;
    if (Comp_big(m1,m2) > 1){
        h = Montgomery_Mult(key.qinv,Sub_big((Add_big(m1,key.p)),m2),key.p);
    } else {
        h = Montgomery_Mult(key.qinv,Sub_big(m1,m2),key.p);
    }
    bignmb res = Add_big(m2,Mult_big(h,key.q));
    free_big(m1);
    free_big(m2);
    free_big(h);
    return res ;
}

// Convertit un tableau d'octets (Big-Endian) en bignmb
bignmb bytes_to_big(const uint8_t* data, size_t len) {
    size_t num_words = (len + 7) / 8;
    bignmb res = new_big(0);
    res[0] = num_words;
    
    int byte_idx = len - 1;
    for (size_t i = 1; i <= num_words; i++) {
        uint64_t word = 0;
        for (int j = 0; j < 8; j++) {
            if (byte_idx >= 0) {
                word |= ((uint64_t)data[byte_idx]) << (j * 8);
                byte_idx--;
            }
        }
        res[i] = word;
    }
    
    // Nettoyer les zéros de tête
    while(res[0] > 1 && res[res[0]] == 0) res[0]--;
    return res;
}

// Applique le Padding PKCS#1 v1.5
void appliquer_padding_pkcs1(uint8_t *dest, const uint8_t *data_brute, size_t taille_lue) {
    memset(dest, 0, RSA_KEY_SIZE);
    dest[0] = 0x00;
    dest[1] = 0x02;
    
    int fin_ps = RSA_KEY_SIZE - taille_lue - 1;
    for (int i = 2; i < fin_ps; i++) {
        uint8_t r;
        do {
            r = rand() % 256; // Idéalement, utilisez /dev/urandom ici pour plus de sécurité
        } while (r == 0); // Pas de 0 dans le padding !
        dest[i] = r;
    }
    
    dest[fin_ps] = 0x00; // Séparateur
    memcpy(dest + fin_ps + 1, data_brute, taille_lue);
}

// Retire le Padding PKCS#1 v1.5 et récupère les données
int retirer_padding_pkcs1(const uint8_t *padded_data, size_t padded_len, uint8_t *dest_brute) {
    int idx = 0;
    
    // big_to_bytes enlève souvent le 0x00 de tête, on s'adapte :
    if (padded_data[idx] == 0x00) idx++;
    
    if (padded_data[idx] != 0x02) return -1; // Erreur de format
    idx++;
    
    // Avancer jusqu'au séparateur 0x00
    while (idx < (int)padded_len && padded_data[idx] != 0x00) {
        idx++;
    }
    
    if (idx == (int)padded_len) return -1; // Séparateur non trouvé
    idx++; // Sauter le 0x00
    
    int data_len = padded_len - idx;
    memcpy(dest_brute, padded_data + idx, data_len);
    
    return data_len;
}


// Chiffrement d'un fichier complet par blocs
int chiffrer_fichier(const char* fichier_in, const char* fichier_out, bignmb e, bignmb n) {
    FILE* f_in = fopen(fichier_in, "rb");
    FILE* f_out = fopen(fichier_out, "wb");
    
    if (!f_in || !f_out) {
        perror("Erreur d'ouverture des fichiers");
        if (f_in) fclose(f_in);
        if (f_out) fclose(f_out);
        return -1;
    }

    bignmb R2 = Calculer_R2_Mod_N(n); // Précalcul une seule fois pour tout le fichier !
    
    uint8_t buffer_lu[BLOC_SIZE];
    uint8_t bloc_padde[RSA_KEY_SIZE];
    size_t octets_lus;


    // Streaming : Boucle de lecture du fichier
    while ((octets_lus = fread(buffer_lu, 1, BLOC_SIZE, f_in)) > 0) {
        
        // Ajouter le padding
        appliquer_padding_pkcs1(bloc_padde, buffer_lu, octets_lus);
        
        // Convertir en BigInt
        bignmb M = bytes_to_big(bloc_padde, RSA_KEY_SIZE);
        
        // Chiffrement mathématique
        bignmb C = chiffrer_bloc(M, e, n, R2);
        
        // Convertir le chiffré en octets
        uint8_t* c_bytes = NULL;
        size_t c_len = big_to_bytes(C, &c_bytes);
        
        // Ajuster la taille pour toujours écrire RSA_KEY_SIZE octets
        uint8_t bloc_final[RSA_KEY_SIZE];
        memset(bloc_final, 0, RSA_KEY_SIZE);
        if (c_len <= RSA_KEY_SIZE) {
            memcpy(bloc_final + (RSA_KEY_SIZE - c_len), c_bytes, c_len);
        }

        // Écriture dans le fichier
        fwrite(bloc_final, 1, RSA_KEY_SIZE, f_out);
        
        // Nettoyage mémoire du bloc
        free_big(M);
        free_big(C);
        free(c_bytes);
    }

    free_big(R2);
    fclose(f_in);
    fclose(f_out);
    return 0;
}


// Déchiffrement d'un fichier complet par blocs
int dechiffrer_fichier(const char* fichier_in, const char* fichier_out, rsa_keys cles_privees) {
    FILE* f_in = fopen(fichier_in, "rb");
    FILE* f_out = fopen(fichier_out, "wb");
    
    if (!f_in || !f_out) {
        perror("Erreur d'ouverture des fichiers");
        if (f_in) fclose(f_in);
        if (f_out) fclose(f_out);
        return -1;
    }

    // Le fichier chiffré contient des blocs de taille stricte RSA_KEY_SIZE
    uint8_t buffer_lu[RSA_KEY_SIZE];
    uint8_t buffer_clair[RSA_KEY_SIZE];
    size_t octets_lus;


    // Lecture par blocs de la taille de la clé (256 octets)
    while ((octets_lus = fread(buffer_lu, 1, RSA_KEY_SIZE, f_in)) > 0) {
        
        if (octets_lus != RSA_KEY_SIZE) {
            fprintf(stderr, "Erreur : Fichier chiffré corrompu (bloc incomplet).\n");
            break;
        }

        // onvertir en BigInt
        bignmb C = bytes_to_big(buffer_lu, RSA_KEY_SIZE);
        
        // Déchiffrement mathématique (Optimisé CRT)
        bignmb M = dechiffrement_bloc(C, cles_privees);
        
        // Convertir le clair en octets
        uint8_t* m_bytes = NULL;
        size_t m_len = big_to_bytes(M, &m_bytes);
        
        // Retirer le padding
        int donnees_reelles_len = retirer_padding_pkcs1(m_bytes, m_len, buffer_clair);
        
        if (donnees_reelles_len < 0) {
            fprintf(stderr, "Erreur de déchiffrement : Padding invalide.\n");
            free_big(C); free_big(M); free(m_bytes);
            break; // Clé incorrecte ou fichier corrompu
        }

        // Écrire le texte clair dans le fichier de sortie
        fwrite(buffer_clair, 1, donnees_reelles_len, f_out);
        
        // Nettoyage mémoire du bloc
        free_big(C);
        free_big(M);
        free(m_bytes);
    }

    fclose(f_in);
    fclose(f_out);
    return 0;
}


// Décodeur Base64 (Ignore les retours à la ligne du PEM)
static int b64_index(char c) {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (c >= '0' && c <= '9') return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
}

size_t base64_decode(const char* input, uint8_t* output) {
    size_t in_len = strlen(input);
    size_t out_len = 0;
    uint32_t triple = 0;
    int count = 0;

    for (size_t i = 0; i < in_len; i++) {
        if (input[i] == '=' || input[i] == '\n' || input[i] == '\r') continue;
        
        int val = b64_index(input[i]);
        if (val == -1) continue;

        triple = (triple << 6) | val;
        count++;

        if (count == 4) {
            output[out_len++] = (triple >> 16) & 0xFF;
            output[out_len++] = (triple >> 8) & 0xFF;
            output[out_len++] = triple & 0xFF;
            triple = 0;
            count = 0;
        }
    }

    if (count == 3) {
        triple <<= 6;
        output[out_len++] = (triple >> 16) & 0xFF;
        output[out_len++] = (triple >> 8) & 0xFF;
    } else if (count == 2) {
        triple <<= 12;
        output[out_len++] = (triple >> 16) & 0xFF;
    }
    return out_len;
}


// Lit un uint32_t en Big-Endian (pour format SSH)
static uint32_t read_u32_be(const uint8_t* buf) {
    return (buf[0] << 24) | (buf[1] << 16) | (buf[2] << 8) | buf[3];
}

// Lit la longueur d'un champ ASN.1 DER (pour format PEM)
static size_t read_der_len(const uint8_t* buf, size_t* offset) {
    uint8_t first = buf[(*offset)++];
    if (first < 128) return first;
    
    int num_bytes = first & 0x7F;
    size_t len = 0;
    for (int i = 0; i < num_bytes; i++) {
        len = (len << 8) | buf[(*offset)++];
    }
    return len;
}


// Importe la clé publique depuis un fichier OpenSSH
int import_public_key_ssh(const char* filename, bignmb* e_out, bignmb* n_out) {
    FILE* f = fopen(filename, "r");
    if (!f) return -1;

    char b64_buffer[8192] = {0};
    
    // On ignore le premier mot ("ssh-rsa") et on lit le Base64
    if (fscanf(f, "%*s %8191s", b64_buffer) != 1) {
        fclose(f);
        return -1;
    }
    fclose(f);

    uint8_t blob[8192];
    size_t blob_len = base64_decode(b64_buffer, blob);
    size_t offset = 0;

    // Lire et ignorer le type "ssh-rsa"
    uint32_t type_len = read_u32_be(blob + offset); offset += 4 + type_len;
    
    // Lire l'exposant 'e'
    uint32_t e_len = read_u32_be(blob + offset); offset += 4;
    *e_out = bytes_to_big(blob + offset, e_len); offset += e_len;

    // Lire le modulo 'n'
    uint32_t n_len = read_u32_be(blob + offset); offset += 4;
    *n_out = bytes_to_big(blob + offset, n_len);

    return 0;
}

// Fonction locale pour lire un INTEGER dans l'ASN.1
    bignmb read_next_integer(const uint8_t* buffer, size_t* off) {
        if (buffer[(*off)++] != 0x02) return NULL; // Vérifier le tag INTEGER
        size_t len = read_der_len(buffer, off);
        bignmb res = bytes_to_big(buffer + *off, len);
        *off += len;
        return res;
    }
    


// Importe la clé privée depuis un fichier PEM
rsa_keys import_private_key_pem(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) exit(0);

    char b64_buffer[16384] = {0};
    char line[256];
    
    // Lecture du fichier en ignorant les entêtes -----BEGIN... et -----END...
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '-') {
            strncat(b64_buffer, line, sizeof(b64_buffer) - strlen(b64_buffer) - 1);
        }
    }
    fclose(f);

    uint8_t der[16384];
    size_t der_len = base64_decode(b64_buffer, der);
    size_t offset = 0;

    // Vérifier le tag SEQUENCE (0x30)
    if (der[offset++] != 0x30) exit(0);
    read_der_len(der, &offset); // On passe la longueur globale


    // Extraction dans l'ordre strict de la norme PKCS#1
    bignmb version = read_next_integer(der, &offset); 
    free_big(version);
    rsa_keys keys ;
    
    keys.n    = read_next_integer(der, &offset);
    keys.e    = read_next_integer(der, &offset);
    keys.d    = read_next_integer(der, &offset);
    keys.p    = read_next_integer(der, &offset);
    keys.q    = read_next_integer(der, &offset);
    keys.dp   = read_next_integer(der, &offset);
    keys.dq   = read_next_integer(der, &offset);
    keys.qinv = read_next_integer(der, &offset);

    return keys;
}