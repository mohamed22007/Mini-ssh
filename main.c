#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rsa.h"

void print_usage(const char* prog_name) {
    fprintf(stderr, "Erreur de syntaxe.\n");
    fprintf(stderr, "Utilisation de l'outil %s :\n\n", prog_name);
    fprintf(stderr, "  1. Générer des clés :\n");
    fprintf(stderr, "     %s keygen <cle_publique_out> <cle_privee_out>\n", prog_name);
    fprintf(stderr, "     Exemple : %s keygen id_rsa.pub id_rsa\n\n", prog_name);
    
    fprintf(stderr, "  2. Chiffrer un fichier :\n");
    fprintf(stderr, "     %s encrypt <cle_publique_in> <fichier_source> <fichier_destination>\n", prog_name);
    fprintf(stderr, "     Exemple : %s encrypt id_rsa.pub message.txt message.enc\n\n", prog_name);
    
    fprintf(stderr, "  3. Déchiffrer un fichier :\n");
    fprintf(stderr, "     %s decrypt <cle_privee_in> <fichier_source> <fichier_destination>\n", prog_name);
    fprintf(stderr, "     Exemple : %s decrypt id_rsa message.enc message_decrypte.txt\n\n", prog_name);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    const char *command = argv[1];

    // ==========================================
    // COMMANDE 1 : KEYGEN (Génération de clés)
    // ==========================================
    if (strcmp(command, "keygen") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Erreur : arguments manquants pour 'keygen'.\n");
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        const char *pub_file = argv[2];
        const char *priv_file = argv[3];

        printf("Génération des clés RSA en cours (cela peut prendre quelques secondes)...\n");
        rsa_keys keys = generate_rsa_keys();
        printf("Clés générées avec succès !\n");

        if (export_public_key_ssh(&keys, pub_file, "Mssh_key") == 0) {
            printf(" -> Clé publique exportée vers : %s\n", pub_file);
        } else {
            fprintf(stderr, "Erreur lors de l'exportation de la clé publique.\n");
        }

        if (export_private_key_pem(&keys, priv_file) == 0) {
            printf(" -> Clé privée exportée vers : %s\n", priv_file);
        } else {
            fprintf(stderr, "Erreur lors de l'exportation de la clé privée.\n");
        }

        free_rsa_keys(&keys);
        return EXIT_SUCCESS;
    }

    // ==========================================
    // COMMANDE 2 : ENCRYPT (Chiffrement)
    // ==========================================
    else if (strcmp(command, "encrypt") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Erreur : arguments manquants pour 'encrypt'.\n");
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        const char *pub_key_file = argv[2];
        const char *src_file = argv[3];
        const char *dest_file = argv[4];

        bignmb e, n;
        if (import_public_key_ssh(pub_key_file, &e, &n) != 0) {
            fprintf(stderr, "Erreur : Impossible de lire la clé publique %s\n", pub_key_file);
            return EXIT_FAILURE;
        }

        if (chiffrer_fichier(src_file, dest_file, e, n) != 0) {
            fprintf(stderr, "Échec du chiffrement.\n");
            free_big(e); free_big(n);
            return EXIT_FAILURE;
        }

        free_big(e); free_big(n);
        return EXIT_SUCCESS;
    }

    // ==========================================
    // COMMANDE 3 : DECRYPT (Déchiffrement)
    // ==========================================
    else if (strcmp(command, "decrypt") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Erreur : arguments manquants pour 'decrypt'.\n");
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        const char *priv_key_file = argv[2];
        const char *src_file = argv[3];
        const char *dest_file = argv[4];

        rsa_keys keys;
        // On initialise à NULL pour que free_rsa_keys fonctionne même si l'import échoue partiellement
        memset(&keys, 0, sizeof(rsa_keys)); 

        if (import_private_key_pem(priv_key_file, &keys) != 0) {
            fprintf(stderr, "Erreur : Impossible de lire la clé privée %s\n", priv_key_file);
            return EXIT_FAILURE;
        }

        if (dechiffrer_fichier(src_file, dest_file, keys) != 0) {
            fprintf(stderr, "Échec du déchiffrement.\n");
            free_rsa_keys(&keys);
            return EXIT_FAILURE;
        }

        free_rsa_keys(&keys);
        return EXIT_SUCCESS;
    }

    // ==========================================
    // COMMANDE INCONNUE
    // ==========================================
    else {
        fprintf(stderr, "Erreur : Commande inconnue '%s'.\n", command);
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
}