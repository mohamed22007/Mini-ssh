#include <stdio.h>
#include <stdlib.h>
#include "rsa.h"

int main(int argc, char *argv[]) {
    // Vérification du nombre d'arguments
    if (argc != 3) {
        fprintf(stderr, "Erreur de syntaxe.\n");
        fprintf(stderr, "Usage: %s <fichier_cle_publique> <fichier_cle_privee>\n", argv[0]);
        fprintf(stderr, "Exemple: %s id_rsa.pub id_rsa\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *pub_file = argv[1];
    const char *priv_file = argv[2];

    printf("Génération des clés RSA en cours (cela peut prendre quelques secondes)...\n");
    
    // Appel à ta fonction principale
    rsa_keys keys = generate_rsa_keys();
    
    printf("Clés générées avec succès !\n");

    // Exportation de la clé publique
    printf("Exportation de la clé publique vers : %s\n", pub_file);
    if (export_public_key_ssh(&keys, pub_file, "mon_rsa_custom") == 0) {
        printf(" -> Clé publique exportée.\n");
    } else {
        fprintf(stderr, "Erreur lors de la création du fichier %s\n", pub_file);
    }

    // Exportation de la clé privée
    printf("Exportation de la clé privée vers : %s\n", priv_file);
    if (export_private_key_pem(&keys, priv_file) == 0) {
        printf(" -> Clé privée exportée.\n");
    } else {
        fprintf(stderr, "Erreur lors de la création du fichier %s\n", priv_file);
    }

    // Nettoyage de la mémoire
    free_rsa_keys(&keys);
    
    printf("Programme terminé avec succès.\n");
    return EXIT_SUCCESS;
}