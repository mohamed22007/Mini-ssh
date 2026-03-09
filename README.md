# Mini-ssh

Mini-ssh est un outil cryptographique écrit en C permettant de simuler une communication chiffrée de bout en bout. Il implémente son propre moteur de gestion de grands nombres (BigInt) et l'algorithme RSA depuis zéro, y compris la génération de nombres premiers robustes.

##  Fonctionnalités

1. **Génération de clés RSA** : Génère des paires de clés (publique/privée) avec export aux formats standards (OpenSSH pour la clé publique, PEM/PKCS#1 pour la clé privée).
2. **Chiffrement de fichiers** : Chiffrement par blocs (PKCS#1 v1.5) d'un fichier en utilisant une clé publique.
3. **Déchiffrement de fichiers** : Déchiffrement d'un fichier chiffré `.enc` en utilisant la clé privée correspondante.
4. **Moteur BigInt** : Bibliothèque mathématique développée sur mesure en C et en Assembleur x86_64 pour optimiser l'arithmétique des grands nombres (Addition, Soustraction, Multiplication de Montgomery).

##  Compilation

Le projet utilise `make`. Pour compiler le programme, exécutez simplement :

```bash
make