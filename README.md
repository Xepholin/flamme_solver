# Projet Flamme - Solveur d'équation de diffusion non linéaire

Ce projet vise à résoudre une équation de diffusion non linéaire, spécifiquement modélisant une flamme, en utilisant deux méthodes numériques : un schéma implicite linéarisé et la méthode de Newton. Le projet est implémenté en C et utilise la bibliothèque HYPRE pour la résolution de systèmes linéaires.

## Dépendances

- **HYPRE** : Pour la résolution de systèmes linéaires.
- **MPI** : Nécessaire pour HYPRE.

## Compilation

Pour compiler le projet, assurez-vous d'avoir installé les dépendances nécessaires (MPI et HYPRE). Ensuite, suivez les étapes suivantes :
```bash
mkdir build
cd build
cmake ..
make
```

## Exécution

Après la compilation, un exécutable nommé flamme sera généré dans le répertoire build. Vous pouvez l'exécuter avec la commande suivante :
```bash
mpirun -np <nombre_de_processus> ./flamme [fichier_ils] [fichier_newton]

    <nombre_de_processus> : Nombre de processus MPI à utiliser.

    [fichier_ils] : Chemin du fichier de sortie pour la solution du schéma implicite linéarisé (optionnel, par défaut ../solution_ils.dat).

    [fichier_newton] : Chemin du fichier de sortie pour la solution de la méthode de Newton (optionnel, par défaut ../solution_newton.dat).
```
## Résultats

Les solutions sont sauvegardées dans les fichiers spécifiés (ou par défaut dans solution_ils.dat et solution_newton.dat). Ces fichiers contiennent les valeurs de la solution à chaque point de la grille.

## Structure du projet

- **src/main.c** : Point d'entrée du programme. Initialise les paramètres, résout l'équation de diffusion avec les deux méthodes, et sauvegarde les résultats.
- **src/solver.c** : Contient les fonctions pour initialiser la solution, résoudre l'équation de diffusion avec les deux méthodes, et sauvegarder les résultats.
- **src/matrix.c** : Implémente une fonction pour résoudre un système tridiagonal en utilisant HYPRE.
- **include/const.h** : Définit les constantes du projet (dimension, paramètres physiques, etc.).
- **include/solver.h** : Déclare les fonctions de résolution et d'initialisation.
- **include/matrix.h** : Déclare la fonction de résolution de systèmes tridiagonaux.
- **CMakeLists.txt** : Fichier de configuration CMake pour compiler le projet.
