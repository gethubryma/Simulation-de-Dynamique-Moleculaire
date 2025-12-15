
# Simulation de Dynamique Moléculaire Lennard-Jones (A/B)

## 1. Présentation générale

Ce projet implémente un **moteur de dynamique moléculaire (Molecular Dynamics, MD)** en trois dimensions destiné à l’étude d’un système de particules interagissant par un **potentiel de Lennard-Jones**.  
Le système simulé est **binaire** (particules de type A et une particule de type B), confiné dans une **boîte cubique périodique**, avec contrôle thermique assuré par un **thermostat de Berendsen**.

Le code est écrit en **C++**, exploite le **parallélisme OpenMP** pour le calcul des forces, et génère des fichiers **PDB** permettant la visualisation directe des trajectoires atomiques.

---

## 2. Objectifs scientifiques et numériques

- Implémenter une dynamique moléculaire classique conforme aux standards.
- Modéliser des interactions **A–A** et **B–A** par un potentiel de Lennard-Jones.
- Étudier l’évolution temporelle des **énergies** et de la **température**.
- Vérifier la **conservation du momentum** du centre de masse.
- Exploiter le calcul parallèle pour améliorer les performances.
- Produire des trajectoires exploitables pour l’analyse structurale.

---

## 3. Modèle physique

### 3.1 Système simulé

- Nombre total de particules : $( N )$
- Espace : $( \mathbb{R}^3 )$
- Boîte cubique de taille :

$$
L = 42.0
$$

- Conditions aux limites : **périodiques**
- Degrés de liberté :

$$
N_{dof} = 3N - 3
$$

(les trois degrés correspondant au mouvement global du centre de masse sont retirés)

---

### 3.2 Potentiel de Lennard-Jones

Les interactions entre particules sont décrites par le potentiel :

$$
U(r) = 4\varepsilon \left[\left(\frac{r^\*}{r}\right)^{12} - 2\left(\frac{r^\*}{r}\right)^6 \right]
$$

où :
- $( r )$ est la distance inter-particulaire,
- $( \varepsilon )$ contrôle la profondeur du puits de potentiel,
- $( r^\* )$ est la distance caractéristique.

#### Paramètres selon le type d’interaction

- **A–A** :

$$
\varepsilon = \varepsilon^\*, \quad r^\* = r^\*
$$

- **B–A** :

$$
\varepsilon = 1.0, \quad r^\* = 3.5
$$

Un **rayon de coupure** est appliqué :

$$
r \leq R_{cut} = 10.0
$$

---

### 3.3 Forces inter-particulaires

La force exercée sur une particule est donnée par le gradient du potentiel :

$$
\vec{F}_{ij} = -\nabla U(r_{ij})
$$

Dans le cas du potentiel de Lennard-Jones :

$$
\vec{F}_{ij} = 24\varepsilon \left[ 2\left(\frac{r^\*}{r}\right)^{12} - \left(\frac{r^\*}{r}\right)^6 \right] \frac{\vec{r}_{ij}}{r^2}
$$

---

## 4. Conditions aux limites périodiques

Les positions sont repliées dans la boîte selon :

$$
x \leftarrow x - \lfloor x / L \rfloor \cdot L
$$

Le calcul des forces prend en compte les **images périodiques** à l’aide d’une maille :

$$
3 \times 3 \times 3
$$

afin de respecter la symétrie du système.

---

## 5. Dynamique et intégration temporelle

### 5.1 Équations du mouvement

Les équations de la mécanique newtonienne sont utilisées :

$$
m \frac{d^2 \vec{r}_i}{dt^2} = \vec{F}_i
$$

---

### 5.2 Schéma Velocity-Verlet

L’intégration temporelle est assurée par le schéma **Velocity-Verlet** :

1. Mise à jour intermédiaire des vitesses :

$$
\vec{v}_i\left(t + \frac{\Delta t}{2}\right) =
\vec{v}_i(t) + \frac{\Delta t}{2m}\vec{F}_i(t)
$$

2. Mise à jour des positions :

$$
\vec{r}_i(t + \Delta t) =
\vec{r}_i(t) + \vec{v}_i\left(t + \frac{\Delta t}{2}\right)\Delta t
$$

3. Mise à jour finale des vitesses :

$$
\vec{v}_i(t + \Delta t) =
\vec{v}_i\left(t + \frac{\Delta t}{2}\right) +
\frac{\Delta t}{2m}\vec{F}_i(t + \Delta t)
$$

Le pas de temps utilisé est :

$$
\Delta t = 0.01
$$

---

## 6. Température et thermostat

### 6.1 Énergie cinétique

L’énergie cinétique totale est donnée par :

$$
E_{kin} = \sum_{i=1}^{N} \frac{p_i^2}{2m}
$$

---

### 6.2 Température instantanée

La température est calculée à partir de l’équipartition de l’énergie :

$$
T = \frac{E_{kin}}{N_{dof} R}
$$

où $( R )$ est la constante des gaz parfaits.

---

### 6.3 Thermostat de Berendsen

Le thermostat ajuste les momenta par le facteur :

$$
\lambda = \gamma \left( \frac{T_{cible}}{T_{actuelle}} - 1 \right)
$$

Les momenta sont ensuite modifiés selon :

$$
\vec{p}_i \leftarrow (1 + \lambda)\vec{p}_i
$$

Une phase de **chauffage progressif** est appliquée durant les premiers pas de simulation.

---

## 7. Initialisation du système

- Lecture des positions initiales depuis `particule.xyz`
- Génération aléatoire des momenta initiaux
- Ajustement des momenta pour imposer une température initiale
- Suppression du mouvement global du centre de masse :

$$
\sum_i \vec{p}_i = 0
$$

---

## 8. Calculs énergétiques spécifiques

### 8.1 Énergie moyenne A–A

Pour les particules de type A: 

$$
E_{AA}^{\mathrm{moy}} = \frac{2}{N_A} \sum_{i < j} U_{ij}
$$


---

### 8.2 Énergie B–A

Pour la particule B unique :

$$
E_{BA} = \sum_{i \in A} U_{Bi}
$$

---

## 9. Sorties et visualisation

### 9.1 Affichage console

- Énergie potentielle
- Énergie cinétique
- Énergie totale
- Température
- Momentum du centre de masse

---

### 9.2 Fichier PDB

- Écriture périodique des configurations
- Visualisation directe des trajectoires atomiques
- Analyse structurale et dynamique

---

## 10. Structure du projet

```

.
├── main.cpp
├── Simulation.h
├── Simulation.cpp
├── particule.xyz
└── Fichier_sortie.pdb

````

---

## 11. Compilation

```bash
g++ -O2 -fopenmp main.cpp Simulation.cpp -o simulation
````

---

## 12. Exécution

```bash
./simulation
```

Le programme demande :

* le nombre de particules,
* le nombre de pas de simulation.

---

## 13. Conclusion

Ce projet constitue une implémentation rigoureuse et complète d’une dynamique moléculaire classique.
Il met en œuvre les équations fondamentales de la mécanique statistique, un schéma d’intégration robuste, un contrôle thermique maîtrisé et un calcul parallèle efficace, fournissant ainsi une base solide pour l’étude numérique de systèmes particulaires et le développement de simulations plus avancées.

