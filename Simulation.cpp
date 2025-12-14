#include "Simulation.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <omp.h>


const double Simulation::L = 42.0;
const double Simulation::Rcut = 10.0;
const double Simulation::dt = 0.01;
const double Simulation::conversionF = 0.0001 * 4.186;
const double Simulation::mass = 18.0;
const double Simulation::Rcons = 0.00199;
const double Simulation::eStar = 0.2;
const double Simulation::rStar = 3.0;
const double Simulation::gamma = 0.5;



Simulation::Simulation(int nbParticules, int nbPas)
    : N(nbParticules),
      nbPas(nbPas),
      NdL(3 * nbParticules - 3),
      kinEnergy(0.0),
      potEnergy(0.0),
      T_actuel(0.0),
      pdbFile("Fichier_sortie.pdb")
{
    warmUp = static_cast<int>(0.1 * nbPas);
    if (warmUp < 1)
        warmUp = 1;

    pos.resize(N, std::vector<double>(3, 0.0));
    v.resize(N, std::vector<double>(3, 0.0));
    force.resize(N, std::vector<double>(3, 0.0));
    mom.resize(N, std::vector<double>(3, 0.0));

    //stocker le type de chaque particule
    types.resize(N);
}



void Simulation::chargerPositions(const std::string &nomFichier)
{
    std::ifstream fichier(nomFichier);
    if (!fichier.is_open()) {
        std::cerr << "Fichier introuvable : " << nomFichier << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // on lit un entier représentant le type ()A/B)
    int dummy;
    for (int i = 0; i < N; ++i) {
        if (!fichier.good())
            break;
        fichier >> dummy >> pos[i][0] >> pos[i][1] >> pos[i][2];
        types[i] = dummy;
    }
    fichier.close();
}


void Simulation::initialiserMomenta()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            double c = dist(gen);
            double s = dist(gen);
            mom[i][d] = ((s < 0.5) ? 1.0 : -1.0) * c;
        }
    }
}

void Simulation::ajusterTempInit(double Tinit)
{
    double kinInit = 0.0;
    for (int i = 0; i < N; ++i) {
        double px = mom[i][0];
        double py = mom[i][1];
        double pz = mom[i][2];
        kinInit += (px * px + py * py + pz * pz) / mass;
    }
    kinInit /= (2.0 * conversionF);

    double facteur = std::sqrt((NdL * Rcons * Tinit) / kinInit);
    for (int i = 0; i < N; ++i) {
        mom[i][0] *= facteur;
        mom[i][1] *= facteur;
        mom[i][2] *= facteur;
    }
}

void Simulation::enleverDeriveCM()
{
    double pxCM = 0.0, pyCM = 0.0, pzCM = 0.0;
    for (int i = 0; i < N; ++i) {
        pxCM += mom[i][0];
        pyCM += mom[i][1];
        pzCM += mom[i][2];
    }
    pxCM /= N;
    pyCM /= N;
    pzCM /= N;

    for (int i = 0; i < N; ++i) {
        mom[i][0] -= pxCM;
        mom[i][1] -= pyCM;
        mom[i][2] -= pzCM;
    }
}

void Simulation::mettreDansBoite()
{
    for (int i = 0; i < N; ++i) {
        pos[i][0] -= std::floor(pos[i][0] / L) * L;
        pos[i][1] -= std::floor(pos[i][1] / L) * L;
        pos[i][2] -= std::floor(pos[i][2] / L) * L;
    }
}


void Simulation::calculerForcesEtEnergie()
{
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        force[i][0] = force[i][1] = force[i][2] = 0.0;
    }
    potEnergy = 0.0;

    // Définir les décalages de symétrie pour une maille 3×3×3
    std::vector<double> symX, symY, symZ;
    int nSym = 0;
    for (int ix = -1; ix <= 1; ix++) {
        for (int iy = -1; iy <= 1; iy++) {
            for (int iz = -1; iz <= 1; iz++) {
                symX.push_back(ix * L);
                symY.push_back(iy * L);
                symZ.push_back(iz * L);
                nSym++;
            }
        }
    }

    // Boucle sur les images de symétrie parallélisée
    #pragma omp parallel for reduction(+:potEnergy) schedule(dynamic)
    for (int iSym = 0; iSym < nSym; iSym++) {
        double dxSym = symX[iSym];
        double dySym = symY[iSym];
        double dzSym = symZ[iSym];

        for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = pos[i][0] - (pos[j][0] + dxSym);
                double dy = pos[i][1] - (pos[j][1] + dySym);
                double dz = pos[i][2] - (pos[j][2] + dzSym);

                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > Rcut * Rcut)
                    continue;

                double r = std::sqrt(r2);
                if (r < 1e-12)
                    continue; 
                
                //on choisit eLocal, rLocal selon le type des particules 
                double eLocal, rLocal;
                if (types[i] == 2 && types[j] == 2) {
                    eLocal = eStar;
                    rLocal = rStar;
                } else {
                   
                    eLocal = 1.0;
                    rLocal = 3.5;
                }

                // Utilisation de eLocal et rLocal pour le calcul du potentiel
                double s6  = std::pow(rLocal / r, 6);
                double s12 = s6 * s6;
                double uij = 4.0 * eLocal * (s12 - 2.0 * s6);
                potEnergy += uij;


                double fVal = 24.0 * eStar * (2.0 * s12 - s6) / r2;
                double Fxij = fVal * dx;
                double Fyij = fVal * dy;
                double Fzij = fVal * dz;

                // Mises à jour atomiques pour éviter les race conditions
                #pragma omp atomic
                force[i][0] += Fxij;
                #pragma omp atomic
                force[i][1] += Fyij;
                #pragma omp atomic
                force[i][2] += Fzij;

                #pragma omp atomic
                force[j][0] -= Fxij;
                #pragma omp atomic
                force[j][1] -= Fyij;
                #pragma omp atomic
                force[j][2] -= Fzij;
            }
        }
    }
}

// Mise à jour des vitesses 
void Simulation::majVitessesAvant()
{
    for (int i = 0; i < N; ++i) {
        v[i][0] += 0.5 * (force[i][0] / mass) * dt;
        v[i][1] += 0.5 * (force[i][1] / mass) * dt;
        v[i][2] += 0.5 * (force[i][2] / mass) * dt;
    }
}

// Mise à jour des positions
void Simulation::majPositions()
{
    for (int i = 0; i < N; ++i) {
        pos[i][0] += v[i][0] * dt;
        pos[i][1] += v[i][1] * dt;
        pos[i][2] += v[i][2] * dt;
    }
}

// Mise à jour des vitesses 
void Simulation::majVitessesApres()
{
    for (int i = 0; i < N; ++i) {
        v[i][0] += 0.5 * (force[i][0] / mass) * dt;
        v[i][1] += 0.5 * (force[i][1] / mass) * dt;
        v[i][2] += 0.5 * (force[i][2] / mass) * dt;
    }
    // Mise à jour des moment
    for (int i = 0; i < N; ++i) {
        mom[i][0] = mass * v[i][0];
        mom[i][1] = mass * v[i][1];
        mom[i][2] = mass * v[i][2];
    }
}

// Appliquer le thermostat de Berendsen pour ajuster la température
void Simulation::appliquerThermostatBerendsen(double T_cibl)
{
    kinEnergy = 0.0;
    for (int i = 0; i < N; ++i) {
        mom[i][0] = mass * v[i][0];
        mom[i][1] = mass * v[i][1];
        mom[i][2] = mass * v[i][2];
        double px = mom[i][0];
        double py = mom[i][1];
        double pz = mom[i][2];
        kinEnergy += (px * px + py * py + pz * pz) / mass;
    }
    kinEnergy /= (2.0 * conversionF);
    T_actuel = kinEnergy / (NdL * Rcons);

    double lambda = gamma * (T_cibl / T_actuel - 1.0);
    for (int i = 0; i < N; ++i) {
        mom[i][0] *= (1.0 + lambda);
        mom[i][1] *= (1.0 + lambda);
        mom[i][2] *= (1.0 + lambda);
        v[i][0] = mom[i][0] / mass;
        v[i][1] = mom[i][1] / mass;
        v[i][2] = mom[i][2] / mass;
    }
}

// Écrire un cadre PDB pour visualiser la structure
void Simulation::ecrireFramePDB(int step)
{
    static std::ofstream sortie(pdbFile, std::ios::out | std::ios::app);
    if (!sortie.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le PDB : " << pdbFile << std::endl;
        return;
    }

    sortie << "CRYST1"
           << std::right << std::setw(9) << static_cast<int>(L)
           << std::right << std::setw(9) << static_cast<int>(L)
           << std::right << std::setw(9) << static_cast<int>(L)
           << "  90.00  90.00  90.00 P\n";
    sortie << "MODEL" << std::setw(10) << step << "\n";

    for (int i = 0; i < N; ++i) {
        sortie << "ATOM  " << std::setw(6) << i + 1
               << "  Xx          0   "
               << std::fixed << std::setw(8) << std::setprecision(3) << pos[i][0]
               << std::fixed << std::setw(8) << std::setprecision(3) << pos[i][1]
               << std::fixed << std::setw(8) << std::setprecision(3) << pos[i][2]
               << "                  MYRES\n";
    }
    sortie << "TER\nENDMDL\n";
}

// fonction pour calculer la moyenne 

// (Recherche toutes les particules de type A et calcule leurs paires (i<j).
// (Somme l’énergie de Lennard-Jones de chaque paire A–A et obtient un total.
// (Retourne l’énergie moyenne par particule A en divisant par le nombre de particules A.

double Simulation::calculerEnergieMoyenneAA() 
{
    std::vector<int> indicesA; 
    for (int i = 0; i < N; ++i) {
        if (types[i] == 2) {
            indicesA.push_back(i);
        }
    }
    int NA = (int)indicesA.size();

    double totalEnergyAA = 0.0; 
    for (size_t k = 0; k < indicesA.size(); ++k) {
        int i = indicesA[k];
        for (size_t l = k + 1; l < indicesA.size(); ++l) {
            int j = indicesA[l];
            double dx = pos[i][0] - pos[j][0];
            double dy = pos[i][1] - pos[j][1];
            double dz = pos[i][2] - pos[j][2];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            double s6  = std::pow(rStar / r, 6);
            double s12 = s6 * s6;
            double uij = 4.0 * eStar * (s12 - 2.0 * s6);
            totalEnergyAA += uij;
        }
    }

    double energieMoyA = 0.0;
    if (NA > 0) {
        energieMoyA = (2.0 * totalEnergyAA) / NA;
    }
    return energieMoyA;
}

// (1Identifie l’unique particule B (type=1).
// (Calcule l’énergie de Lennard-Jones pour chaque paire (B, A) et cumule le résultat.
// (Retourne l’énergie totale de B par rapport à toutes les particules A.

double Simulation::calculerEnergieBA() 
{
    int idxB = -1;
    for (int i = 0; i < N; ++i) {
        if (types[i] == 1) {
            idxB = i;
            break;
        }
    }
    if (idxB < 0) {
        std::cerr << "Aucune particule B trouvée !" << std::endl;
        return 0.0;
    }
    double totalEnergyBA = 0.0;
    for (int i = 0; i < N; ++i) {
        if (types[i] == 2) {
            double dx = pos[idxB][0] - pos[i][0];
            double dy = pos[idxB][1] - pos[i][1];
            double dz = pos[idxB][2] - pos[i][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            double eLocal = 1.0;
            double rLocal = 3.5;

            double s6  = std::pow(rLocal / r, 6);
            double s12 = s6 * s6;
            double uij = 4.0 * eLocal * (s12 - 2.0 * s6);
            totalEnergyBA += uij;
        }
    }
    return totalEnergyBA;
}



void Simulation::lancerSimulation()
{
    chargerPositions("particule.xyz");
    initialiserMomenta();

    double Tinit = 50.0;
    ajusterTempInit(Tinit);
    enleverDeriveCM();
    ajusterTempInit(Tinit);

    // Initialiser les vitesses
    for (int i = 0; i < N; ++i) {
        v[i][0] = mom[i][0] / mass;
        v[i][1] = mom[i][1] / mass;
        v[i][2] = mom[i][2] / mass;
    }

    // Calcul des forces 
    calculerForcesEtEnergie();

    double T0 = 300.0;
   
    for (int step = 1; step <= nbPas; ++step) {
        // Intégration Velocity–Verlet
        majVitessesAvant();
        majPositions();
        mettreDansBoite();
        calculerForcesEtEnergie();
        majVitessesApres();

        // Régulation de la température (montée en température durant warm‑up)
        double Tcourante = (step <= warmUp) ?
                           (50.0 + (T0 - 50.0) * (static_cast<double>(step) / warmUp)) : T0;
        appliquerThermostatBerendsen(Tcourante);

        if (step % 100 == 0) {
            double pxCM = 0.0, pyCM = 0.0, pzCM = 0.0;
            for (int i = 0; i < N; ++i) {
                pxCM += mom[i][0];
                pyCM += mom[i][1];
                pzCM += mom[i][2];
            }
            double energieTotale = potEnergy + kinEnergy;
            std::cout << "---- Etape " << step << " ----\n";
            std::cout << "Energie potentielle = " << potEnergy << "\n";
            std::cout << "Energie cinetique   = " << kinEnergy << "\n";
            std::cout << "Temperature = " << T_actuel << "\n";
            std::cout << "Energie totale       = " << energieTotale << "\n";
            std::cout << "Momentum CM = " << pxCM << "  " << pyCM << "  " << pzCM << "\n";
            
            std::cout << "---------------------------------------\n";

            ecrireFramePDB(step);
        }
    }
    double eMoyAA = calculerEnergieMoyenneAA();
    double eBA    = calculerEnergieBA(); 
    std::cout << "Energie moyenne d'une particule A vis-a-vis des autres A = " << eMoyAA << std::endl; 
    std::cout << "Energie de la particule B vis-a-vis de toutes les A = "       << eBA    << std::endl;
    std::cout << "----------------------------------------------------\n";
    std::cout << "temperature finale: " << T_actuel << " K\n";
    std::cout << "----------------------------------------------------\n";
}
