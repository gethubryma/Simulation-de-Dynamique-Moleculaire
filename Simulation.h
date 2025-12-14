#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <vector>


class Simulation {

    
private:
    
    int N;          // Nombre de particules
    int nbPas;      // Nombre de pas de simulation
    int NdL;        // (3*N - 3)
    int warmUp;     // Nombre de pas pour la montée en température

    // Constantes
    static const double L;             // Taille de la boîte
    static const double Rcut;          // Rayon de coupure
    static const double dt;            // Pas de temps
    static const double conversionF;   
    static const double mass;          
    static const double Rcons;         
    static const double eStar;         
    static const double rStar;         
    static const double gamma;     
    


    // type pour chaque particule 
    std::vector<int> types;

    double kinEnergy;
    double potEnergy;
    double T_actuel;

    std::vector<std::vector<double>> pos;
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> force;
    std::vector<std::vector<double>> mom;
    std::string pdbFile;
public:
    
    Simulation(int nbParticules, int nbPas);

    
    void lancerSimulation();
    //pour le calcul de la moyenne et l'energie 
    double calculerEnergieMoyenneAA();
    double calculerEnergieBA();

private:
    
    void chargerPositions(const std::string &nomFichier);
    void initialiserMomenta();         
    void ajusterTempInit(double Tinit);
    void enleverDeriveCM();
    void mettreDansBoite();
    void calculerForcesEtEnergie();   
    void majVitessesAvant();
    void majPositions();
    void majVitessesApres();
    void appliquerThermostatBerendsen(double Tcible);
    void ecrireFramePDB(int step);

};

#endif // SIMULATION_H
