#include "Simulation.h"
#include <iostream>

int main()
{
    int nbParticules, nbPas;
    std::cout << "Entrez le nombre de particules : ";
    std::cin >> nbParticules;

    if (nbParticules > 1000) {
        std::cerr << "Erreur :  " << std::endl;
        return 1;
    }

    std::cout << "Entrez le nombre de pas : ";
    std::cin >> nbPas;

    // Création et lancement de la simulation
    Simulation sim(nbParticules, nbPas);
    sim.lancerSimulation();

    return 0;
}
