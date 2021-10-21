#include <cmath>
#include <random>
#include "config.h"


using namespace std;
vector<Crystal> GenerateCrystal(Input *input)
{
    vector<Crystal> crystal_vector;

    random_device rd;
    mt19937 gen(rd());
    int random_seed = input->GetRandomSeed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> dis(1, 230);

    vector<unsigned int> atoms;
    string comp = "";
    vector<string> element = input->GetElement();
    vector<int> composition = input->GetComposition();
    for (unsigned int i = 0; i < element.size(); ++i) {
        comp += element[i];
        comp += to_string(composition[i] * input->GetZNumber());
    }
    ElemInfo::readComposition(comp, atoms);

    double initial_volume = input->GetVolume();
    double min_length = pow(initial_volume, 1.0/3) / 3;
    double max_length = pow(initial_volume, 1.0/3) * 3;
    latticeStruct lmin = latticeStruct(min_length,
                                       min_length,
                                       min_length,
                                       30.0, 30.0, 30.0);
    latticeStruct lmax = latticeStruct(max_length,
                                       max_length,
                                       max_length,
                                       150.0, 150.0, 150.0);
    double IADScalingFactor = 1.0;
    double minRadius = 0.0;
    // TODO: pairwise constraint
    vector<pair<unsigned int, double>> manualAtomicRadii;
    double minVolume = initial_volume * 0.9;
    double maxVolume = initial_volume * 1.1;
    vector<pair<unsigned int, char>> forcedWyckAssignments;
    char verbosity = 'n';
    int maxAttempts = 100;
    bool forceMostGeneralWyckPos = false;

    int n_population = 0;
    while (n_population < input->GetPopulation()) {
        randSpgInput tmp_input(dis(gen), atoms, lmin, lmax,
                               IADScalingFactor, minRadius, manualAtomicRadii,
                               minVolume, maxVolume, forcedWyckAssignments,
                               verbosity, maxAttempts, forceMostGeneralWyckPos);
        /* Caution: sorted by # of types */
        Crystal crystal = RandSpg::randSpgCrystal(tmp_input);
        if (crystal.getVolume() > 0) {
            SortCrystal(&crystal, atoms);
            crystal_vector.push_back(crystal);
            n_population++;
        }
    }

    return crystal_vector;
}

void SortCrystal(Crystal *crystal, vector<unsigned int> atoms)
{
    vector<atomStruct> new_atoms;
    vector<atomStruct> old_atoms = crystal->getAtoms();
    for (auto i : atoms) {
        for (unsigned int j = 0; j < old_atoms.size(); ++j) {
            if (old_atoms[j].atomicNum == i) {
                new_atoms.push_back(old_atoms[j]);
                old_atoms.erase(old_atoms.begin() + j);
                break;
            } 
        }
    }
    crystal->setAtoms(new_atoms);
}
