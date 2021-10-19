#include <random>
#include "config.h"

using namespace std;
vector<Crystal> GenerateCrystal(Input *input)
{
    string comp = "";
    vector<string> element = input->GetElement();
    vector<int> composition = input->GetComposition();
    for (unsigned int i = 0; i < element.size(); ++i) {
        comp += element[i];
        comp += to_string(composition[i] * input->GetZNumber());
    }

    vector<uint> atoms;
    ElemInfo::readComposition(comp, atoms);

    latticeStruct lmin = latticeStruct(3.0, 3.0, 3.0, 60.0, 60.0, 60.0);
    latticeStruct lmax = latticeStruct(10.0, 10.0, 10.0, 120.0, 120.0, 120.0);

    /* randSpgInput
     * spg
     * atoms
     * latticeMins
     * latticeMaxes
     * IADScalingFactor
     * minRadius
     * manualAtomicRadii
     * customMinIADs
     * minVolume
     * maxVolume
     * forcedWyckAssignments
     * verbosity
     * maxAttempts 
     * forceMostGeneralWyckPos
     */

    // TODO: Add input parameters
    vector<Crystal> crystal_vector;
    random_device rd;
    mt19937 gen(rd());
    int random_seed = input->GetRandomSeed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> dis(1, 230);
    int num_attempt = 0;
    while (num_attempt < input->GetPopulation()) {
        randSpgInput randspg_input(dis(gen), atoms, lmin, lmax);
        /* Caution: sorted by # of types */
        Crystal crystal = RandSpg::randSpgCrystal(randspg_input);
        if (crystal.getVolume() > 0) {
            SortCrystal(&crystal, atoms);
            crystal_vector.push_back(crystal);
            num_attempt++;
        }
    }

    return crystal_vector;
}

void SortCrystal(Crystal *crystal, vector<uint> atoms)
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
