#include <limits>
#include <random>
#include "config.h"


#define MAXLINE 256
using namespace std;
vector<Crystal> generate(Input input)
{
    string comp = "";
    vector<string> element = input.get_element();
    vector<uint> composition = input.get_composition();
    for (unsigned int i = 0; i < element.size(); ++i) {
        comp += element[i];
        comp += to_string(composition[i] * input.get_z_number());
    }

    vector<uint> atoms;
    ElemInfo::readComposition(comp, atoms);

    latticeStruct lmin = latticeStruct(3.0, 3.0, 3.0, 60.0, 60.0, 60.0);
    latticeStruct lmax = latticeStruct(10.0, 10.0, 10.0, 120.0, 120.0, 120.0);

    /* randSpgInput
    spg
    atoms
    latticeMins
    latticeMaxes
    IADScalingFactor
    minRadius
    manualAtomicRadii
    customMinIADs
    minVolume
    maxVolume
    forcedWyckAssignments
    verbosity
    maxAttempts 
    forceMostGeneralWyckPos
    */

    // TODO: Add input parameters
    vector<Crystal> c_vector;
    random_device rd;
    mt19937 gen(rd());
    int random_seed = input.get_random_seed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> dis(1, 230);
    for (int i = 0; i < input.get_population(); ++i) {
        randSpgInput randspg_input(dis(gen), atoms, lmin, lmax);
        /* caution: sort automatically */
        c_vector.push_back(RandSpg::randSpgCrystal(randspg_input));
    }

    return c_vector;
}
