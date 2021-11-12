#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include "elemInfo.h"
#include "randSpg.h"
#include "config.h"


using namespace std;
void RandomGeneration(Input *input, Crystal *crystal_list, int begin, int end)
{
    random_device rd;
    mt19937 gen(rd());
    int random_seed = input->GetRandomSeed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> spg(1, 230);

    vector<unsigned int> atoms;
    string tmp_comp = "";
    int z_number = input->GetZNumber();
    vector<string> element = input->GetElement();
    vector<int> composition = input->GetComposition();
    for (int i = 0; i < input->GetNelement(); ++i) {
        tmp_comp += element[i];
        tmp_comp += to_string(composition[i] * z_number);
    }
    ElemInfo::readComposition(tmp_comp, atoms);

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
    vector<pair<unsigned int, double>> manualAtomicRadii;
    double minVolume = initial_volume * 0.9;
    double maxVolume = initial_volume * 1.1;
    vector<pair<unsigned int, char>> forcedWyckAssignments;
    char verbosity = 'n';
    int maxAttempts = 100;
    bool forceMostGeneralWyckPos = false;

    for (int n_crystal = 0; n_crystal < end - begin; ++n_crystal) {
        Crystal crystal;
        do {
            randSpgInput tmp_input(spg(gen), atoms, lmin, lmax,
                                   IADScalingFactor, minRadius,
                                   manualAtomicRadii, minVolume, maxVolume,
                                   forcedWyckAssignments, verbosity,
                                   maxAttempts, forceMostGeneralWyckPos);
            /* Caution: sorted by # of types */
            crystal = RandSpg::randSpgCrystal(tmp_input);
        } while (crystal.getVolume() == 0);
        SortCrystal(&crystal, atoms);
        crystal_list[begin + n_crystal] = crystal; 
    }
}

void Crossover(Input *input, Crystal *parent_list, int n_gene,
               Crystal *crystal_list, int begin, int end)
{
}

void Permutation(Input *input, Crystal *parent_list, int n_gene,
                 Crystal *crystal_list, int begin, int end)
{
    vector<unsigned int> atoms;
    string tmp_comp = "";
    int z_number = input->GetZNumber();
    int nelement = input->GetNelement();
    vector<string> element = input->GetElement();
    vector<int> composition = input->GetComposition();
    for (int i = 0; i < nelement; ++i) {
        tmp_comp += element[i];
        tmp_comp += to_string(composition[i] * z_number);
    }
    ElemInfo::readComposition(tmp_comp, atoms);

    random_device rd;
    mt19937 gen(rd());
    int random_seed = input->GetRandomSeed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> uni(0, n_gene - 1);
    uniform_int_distribution<int> type(0, nelement - 1);

    for (int n_crystal = 0; n_crystal < end - begin; ++n_crystal) {
        Crystal crystal = crystal_list[uni(gen)];
        vector<atomStruct> as = crystal.getAtoms();
        int i;
        int j;
        /* swap element type */
        do {
            i = type(gen); 
            j = type(gen); 
        } while (i == j);

        /* the number of swap atom */
        int max_swap;
        if (composition[i] > composition[j]) {
            max_swap = composition[j] * z_number / 2;
        } else {
            max_swap = composition[i] * z_number / 2;
        }
        uniform_int_distribution<int> swap(0, max_swap);
        int n_swap = swap(gen);

        /* find pair */
        int *acc_index = new int[nelement]();
        for (int k = 1; k < nelement; ++k) {
            acc_index[k] = acc_index[k - 1] + composition[k - 1];
        }
        int *i_index = new int[composition[i]];
        iota(i_index, i_index + composition[i], acc_index[i]);
        shuffle(i_index, i_index + composition[i], gen);
        int *j_index = new int[composition[j]];
        iota(j_index, j_index + composition[j], acc_index[j]);
        shuffle(j_index, j_index + composition[j], gen);

        /* swap */ 
        for (int k = 0; k < n_swap; ++k) {
            iter_swap(as.begin() + i_index[k], as.begin() + j_index[k]); 
        }

        crystal.setAtoms(as);
        SortCrystal(&crystal, atoms);
        crystal_list[begin + n_crystal] = crystal;

        delete []acc_index;
        delete []i_index;
        delete []j_index;
    }
}

void LatticeMutation(Input *input, Crystal *parent_list, int n_gene,
                     Crystal *crystal_list, int begin, int end)
{
    random_device rd;
    mt19937 gen(rd());
    int random_seed = input->GetRandomSeed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> uni(0, n_gene - 1);
    normal_distribution<double> norm(0, 0.1);

    double epsilon;
    for (int n_crystal = 0; n_crystal < end - begin; ++n_crystal) {
        Crystal crystal = parent_list[uni(gen)];
        latticeStruct ls = crystal.getLattice();
        /* a */
        epsilon = norm(gen);
        while ((epsilon > 0.5) || (epsilon < -0.5)) {
            epsilon = norm(gen);
        }
        ls.a *= (1 + epsilon);
        /* b */
        epsilon = norm(gen);
        while ((epsilon > 0.5) || (epsilon < -0.5)) {
            epsilon = norm(gen);
        }
        ls.b *= (1 + epsilon);
        /* c */
        epsilon = norm(gen);
        while ((epsilon > 0.5) || (epsilon < -0.5)) {
            epsilon = norm(gen);
        }
        ls.c *= (1 + epsilon);
        /* alpha */
        epsilon = norm(gen) * 180;
        while ((ls.alpha + epsilon > 180) || (ls.alpha + epsilon < 0)) {
            epsilon = norm(gen);
        }
        ls.alpha += epsilon;
        /* beta */
        epsilon = norm(gen) * 180;
        while ((ls.alpha + epsilon > 180) || (ls.alpha + epsilon < 0)) {
            epsilon = norm(gen);
        }
        ls.beta += epsilon;
        /* gamma */
        epsilon = norm(gen) * 180;
        while ((ls.alpha + epsilon > 180) || (ls.alpha + epsilon < 0)) {
            epsilon = norm(gen);
        }
        ls.gamma += epsilon;

        /* insert crystal */
        crystal_list[begin + n_crystal] = crystal;
    }
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
