#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include "elemInfo.h"
#include "randSpg.h"
#include "config.h"
#include "utils.h"


using namespace std;
void RandomGeneration(Input *input, Crystal *tmp_crystal_list, int begin, int end)
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
        } while (!CheckLattice(crystal.getLattice()));
        SortCrystal(&crystal, atoms);
        tmp_crystal_list[begin + n_crystal] = crystal; 
    }
}

void Crossover(Input *input, Crystal *crystal_list, int n_gene,
               Crystal *tmp_crystal_list, int begin, int end)
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
    uniform_int_distribution<int> gene(0, n_gene - 1);

    for (int n_crystal = 0; n_crystal < end - begin; ++n_crystal) {
        /* select old_crystals */
        int i;
        int j;
        do {
            i = gene(gen);    
            j = gene(gen);    
        } while (i == j);
        Crystal crystal_i = crystal_list[i];
        Crystal crystal_j = crystal_list[j];
        /* make slab */
        for (int k = 0; k < 10; ++k) {
            for (int l = 0; l < 3; ++l) {
                /* copy crystal */
            }
        }
    }
}

void Permutation(Input *input, Crystal *crystal_list, int n_gene,
                 Crystal *tmp_crystal_list, int begin, int end)
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
    uniform_int_distribution<int> gene(0, n_gene - 1);
    uniform_int_distribution<int> type(0, nelement - 1);

    for (int n_crystal = 0; n_crystal < end - begin; ++n_crystal) {
        Crystal crystal = crystal_list[gene(gen)];
        vector<atomStruct> as = crystal.getAtoms();
        /* select element type */
        int i;
        int j;
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
        if (max_swap < 1) {
            max_swap = 1;
        }
        uniform_int_distribution<int> swap(1, max_swap);
        int n_swap = swap(gen);

        /* find pair */
        int *acc_index = new int[nelement]();
        for (int k = 1; k < nelement; ++k) {
            acc_index[k] = acc_index[k - 1] + composition[k - 1] * z_number;
        }

        int *i_index = new int[composition[i] * z_number];
        iota(i_index, i_index + composition[i] * z_number, acc_index[i]);
        shuffle(i_index, i_index + composition[i] * z_number, gen);
        int *j_index = new int[composition[j] * z_number];
        iota(j_index, j_index + composition[j] * z_number, acc_index[j]);
        shuffle(j_index, j_index + composition[j] * z_number, gen);

        /* swap */ 
        for (int k = 0; k < n_swap; ++k) {
            iter_swap(as.begin() + i_index[k], as.begin() + j_index[k]); 
        }

        crystal.setAtoms(as);
        SortCrystal(&crystal, atoms);
        tmp_crystal_list[begin + n_crystal] = crystal;

        delete []acc_index;
        delete []i_index;
        delete []j_index;
    }
}

void LatticeMutation(Input *input, Crystal *crystal_list, int n_gene,
                     Crystal *tmp_crystal_list, int begin, int end)
{
    random_device rd;
    mt19937 gen(rd());
    int random_seed = input->GetRandomSeed();
    if (random_seed != -1) {
        gen.seed(random_seed);
    }
    uniform_int_distribution<int> gene(0, n_gene - 1);
    normal_distribution<double> norm(0, 0.1);

    double epsilon;
    for (int n_crystal = 0; n_crystal < end - begin; ++n_crystal) {
        Crystal crystal = crystal_list[gene(gen)];
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
        double alpha = ls.alpha;
        double beta = ls.beta;
        double gamma = ls.gamma;
        while (1) {
            /* alpha */
            epsilon = norm(gen) * 180;
            while ((alpha + epsilon > 180) || (alpha + epsilon < 0)) {
                epsilon = norm(gen);
            }
            ls.alpha = alpha + epsilon;
            /* beta */
            epsilon = norm(gen) * 180;
            while ((beta + epsilon > 180) || (beta + epsilon < 0)) {
                epsilon = norm(gen);
            }
            ls.beta = beta + epsilon;
            /* gamma */
            epsilon = norm(gen) * 180;
            while ((gamma + epsilon > 180) || (gamma + epsilon < 0)) {
                epsilon = norm(gen);
            }
            ls.gamma = gamma + epsilon;
            if (CheckLattice(ls)) {
                crystal.setLattice(ls);
                break;
            } 
        }
        /* insert crystal */
        vector<atomStruct> as = crystal.getAtoms();
        tmp_crystal_list[begin + n_crystal] = crystal;
    }
}
