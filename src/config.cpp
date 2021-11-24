#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include "elemInfo.h"
#include "randSpg.h"
#include "calculator.h"
#include "config.h"
#include "utils.h"


using namespace std;
void RandomGeneration(Input *input, Crystal *tmp_crystal_list, int n_crystal)
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

    for (int n = 0; n < n_crystal; ++n) {
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
        tmp_crystal_list[n] = crystal; 
    }
}

void Crossover(Input *input, Crystal *crystal_list, int n_gene,
               Crystal *tmp_crystal_list, int n_crystal, MPI_Comm comm)
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

    for (int n = 0; n < n_crystal; ++n) {
        /* select crystals */
        int i;
        int j;
        do {
            i = gene(gen);    
            j = gene(gen);    
        } while (i == j);
        array<double, 2> slab_prob;
        array<int, 2> slab_index;
        array<int, 2> parents = {i, j};
        /* for variable n_grid */
        int n_grid = 10;
        for (int k = 0; k < 2; ++k) {
            double *prob = new double[n_grid * 3];
            vector<atomStruct> as = crystal_list[parents[k]].getAtoms();
            latticeStruct ls = crystal_list[parents[k]].getLattice();
            double lo = 0;
            double hi = lo + 0.5;
            for (int g = 0; g < n_grid; ++g) {
                for (int l = 0; l < 3; ++l) {
                    vector<atomStruct> tmp_as;
                    for (auto v : as) {
                        if ((l == 0) && (v.x >= lo) && (v.x < hi)) {
                            tmp_as.push_back(v);
                        } else if ((l == 1) && (v.y >= lo) && (v.y < hi)) {
                            tmp_as.push_back(v);
                        } else if ((l == 2) && (v.z >= lo) && (v.z < hi)) {
                            tmp_as.push_back(v);
                        } else {
                            continue;
                        }
                    } 
                    Crystal tmp_crystal;
                    tmp_crystal.setAtoms(tmp_as);
                    tmp_crystal.setLattice(ls);
                    /* oneshot for slab by function overloading */
                    double pe = Oneshot(input, &tmp_crystal, l, &comm);
                    prob[g * 3 + l] = exp(-pe / (298 * 8.61733 * 1e-5));
                }
                lo += 0.5 / 10;
                hi += 0.5 / 10;
            }
            double *acc_prob = new double[n_grid * 3];
            acc_prob[0] = prob[0];
            for (int g = 1; g < n_grid * 3; ++g) {
                acc_prob[g] = acc_prob[g - 1] + prob[g];
            }
            uniform_real_distribution<double> real(0, acc_prob[n_grid * 3 - 1]);
            double random = real(gen);
            int index;
            for (index = 0; index < n_grid * 3; ++index) {
                if (random < acc_prob[index]) {
                    break;
                }
            }
            slab_index[k] = index;
            slab_prob[k] = prob[index];
            delete []prob;
            delete []acc_prob;
        }
        //tmp_crystal_list[begin + n_crystal] = tmp_crystal;
    }
}

void Permutation(Input *input, Crystal *crystal_list, int n_gene,
                 Crystal *tmp_crystal_list, int n_crystal)
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

    for (int n = 0; n < n_crystal; ++n) {
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
        tmp_crystal_list[n] = crystal;

        delete []acc_index;
        delete []i_index;
        delete []j_index;
    }
}

void LatticeMutation(Input *input, Crystal *crystal_list, int n_gene,
                     Crystal *tmp_crystal_list, int n_crystal)
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
    for (int n = 0; n < n_crystal; ++n) {
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
        tmp_crystal_list[n] = crystal;
    }
}
