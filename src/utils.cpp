#include <algorithm>
#include "utils.h"

using namespace std;
void SortCrystal(Crystal *crystal, vector<unsigned int> atoms)
{
    vector<atomStruct> as = crystal->getAtoms();
    auto compare = [&atoms](atomStruct i, atomStruct j)
                    {return find(atoms.begin(), atoms.end(), i.atomicNum) <
                            find(atoms.begin(), atoms.end(), j.atomicNum);};
    stable_sort(as.begin(), as.end(), compare);
    crystal->setAtoms(as);
}


void WriteCrystal(Crystal *crystal, string filename)
{
}


void GathervCrystal(Crystal *local_crystal_list,
                    Crystal *global_crystal_list,
                    MPI_Datatype mpi_latticeStruct,
                    MPI_Datatype mpi_atomStruct,
                    int local_crystal,
                    MPI_Comm comm)
{
    int rank;
    int size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int n_atoms;
    if (rank == 0) {
        n_atoms = (int)local_crystal_list[0].numAtoms();
    }
    MPI_Bcast(&n_atoms, 1, MPI_INT, 0, comm);

    int *global_crystal = new int[size];
    MPI_Gather(&local_crystal, 1, MPI_INT,
               global_crystal, 1, MPI_INT, 0, comm);

    int local_atom = local_crystal * n_atoms;
    int *global_atom = new int[size];
    MPI_Gather(&local_atom, 1, MPI_INT,
               global_atom, 1, MPI_INT, 0, comm);

    int total_crystal = 0;
    for (int i = 0; i < size; ++i) {
        total_crystal += global_crystal[i];
    }
    MPI_Bcast(&total_crystal, 1, MPI_INT, 0, comm);

    /* split Crystal into latticeStruct and atomStruct */
    latticeStruct *local_ls = new latticeStruct[local_crystal];
    atomStruct *local_as = new atomStruct[n_atoms * local_crystal];
    for (int i = 0; i < local_crystal; ++i) {
        local_ls[i] = local_crystal_list[i].getLattice();
        vector<atomStruct> as = local_crystal_list[i].getAtoms();
        for (int j = 0; j < n_atoms; ++j) {
            local_as[i * n_atoms + j] = as[j];
        }
    }

    /* gatherv latticeStruct */
    int *disp = new int[size]();
    for (int i = 1; i < size; ++i) {
        disp[i] = disp[i - 1] + global_crystal[i - 1]; 
    }
    latticeStruct *global_ls = new latticeStruct[total_crystal];
    MPI_Gatherv(local_ls, local_crystal, mpi_latticeStruct,
                global_ls, global_crystal, disp, mpi_latticeStruct,
                0, comm);

    /* gatherv atomStruct */
    for (int i = 1; i < size; ++i) {
        disp[i] = disp[i - 1] + global_atom[i - 1];
    }
    atomStruct *global_as = new atomStruct[total_crystal * n_atoms];
    MPI_Gatherv(local_as, local_atom, mpi_atomStruct,
                global_as, global_atom, disp, mpi_atomStruct,
                0, comm);

    /* construct Crystal */
    if (rank == 0) {
        for (int i = 0; i < total_crystal; ++i) {
            global_crystal_list[i].setLattice(global_ls[i]);
            vector<atomStruct> as;
            as.reserve(n_atoms);
            for (int j = 0; j < n_atoms; ++j) {
                as.push_back(global_as[i * n_atoms + j]);
            }
            global_crystal_list[i].setAtoms(as);
        }
    }

    delete []disp;
    delete []global_crystal;
    delete []local_ls;
    delete []local_as;
    delete []global_ls;
    delete []global_as;
}


void BcastCrystal(Crystal *crystal_list,
                  MPI_Datatype mpi_latticeStruct,
                  MPI_Datatype mpi_atomStruct,
                  int index,
                  MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    latticeStruct ls;
    if (rank == 0) {
        ls = crystal_list[index].getLattice();
    }
    MPI_Bcast(&ls, 1, mpi_latticeStruct, 0, comm);

    int n_atoms;
    if (rank == 0) {
        n_atoms = (int)crystal_list[index].numAtoms();
    }
    MPI_Bcast(&n_atoms, 1, MPI_INT, 0, comm);

    vector<atomStruct> as;
    if (rank == 0) {
        as = crystal_list[index].getAtoms();
    } else {
        as.resize(n_atoms);
    }
    MPI_Bcast(&as[0], n_atoms, mpi_atomStruct, 0, comm);

    if (rank != 0) {
        crystal_list[index] = Crystal(ls, as);
    }
}

double GetMassFromSymbol(string symbol)
{
    const vector<pair<string, double>> mass =
    {
        {"H", 1.008},
        {"He", 4.003},
        {"Li", 6.941},
        {"Be", 9.012},
        {"B", 10.811},
        {"C", 12.011},
        {"N", 14.007},
        {"O", 15.999},
        {"F", 18.998},
        {"Ne", 20.180},
        {"Na", 22.990},
        {"Mg", 24.305},
        {"Al", 26.982},
        {"Si", 28.086},
        {"P", 30.974},
        {"S", 32.065},
        {"Cl", 35.453},
        {"Ar", 39.948},
        {"K", 39.098},
        {"Ca", 40.078},
        {"Sc", 44.956},
        {"Ti", 47.867},
        {"V", 50.942},
        {"Cr", 51.996},
        {"Mn", 54.938},
        {"Fe", 55.845},
        {"Co", 58.933},
        {"Ni", 58.693},
        {"Cu", 63.546},
        {"Zn", 65.380},
        {"Ga", 69.723},
        {"Ge", 72.640},
        {"As", 74.922},
        {"Se", 78.960},
        {"Br", 79.904},
        {"Kr", 83.798},
        {"Rb", 85.468},
        {"Sr", 87.620},
        {"Y", 88.906},
        {"Zr", 91.224},
        {"Nb", 92.906},
        {"Mo", 95.960},
        {"Tc", 98.000},
        {"Ru", 101.07},
        {"Rh", 102.91},
        {"Pd", 106.42},
        {"Ag", 107.87},
        {"Cd", 112.41},
        {"In", 114.82},
        {"Sn", 118.71},
        {"Sb", 121.76},
        {"Te", 127.60},
        {"I", 126.90},
        {"Xe", 131.29},
        {"Cs", 132.91},
        {"Ba", 137.33},
        {"La", 138.91},
        {"Ce", 140.12},
        {"Pr", 140.91},
        {"Nd", 144.24},
        {"Pm", 145.00},
        {"Sm", 150.36},
        {"Eu", 151.96},
        {"Gd", 157.25},
        {"Tb", 158.93},
        {"Dy", 162.50},
        {"Ho", 164.93},
        {"Er", 167.26},
        {"Tm", 168.93},
        {"Yb", 173.05},
        {"Lu", 174.97},
        {"Hf", 178.49},
        {"Ta", 180.95},
        {"W", 183.84},
        {"Re", 186.21},
        {"Os", 190.23},
        {"Ir", 192.22},
        {"Pt", 195.08},
        {"Au", 196.97},
        {"Hg", 200.59},
        {"Tl", 204.38},
        {"Pb", 207.20},
        {"Bi", 208.98},
        {"Po", 209.00},
        {"At", 210.00},
        {"Rn", 222.00},
        {"Fr", 223.00},
        {"Ra", 226.00},
        {"Ac", 227.00},
        {"Th", 232.04},
        {"Pa", 231.04},
        {"U", 238.03}
    };
    for (auto i : mass) {
        if (i.first == symbol) {
            return i.second;
        }
    }
    return -1;
}
