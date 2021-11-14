#include "utils.h"

using namespace std;
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


void WriteCrystal(Crystal *crystal, string filename)
{
}


void GatherCrystal(Crystal *crystal_list,
                   MPI_Datatype mpi_latticeStruct,
                   MPI_Datatype mpi_atomStruct,
                   int index,
                   MPI_Comm *comm)
{
}


void BcastCrystal(Crystal *crystal_list,
                  MPI_Datatype mpi_latticeStruct,
                  MPI_Datatype mpi_atomStruct,
                  int index,
                  MPI_Comm *comm)
{
    int rank;
    MPI_Comm_rank(*comm, &rank);

    latticeStruct ls;
    if (rank == 0) {
        ls = crystal_list[index].getLattice();
    }
    MPI_Bcast(&ls, 1, mpi_latticeStruct, 0, *comm);

    int n_atoms;
    if (rank == 0) {
        n_atoms = (int)crystal_list[index].numAtoms();
    }
    MPI_Bcast(&n_atoms, 1, MPI_INT, 0, *comm);

    vector<atomStruct> as;
    if (rank == 0) {
        as = crystal_list[index].getAtoms();
    } else {
        as.resize(n_atoms);
    }
    MPI_Bcast(&as[0], n_atoms, mpi_atomStruct, 0, *comm);

    if (rank != 0) {
        crystal_list[index] = Crystal(ls, as);
    }
}
