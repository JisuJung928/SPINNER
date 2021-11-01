#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <mpi.h>

#include "calculator.h"
#include "config.h"
#include "input.h"
#include "utils.h"

using namespace std;
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* latticeStruct type for MPI */
    MPI_Datatype ls_type[6] =
    {
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE
    };
    int ls_block[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint ls_disp[6];

    ls_disp[0] = offsetof(latticeStruct, a);
    ls_disp[1] = offsetof(latticeStruct, b);
    ls_disp[2] = offsetof(latticeStruct, c);
    ls_disp[3] = offsetof(latticeStruct, alpha);
    ls_disp[4] = offsetof(latticeStruct, beta);
    ls_disp[5] = offsetof(latticeStruct, gamma);

    MPI_Datatype mpi_latticeStruct;
    MPI_Type_create_struct(6, ls_block, ls_disp, ls_type, &mpi_latticeStruct);
    MPI_Type_commit(&mpi_latticeStruct);

    /* atomStruct type for MPI */
    MPI_Datatype as_type[4] =
    {
        MPI_UNSIGNED,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE
    };
    int as_block[4] = {1, 1, 1, 1};
    MPI_Aint as_disp[4];

    as_disp[0] = offsetof(atomStruct, atomicNum);
    as_disp[1] = offsetof(atomStruct, x);
    as_disp[2] = offsetof(atomStruct, y);
    as_disp[3] = offsetof(atomStruct, z);

    MPI_Datatype mpi_atomStruct;
    MPI_Type_create_struct(4, as_block, as_disp, as_type, &mpi_atomStruct);
    MPI_Type_commit(&mpi_atomStruct);

    /* read input */
    Input *input = ReadInput("./INPUT");

    /* local MPI */
    int group_number = size / input->GetNpar();
    int group_index = rank / input->GetNpar();
    int local_rank = rank % input->GetNpar();
    MPI_Comm lammps_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group_index, rank, &lammps_comm);

    /* One-sided communication */
    int header;
    if (local_rank == 0) {
        header = 1;
    } else {
        header = 0;
    }
    MPI_Comm header_comm;
    MPI_Comm_split(MPI_COMM_WORLD, header, rank, &header_comm);

    MPI_Win win;
    int *global_index;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_index, &win);
    *global_index = 0;

    /* generate crystal */
    for (int gen = 0; gen < input->GetGeneration(); ++gen) {
        // TODO: population < group_number?
        int population = input->GetPopulation();
        Crystal *crystal_list = new Crystal[population];
        if (gen == 0) {
            if (rank == 0) {
                RandomGeneration(input, crystal_list, 0, population);
            }
        } else {
            // TODO: evolution
            if (rank == 0) {
                //double random_gen = input->GetRandomGen();
                //int random_num = (int)(population * random_gen);
                RandomGeneration(input, crystal_list, 0, population);
            }
        }

        int total_n_atoms = 0;
        for (auto i : input->GetComposition()) {
            total_n_atoms += i * input->GetZNumber();
        }
        total_n_atoms *= population;

        /* global information */
        double *global_energy = new double[population];
        int *global_n_atoms = new int[population];
        int *global_crystal = new int[group_number]; 
        int *global_atom = new int[group_number];
        latticeStruct *global_ls = new latticeStruct[population];
        atomStruct *global_as = new atomStruct[total_n_atoms];

        /* local information*/
        double *local_energy = new double[population];
        int *local_n_atoms = new int[population];
        int local_crystal = 0;
        int local_atom = 0;
        latticeStruct *local_ls = new latticeStruct[population];
        atomStruct *local_as = new atomStruct[total_n_atoms];

        /* broadcast all crystals from root to headers */
        for (int pop = 0; pop < population; ++pop) {
            if (local_rank == 0) {
                DistributeCrystal(
                    crystal_list,
                    mpi_latticeStruct,
                    mpi_atomStruct,
                    pop,
                    &header_comm
                );
            }
        }

        int n_atoms_acc = 0;
        int n_atoms;
        int local_index;
        int one = 1;
        while (1) {
            /* passive syncronization */
            if (local_rank == 0) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
                MPI_Fetch_and_op(&one, &local_index, MPI_INT, 0, (MPI_Aint)0,
                                 MPI_SUM, win);
                MPI_Win_unlock(0, win);
            }

            /* Caution: local rank */
            MPI_Bcast(&local_index, 1, MPI_INT, 0, lammps_comm);
            if (local_index >= population) {
                break;
            }

            /* broadcast the crystal from root of lammps_comm to others */
            DistributeCrystal(
                crystal_list,
                mpi_latticeStruct,
                mpi_atomStruct,
                local_index,
                &lammps_comm
            );
            Crystal crystal = crystal_list[local_index];
            n_atoms = (int)crystal.numAtoms();
            local_n_atoms[local_crystal] = n_atoms;
            local_atom += n_atoms;

            /* relax */
            double energy = Relax(input, &crystal, &lammps_comm);
            local_energy[local_crystal] = energy;

            vector<atomStruct> as = crystal.getAtoms();
            int n_atoms = (int)as.size();
            for (int i = 0; i < n_atoms; ++i) {
                local_as[n_atoms_acc + i] = as[i];
            }
            n_atoms_acc += n_atoms;

            local_ls[local_crystal] = crystal.getLattice();

            /* update local_crystal */
            local_crystal++;
        }

        /* gather */
        if (local_rank == 0) {
            /* send the number of crystals in each node */
            MPI_Allgather(&local_crystal, 1, MPI_INT, 
                          global_crystal, 1, MPI_INT, header_comm);

            /* send the number of atoms in each node  */
            MPI_Gather(&local_atom, 1, MPI_INT,
                       global_atom, 1, MPI_INT, 0, header_comm);

            /* send the number of atoms in each crystal */
            int *disp = new int[group_number]();
            for (int i = 1; i < group_number; ++i) {
                disp[i] = disp[i - 1] + global_crystal[i - 1];
            }

            MPI_Allgatherv(local_n_atoms, local_crystal, MPI_INT,
                           global_n_atoms, global_crystal, disp, MPI_INT,
                           header_comm);

            /* send energy vector */
            MPI_Gatherv(local_energy, local_crystal, MPI_DOUBLE,
                        global_energy, global_crystal, disp, MPI_DOUBLE,
                        0, header_comm);

            /* send latticeStruct vector */
            MPI_Gatherv(local_ls, local_crystal, mpi_latticeStruct,
                        global_ls, global_crystal, disp, mpi_latticeStruct,
                        0, header_comm);

            /* send atomStruct vector */
            for (int i = 1; i < group_number; ++i) {
                disp[i] = disp[i - 1] + global_atom[i - 1];
            }
            MPI_Gatherv(local_as, local_atom, mpi_atomStruct,
                        global_as, global_atom, disp, mpi_atomStruct,
                        0, header_comm);

            delete []disp;
        }

        if (rank == 0) {
            n_atoms_acc = 0;
            for (int i = 0; i < population; ++i) {
                vector<atomStruct> as;
                as.reserve(global_n_atoms[i]);
                for (int j = 0; j < global_n_atoms[i]; ++j) {
                    as[j] = global_as[n_atoms_acc + j];
                }
                n_atoms_acc += global_n_atoms[i];
                crystal_list[i].setAtoms(as);
                crystal_list[i].setLattice(global_ls[i]);
            }
            /* sort by energy */
            int *index = new int[population];
            iota(index, index + population, 0);
            auto compare = [&global_energy](int i, int j)
                           {return global_energy[i] < global_energy[j];};
            stable_sort(index, index + population, compare);
            delete []index;
        }

        /* initialize global_index */
        *global_index = 0;

        delete []crystal_list;

        delete []global_energy;
        delete []global_n_atoms;
        delete []global_crystal;
        delete []global_atom;
        delete []global_ls;
        delete []global_as;

        delete []local_energy;
        delete []local_n_atoms;
        delete []local_ls;
        delete []local_as;
    }

    delete input;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_free(&win);
    MPI_Comm_free(&lammps_comm);
    MPI_Comm_free(&header_comm);
    MPI_Finalize();
}
