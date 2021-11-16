#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
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

    /* log */
    ofstream log;
    if (rank == 0) {
        log.open("LOG");
        log << "***********************************************************************" << endl;
        log << "*       ______   ______   __   __    __   __    __   ______   ______  *" << endl;
        log << "*      / ____/  / _   /  / /  /  |  / /  /  |  / /  / ____/  / __  /  *" << endl;
        log << "*     / /___   / /_/ /  / /  / / | / /  / / | / /  / /___   / /_/ /   *" << endl;
        log << "*    /___  /  / ____/  / /  / /| |/ /  / /| |/ /  / ____/  /   __/    *" << endl;
        log << "*   ____/ /  / /      / /  / / | / /  / / | / /  / /___   / /\\ \\      *" << endl;
        log << "*  /_____/  /_/      /_/  /_/  |__/  /_/  |__/  /_____/  /_/  \\_\\     *" << endl;
        log << "*                                                                     *" << endl;
        log << "***********************************************************************" << endl;
    }

    /* read input */
    Input *input = ReadInput("./INPUT");

    /* local MPI */
    if ((size < input->GetNpar()) || (size % input->GetNpar() != 0)) {
        if (rank == 0) {
            log << "Check # of cores and NPAR in INPUT!" << endl;
            log.close();
        }
        delete input; 
        exit(1);
    }
    int group_number = size / input->GetNpar();
    int group_index = rank / input->GetNpar();
    int local_rank = rank % input->GetNpar();
    MPI_Comm lammps_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group_index, rank, &lammps_comm);

    /* one-sided communication */
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

    /* global variable */
    Crystal *crystal_list = new Crystal[input->GetMaxPopulation()];
    double *energy_list = new double[input->GetMaxPopulation()];
    double *volume_list = new double[input->GetMaxPopulation()];
    double min_energy = 0;
    int n_atoms = 0;
    int n_best = 0;
    int n_gene = 0;
    int population;
    for (auto i : input->GetComposition()) {
        n_atoms += i * input->GetZNumber();
    }

    if (rank == 0) {
        log << "     Generation        Population           Energy           Volume     " << endl;
    }

    /* main loop */
    for (int gen = 0; gen < input->GetGeneration(); ++gen) {
        MPI_Bcast(&n_best, 1, MPI_INT, 0, MPI_COMM_WORLD);
        population = input->GetPopulation() + n_best;
        if (population < group_number) {
            if (rank == 0) {
                log << "Larger NPAR in INPUT is suggested for faster speed." << endl;
            }
        }
        Crystal *tmp_crystal_list = new Crystal[population];
        if (gen == 0) {
            if (rank == 0) {
                RandomGeneration(input, tmp_crystal_list, 0, population);
            }
        } else {
            if (rank == 0) {
                int tmp_population = input->GetPopulation();
                int disp[5] = {
                    0,
                    (int)round(input->GetRandomGen() * tmp_population),
                    (int)round(input->GetCrossover() * tmp_population),
                    (int)round(input->GetPermutation() * tmp_population),
                    (int)round(input->GetLatticeMut() * tmp_population)
                };
                for (int i = 1; i < 5; ++i) {
                    disp[i] = disp[i] + disp[i - 1];
                }
                int diff = tmp_population - disp[4];
                for (int i = 1; i < 5; ++i) {
                    disp[i] += diff;
                }

                /* random generation */
                if (disp[1] - disp[0]) {
                    RandomGeneration(input, tmp_crystal_list, disp[0], disp[1]);
                }
                /* inheritance */
                if (disp[2] - disp[1]) {
                    Crossover(input, crystal_list, n_gene,
                              tmp_crystal_list, disp[1], disp[2]);
                }
                if (disp[3] - disp[2]) {
                    Permutation(input, crystal_list, n_gene,
                                tmp_crystal_list, disp[2], disp[3]);
                }
                if (disp[4] - disp[3]) {
                    LatticeMutation(input, crystal_list, n_gene,
                                    tmp_crystal_list, disp[3], disp[4]);
                }
                /* best structure */
                for (int i = 0; i < n_best; ++i) {
                    tmp_crystal_list[tmp_population + i] = crystal_list[i];
                }
            }
        }

        /* global information */
        double *global_energy = new double[population];
        int *global_crystal = new int[group_number]; 
        int *global_atom = new int[group_number];
        latticeStruct *global_ls = new latticeStruct[population];
        atomStruct *global_as = new atomStruct[n_atoms * population];

        /* local information*/
        double *local_energy = new double[population];
        int local_crystal = 0;
        int local_atom = 0;
        latticeStruct *local_ls = new latticeStruct[population];
        atomStruct *local_as = new atomStruct[n_atoms * population];

        /* broadcast all crystals from root to headers */
        for (int pop = 0; pop < population; ++pop) {
            if (local_rank == 0) {
                BcastCrystal(
                    tmp_crystal_list,
                    mpi_latticeStruct,
                    mpi_atomStruct,
                    pop,
                    &header_comm
                );
            }
        }

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
            BcastCrystal(
                tmp_crystal_list,
                mpi_latticeStruct,
                mpi_atomStruct,
                local_index,
                &lammps_comm
            );

            Crystal crystal = tmp_crystal_list[local_index];

            /* relax */
            double energy = Relax(input, &crystal, &lammps_comm);
            local_energy[local_crystal] = energy;

            vector<atomStruct> as = crystal.getAtoms();
            for (int i = 0; i < n_atoms; ++i) {
                local_as[local_atom + i] = as[i];
            }

            local_ls[local_crystal] = crystal.getLattice();

            /* update local_crystal */
            local_crystal++;
            local_atom += n_atoms;
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
            if (group_number > 1) {
                for (int i = 1; i < group_number; ++i) {
                    disp[i] = disp[i - 1] + global_crystal[i - 1];
                }
            }

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
            for (int i = 0; i < population; ++i) {
                vector<atomStruct> as;
                as.reserve(n_atoms);
                for (int j = 0; j < n_atoms; ++j) {
                    as.push_back(global_as[i * n_atoms + j]);
                }
                tmp_crystal_list[i].setAtoms(as);
                tmp_crystal_list[i].setLattice(global_ls[i]);
            }

            /* sort by energy */
            int *argsort = new int[population];
            iota(argsort, argsort + population, 0);
            auto compare = [&global_energy](int i, int j)
                           {return global_energy[i] < global_energy[j];};
            stable_sort(argsort, argsort + population, compare);
            min_energy = global_energy[argsort[0]];

            /* copy, count and log */
            n_best = 0;
            n_gene = 0;
            char tmp_log[128]; 
            string tmp_string;
            for (int i = 0; i < population; ++i) {
                sprintf(tmp_log,
                        "%12d%18d%20.3f%17.3f",
                        gen + 1,
                        i + 1,
                        global_energy[i],
                        tmp_crystal_list[i].getVolume());
                tmp_string = tmp_log;
                log << tmp_string << endl;
                crystal_list[i] = tmp_crystal_list[argsort[i]];
                energy_list[i] = global_energy[argsort[i]];
                volume_list[i] = crystal_list[i].getVolume();
                if (energy_list[i] - min_energy >
                    input->GetGeneWindow() * n_atoms) {
                    n_gene++;
                    if (energy_list[i] - min_energy >
                        input->GetBestWindow() * n_atoms) {
                        n_best++;
                    }
                }
            }
            delete []argsort;
        }

        /* initialize global_index */
        *global_index = 0;

        delete []tmp_crystal_list;

        delete []global_energy;
        delete []global_crystal;
        delete []global_atom;
        delete []global_ls;
        delete []global_as;

        delete []local_energy;
        delete []local_ls;
        delete []local_as;
    }

    delete []crystal_list;
    delete []energy_list;
    delete []volume_list;
    delete input;

    /* close log */
    if (rank == 0) {
        log.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_free(&win);
    MPI_Comm_free(&lammps_comm);
    MPI_Comm_free(&header_comm);
    MPI_Finalize();
}
