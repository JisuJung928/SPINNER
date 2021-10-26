#include <cmath>
#define _USE_MATH_DEFINES
#include <iostream>

#define LAMMPS_LIB_MPI
#include "library.h"
#include "calculator.h"
#include "crystal.h"


void *LammpsInit(Input *input, Crystal *crystal, MPI_Comm *comm, char *filename)
{
    char cmd[256];
    char *lmpargv[] = {(char *)"liblammps", (char *)"-screen", (char *)"none",
                       (char *)"-log", filename};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    void *lmp = nullptr;
    lmp = lammps_open(lmpargc, lmpargv, *comm, nullptr);

    /* basic */
    const char *cmds[] = {"units metal",
                          "box tilt large",
                          "atom_modify map array sort 0 0.0",
                          "neigh_modify every 1 delay 0 check yes"};
    lammps_commands_list(lmp, sizeof(cmds) / sizeof(const char *), cmds);

    /* box */
    latticeStruct lattice = crystal->getLattice();
    double cell[3][3];
    cell[0][0] = lattice.a;
    cell[0][1] = 0;
    cell[0][2] = 0;
    cell[1][0] = lattice.b * cos(lattice.gamma * M_PI / 180);
    cell[1][1] = lattice.b * sin(lattice.gamma * M_PI / 180);
    cell[1][2] = 0;
    cell[2][0] = lattice.c * cos(lattice.beta * M_PI / 180);
    cell[2][1] = (lattice.b * lattice.c * cos(lattice.alpha * M_PI / 180)
                 - cell[1][0] * cell[2][0]) / cell[1][1];
    cell[2][2] = sqrt(lattice.c * lattice.c
                    - cell[2][0] * cell[2][0]
                    - cell[2][1] * cell[2][1]);

    sprintf(cmd, "region cell prism 0 %f 0 %f 0 %f %f %f %f units box",
            cell[0][0], cell[1][1], cell[2][2],
            cell[1][0], cell[2][0], cell[2][1]);
    lammps_command(lmp, cmd);
    sprintf(cmd, "create_box %d cell", input->GetNelement());
    lammps_command(lmp, cmd);

    /* atoms */
    int *id = new int[crystal->numAtoms()];
    int *type = new int[crystal->numAtoms()];
    double *position = new double[3 * crystal->numAtoms()];

    int z_number = input->GetZNumber();
    int nelement = input->GetNelement();
    vector<int> composition = input->GetComposition();
    vector<atomStruct> atoms = crystal->getAtoms();
    
    int i = 0;
    for (int j = 0; j < nelement; ++j) {
        for (int k = 0; k < composition[j] * z_number; ++k) {
            id[i] = i + 1;
            type[i] = j + 1;
            position[3 * i + 0] = atoms[i].x * cell[0][0]
                                + atoms[i].y * cell[1][0]
                                + atoms[i].z * cell[2][0];
            position[3 * i + 1] = atoms[i].y * cell[1][1]
                                + atoms[i].z * cell[2][1];
            position[3 * i + 2] = atoms[i].z * cell[2][2];
            i++; 
        }
    }
    lammps_create_atoms(lmp, crystal->numAtoms(), id, type, position,
                        nullptr, nullptr, 0);

    /* mass */
    vector<double> mass = input->GetMass();
    for (i = 0; i < nelement; ++i) {
        sprintf(cmd, "mass %d %f", i + 1, mass[i]);
        lammps_command(lmp, cmd);
    }

    /* potential */
    sprintf(cmd, "pair_style %s", input->GetPairStyle());
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->GetPairCoeff());
    lammps_command(lmp, cmd);

    delete []id;
    delete []type;
    delete []position;

    return lmp;
}


double Relax(Input *input, Crystal *crystal, MPI_Comm *comm, int tag)
{
    char cmd[256];

    /* create LAMMPS instance */
    char filename[256];
    sprintf(filename, "log.lammps_%d", tag);
    void *lmp = LammpsInit(input, crystal, comm, filename);

    /* minimize */
    lammps_command(lmp, "fix int all box/relax tri 0.0");
    sprintf(cmd, "minimize 0 %f 1000 10000", input->GetMaxForce());
    lammps_command(lmp, cmd);

    /* update positions */
    double *position = new double[3 * crystal->numAtoms()];
    lammps_gather_atoms(lmp, (char *)"x", 2, 3, position);

    double *boxlo = new double[3];
    double *boxhi = new double[3];
    double xy;
    double yz;
    double xz;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, nullptr, nullptr);

    double cell[3][3];
    cell[0][0] = boxhi[0] - boxlo[0];
    cell[0][1] = 0;
    cell[0][2] = 0;
    cell[1][0] = xy;
    cell[1][1] = boxhi[1] - boxlo[1];
    cell[1][2] = 0;
    cell[2][0] = xz;
    cell[2][1] = yz;
    cell[2][2] = boxhi[2] - boxlo[2];

    double cross[3][3];
    cross[0][0] = cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1];
    cross[0][1] = cell[1][2] * cell[2][0] - cell[1][0] * cell[2][2];
    cross[0][2] = cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0];
    cross[1][0] = cell[2][1] * cell[0][2] - cell[2][2] * cell[0][1];
    cross[1][1] = cell[2][2] * cell[0][0] - cell[2][0] * cell[0][2];
    cross[1][2] = cell[2][0] * cell[0][1] - cell[2][1] * cell[0][0];
    cross[2][0] = cell[0][1] * cell[1][2] - cell[0][2] * cell[1][1];
    cross[2][1] = cell[0][2] * cell[1][0] - cell[0][0] * cell[1][2];
    cross[2][2] = cell[0][0] * cell[1][1] - cell[0][1] * cell[1][0];
    double vol = cross[0][0] * cell[0][0]
               + cross[0][1] * cell[0][1]
               + cross[0][2] * cell[0][2];

    double inv[3][3];
    inv[0][0] = cross[0][0] / vol;
    inv[0][1] = cross[1][0] / vol;
    inv[0][2] = cross[2][0] / vol;
    inv[1][0] = cross[0][1] / vol;
    inv[1][1] = cross[1][1] / vol;
    inv[1][2] = cross[2][1] / vol;
    inv[2][0] = cross[0][2] / vol;
    inv[2][1] = cross[1][2] / vol;
    inv[2][2] = cross[2][2] / vol;

    vector<atomStruct> as = crystal->getAtoms();
    for (unsigned int i = 0; i < crystal->numAtoms(); ++i) {
        as[i].x = position[3 * i + 0] * inv[0][0]
                + position[3 * i + 1] * inv[1][0]
                + position[3 * i + 2] * inv[2][0];
        as[i].y = position[3 * i + 0] * inv[0][1]
                + position[3 * i + 1] * inv[1][1]
                + position[3 * i + 2] * inv[2][1];
        as[i].z = position[3 * i + 0] * inv[0][2]
                + position[3 * i + 1] * inv[1][2]
                + position[3 * i + 2] * inv[2][2];
    }
    crystal->setAtoms(as);

    /* update lattice */
    double a = lammps_get_thermo(lmp, (char *)"cella");
    double b = lammps_get_thermo(lmp, (char *)"cellb");
    double c = lammps_get_thermo(lmp, (char *)"cellc");
    double alpha = lammps_get_thermo(lmp, (char *)"cellalpha");
    double beta = lammps_get_thermo(lmp, (char *)"cellbeta");
    double gamma = lammps_get_thermo(lmp, (char *)"cellgamma");
    latticeStruct ls = latticeStruct(a, b, c, alpha, beta, gamma);
    crystal->setLattice(ls);

    /* oneshot */
    lammps_command(lmp, "run 0");
    double pe = lammps_get_thermo(lmp, "pe");

    /* delete LAMMPS instance */
    lammps_close(lmp);

    delete []position;
    delete []boxlo;
    delete []boxhi;

    return pe;
}
