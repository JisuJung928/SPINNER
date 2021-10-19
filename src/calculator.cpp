#include <cmath>
#include "calculator.h"
#include "crystal.h"


#define _USE_MATH_DEFINES
void *LammpsInit(Input input, Crystal crystal, int lmpargc, char **lmpargv)
{
    /* split comm
    int lammps;
    if (me < nprocs_lammps) lammps = 1;
    else lammps = MPI_UNDEFINED;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
    */

    void *lmp = nullptr;
    lmp = lammps_open(lmpargc, lmpargv, MPI_COMM_WORLD, nullptr);

    char cmd[256];
    /* basic */
    const char *cmds[] = {"units metal",
                          "box tilt large",
                          "neigh_modify every 1 delay 0 check yes"};
    lammps_commands_list(lmp, sizeof(cmds) / sizeof(const char *), cmds);
    /* box */
    latticeStruct lattice = crystal.getLattice();
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
    sprintf(cmd, "create_box %d cell", input.GetNelement());
    lammps_command(lmp, cmd);
    /* atoms */
    int *id = new int[crystal.numAtoms()];
    int *type = new int[crystal.numAtoms()];
    double *position = new double[3 * crystal.numAtoms()];

    int z_number = input.GetZNumber();
    int nelement = input.GetNelement();
    vector<int> composition = input.GetComposition();
    vector<atomStruct> atoms = crystal.getAtoms();
    
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
    lammps_create_atoms(lmp, crystal.numAtoms(), id, type, position,
                        nullptr, nullptr, 0);
    /* mass */
    vector<double> mass = input.GetMass();
    for (i = 0; i < nelement; ++i) {
        sprintf(cmd, "mass %d %f", i + 1, mass[i]);
        lammps_command(lmp, cmd);
    }

    delete []id;
    delete []type;
    delete []position;

    return lmp;
}
