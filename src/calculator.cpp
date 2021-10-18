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
    cell[2][1] = (lattice.b * lattice.c * cos(lattice.alpha * M_PI / 180) - cell[1][0] * cell[2][0]) / cell[1][1];
    cell[2][2] = sqrt(lattice.c * lattice.c - cell[2][0] * cell[2][0] - cell[2][1] * cell[2][1]);

    sprintf(cmd, "region cell prism 0 %f 0 %f 0 %f %f %f %f units box",
            cell[0][0], cell[1][1], cell[2][2],
            cell[1][0], cell[2][0], cell[2][1]);
    lammps_command(lmp, cmd);
    sprintf(cmd, "create_box %d cell", input.GetNelement());
    lammps_command(lmp, cmd);

    return lmp;
}
