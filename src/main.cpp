#include <iostream>
#include <string>
#include <mpi.h>

#define LAMMPS_LIB_MPI
#include "library.h"

#include "elemInfo.h"
#include "randSpg.h"

using namespace std;
int main(int argc, char* argv[])
{
    RandSpgOptions options = RandSpgOptions::readOptions(argv[1]);

    vector<uint> atoms;
    string comp = options.getComposition();

    ElemInfo::readComposition(comp, atoms);

    e_verbosity = options.getVerbosity();

    latticeStruct mins = options.getLatticeMins();
    latticeStruct maxes = options.getLatticeMaxes();

    randSpgInput input(1, atoms, mins, maxes);

    Crystal c = RandSpg::randSpgCrystal(input);
    /* test */
    string filename = "./POSCAR";
    string title = "./TEST_CONFIG";
    c.writePOSCAR(filename, title);

    return 0;
}
