#include <iostream>
#include <string>
#include <mpi.h>

#define LAMMPS_LIB_MPI
#include "library.h"

#include "generate_config.h"
#include "input.h"

using namespace std;
int main(int argc, char* argv[])
{
    Input input;
    input.init();
    input.read("./INPUT");

    generate("randSpg.in");
    return 0;
}
