#include <iostream>
#include <mpi.h>

#define LAMMPS_LIB_MPI
#include "library.h"

#include "config.h"
#include "input.h"

using namespace std;
int main(int argc, char* argv[])
{
    Input input;
    read_input(input, "./INPUT");
    generate(input);
    return 0;
}
