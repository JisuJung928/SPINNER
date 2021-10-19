#include <iostream>
#include <mpi.h>

#include "calculator.h"
#include "config.h"
#include "input.h"

using namespace std;
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Input *input = ReadInput("./INPUT");

    vector<Crystal> crystal_vector = GenerateCrystal(input);
    /* test */
    string config_name = "./POSCAR";
    string title = "TEST_CONFIG";
    for (unsigned int i = 0; i < crystal_vector.size(); ++i) {
        /* writePOSCAR resort again */
        crystal_vector[i].writePOSCAR(config_name + to_string(i), title);
    }

    Relax(input, &crystal_vector[0]);

    //TODO: MPI communication of object (or struct)

    delete input;

    return 0;
}
