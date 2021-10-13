#include <iostream>
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

    vector<string> element = input.get_element();
    for (unsigned int i = 0; i < element.size(); ++i) {
        cout << element[i] << endl;
    }

    generate("randSpg.in");
    return 0;
}
