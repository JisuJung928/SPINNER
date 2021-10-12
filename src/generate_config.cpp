#include "generate_config.h"


#define MAXLINE 256
using namespace std;
Crystal generate(string filename)
{
    RandSpgOptions options = RandSpgOptions::readOptions(filename);

    vector<uint> atoms;
    string comp = options.getComposition();

    ElemInfo::readComposition(comp, atoms);

    e_verbosity = options.getVerbosity();

    latticeStruct mins = options.getLatticeMins();
    latticeStruct maxes = options.getLatticeMaxes();

    randSpgInput input(1, atoms, mins, maxes);

    Crystal c = RandSpg::randSpgCrystal(input);
    /* test */
    string config_name = "./POSCAR";
    string title = "./TEST_CONFIG";
    c.writePOSCAR(config_name, title);

    return c;
}
