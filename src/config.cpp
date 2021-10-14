#include "config.h"


#define MAXLINE 256
using namespace std;
void generate(Input input)
{
    /* original
    RandSpgOptions options = RandSpgOptions::readOptions(filename);

    vector<uint> atoms;
    string comp = options.getComposition();

    ElemInfo::readComposition(comp, atoms);

    latticeStruct mins = options.getLatticeMins();
    latticeStruct maxes = options.getLatticeMaxes();

    randSpgInput input(1, atoms, mins, maxes);

    Crystal c = RandSpg::randSpgCrystal(input);
    */

    string comp = "";
    vector<string> element = input.get_element();
    vector<uint> composition = input.get_composition();
    for (unsigned int i = 0; i < element.size(); ++i) {
        comp += element[i];
        comp += to_string(composition[i]);
    }

    vector<uint> atoms;
    ElemInfo::readComposition(comp, atoms);

    /* test 
    string config_name = "./POSCAR";
    string title = "./TEST_CONFIG";
    c.writePOSCAR(config_name, title);

    return c;
    */
}
