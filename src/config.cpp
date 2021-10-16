#include "config.h"


#define MAXLINE 256
using namespace std;
Crystal generate(Input input)
{
    string comp = "";
    vector<string> element = input.get_element();
    vector<uint> composition = input.get_composition();
    for (unsigned int i = 0; i < element.size(); ++i) {
        comp += element[i];
        comp += to_string(composition[i]);
    }

    vector<uint> atoms;
    ElemInfo::readComposition(comp, atoms);

    latticeStruct lmin = latticeStruct(3.0, 3.0, 3.0, 60.0, 60.0, 60.0);
    latticeStruct lmax = latticeStruct(10.0, 10.0, 10.0, 120.0, 120.0, 120.0);


    randSpgInput spg_input(1, atoms, lmin, lmax);
    /* caution: sort automatically */
    Crystal c = RandSpg::randSpgCrystal(spg_input);

    return c;
}
