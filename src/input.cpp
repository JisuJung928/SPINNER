#include <iostream>
#include <fstream>
#include "input.h"


double GetMassFromSymbol(string symbol)
{
    const vector<string> symbols =
    {
        " ",
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U"
    };
    const vector<double> mass =
    {
        0.000,
        1.008,
        4.003,
        6.941,
        9.012,
        10.811,
        12.011,
        14.007,
        15.999,
        18.998,
        20.180,
        22.990,
        24.305,
        26.982,
        28.086,
        30.974,
        32.065,
        35.453,
        39.948,
        39.098,
        40.078,
        44.956,
        47.867,
        50.942,
        51.996,
        54.938,
        55.845,
        58.933,
        58.693,
        63.546,
        65.380,
        69.723,
        72.640,
        74.922,
        78.960,
        79.904,
        83.798,
        85.468,
        87.620,
        88.906,
        91.224,
        92.906,
        95.960,
        98.000,
        101.07,
        102.91,
        106.42,
        107.87,
        112.41,
        114.82,
        118.71,
        121.76,
        127.60,
        126.90,
        131.29,
        132.91,
        137.33,
        138.91,
        140.12,
        140.91,
        144.24,
        145.00,
        150.36,
        151.96,
        157.25,
        158.93,
        162.50,
        164.93,
        167.26,
        168.93,
        173.05,
        174.97,
        178.49,
        180.95,
        183.84,
        186.21,
        190.23,
        192.22,
        195.08,
        196.97,
        200.59,
        204.38,
        207.20,
        208.98,
        209.00,
        210.00,
        222.00,
        223.00,
        226.00,
        227.00,
        232.04,
        231.04,
        238.03
    };
    for (uint i = 0; i < symbols.size(); ++i) {
        if (symbols[i] == symbol) {
            return mass[i];
        }
    }
    return -1;
}


#define MAXLINE 256
using namespace std;
Input *ReadInput(string filename)
{
    Input *input = new Input();
    char *ptr, line[MAXLINE];

    ifstream fp;
    fp.open(filename);
    if (fp.is_open()) {
        while (!fp.eof()) {
            fp.getline(line, MAXLINE); 
            if (strcmp(line, "") == 0) {
                continue;
            } else {
                ptr = strtok(line, " \n\t");
            }
            if (strncmp(ptr, "#", 1) == 0) {
                continue;
            } else if (strcmp(ptr, "ELEMENT") == 0) {
                strtok(nullptr, " \n\t");
                ptr = strtok(nullptr, " \n\t");
                vector<string> element;
                while (ptr != nullptr) {
                    element.push_back(string(ptr));
                    ptr = strtok(nullptr, " \n");
                }
                input->SetElement(element);
                input->SetNelement(input->GetElement().size());
            } else if (strcmp(ptr, "COMPOSITION") == 0) {
                strtok(nullptr, " \n\t");
                ptr = strtok(nullptr, " \n\t");
                vector<int> composition;
                while (ptr != nullptr) {
                    composition.push_back(atoi(ptr));
                    ptr = strtok(nullptr, " \n");
                }
                input->SetComposition(composition);
            } else if (strcmp(ptr, "Z_NUMBER") == 0) {
                strtok(nullptr, " \n\t");
                input->SetZNumber(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "VOLUME") == 0) {
                strtok(nullptr, " \n\t");
                input->SetVolume(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "PAIR_STYLE") == 0) {
                strtok(nullptr, " \n\t");
                input->SetPairStyle(strtok(nullptr, "\n"));
            } else if (strcmp(ptr, "PAIR_COEFF") == 0) {
                strtok(nullptr, " \n\t");
                input->SetPairCoeff(strtok(nullptr, "\n"));
            } else if (strcmp(ptr, "MAX_FORCE") == 0) {
                strtok(nullptr, " \n\t");
                input->SetMaxForce(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "GENERATION") == 0) {
                strtok(nullptr, " \n\t");
                input->SetGeneration(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "POPULATION") == 0) {
                strtok(nullptr, " \n\t");
                input->SetPopulation(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "RANDOM_SEED") == 0) {
                strtok(nullptr, " \n\t");
                input->SetRandomSeed(atoi(strtok(nullptr, "\n")));
            } else {
                cout << "Check the input tag!" << endl;
            }
        }
    } else {
        cout << "Check the filename!" << endl;
    }
    fp.close();

    /* mass */
    vector<double> mass;
    vector<string> element = input->GetElement();
    for (int i = 0; i < input->GetNelement(); ++i) {
        mass.push_back(GetMassFromSymbol(element[i]));
    }
    input->SetMass(mass);

    return input;
}
