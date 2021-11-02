#include <iostream>
#include <fstream>
#include "input.h"

double GetMassFromSymbol(string symbol)
{
    const vector<pair<string, double>> mass =
    {
        {"H", 1.008},
        {"He", 4.003},
        {"Li", 6.941},
        {"Be", 9.012},
        {"B", 10.811},
        {"C", 12.011},
        {"N", 14.007},
        {"O", 15.999},
        {"F", 18.998},
        {"Ne", 20.180},
        {"Na", 22.990},
        {"Mg", 24.305},
        {"Al", 26.982},
        {"Si", 28.086},
        {"P", 30.974},
        {"S", 32.065},
        {"Cl", 35.453},
        {"Ar", 39.948},
        {"K", 39.098},
        {"Ca", 40.078},
        {"Sc", 44.956},
        {"Ti", 47.867},
        {"V", 50.942},
        {"Cr", 51.996},
        {"Mn", 54.938},
        {"Fe", 55.845},
        {"Co", 58.933},
        {"Ni", 58.693},
        {"Cu", 63.546},
        {"Zn", 65.380},
        {"Ga", 69.723},
        {"Ge", 72.640},
        {"As", 74.922},
        {"Se", 78.960},
        {"Br", 79.904},
        {"Kr", 83.798},
        {"Rb", 85.468},
        {"Sr", 87.620},
        {"Y", 88.906},
        {"Zr", 91.224},
        {"Nb", 92.906},
        {"Mo", 95.960},
        {"Tc", 98.000},
        {"Ru", 101.07},
        {"Rh", 102.91},
        {"Pd", 106.42},
        {"Ag", 107.87},
        {"Cd", 112.41},
        {"In", 114.82},
        {"Sn", 118.71},
        {"Sb", 121.76},
        {"Te", 127.60},
        {"I", 126.90},
        {"Xe", 131.29},
        {"Cs", 132.91},
        {"Ba", 137.33},
        {"La", 138.91},
        {"Ce", 140.12},
        {"Pr", 140.91},
        {"Nd", 144.24},
        {"Pm", 145.00},
        {"Sm", 150.36},
        {"Eu", 151.96},
        {"Gd", 157.25},
        {"Tb", 158.93},
        {"Dy", 162.50},
        {"Ho", 164.93},
        {"Er", 167.26},
        {"Tm", 168.93},
        {"Yb", 173.05},
        {"Lu", 174.97},
        {"Hf", 178.49},
        {"Ta", 180.95},
        {"W", 183.84},
        {"Re", 186.21},
        {"Os", 190.23},
        {"Ir", 192.22},
        {"Pt", 195.08},
        {"Au", 196.97},
        {"Hg", 200.59},
        {"Tl", 204.38},
        {"Pb", 207.20},
        {"Bi", 208.98},
        {"Po", 209.00},
        {"At", 210.00},
        {"Rn", 222.00},
        {"Fr", 223.00},
        {"Ra", 226.00},
        {"Ac", 227.00},
        {"Th", 232.04},
        {"Pa", 231.04},
        {"U", 238.03}
    };
    for (auto i : mass) {
        if (i.first == symbol) {
            return i.second;
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
        // TODO: assertion
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
            } else if (strcmp(ptr, "RELAX_ITER") == 0) {
                strtok(nullptr, " \n\t");
                input->SetRelaxIteration(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "GENERATION") == 0) {
                strtok(nullptr, " \n\t");
                input->SetGeneration(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "POPULATION") == 0) {
                strtok(nullptr, " \n\t");
                input->SetPopulation(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "MAX_POPULATION") == 0) {
                strtok(nullptr, " \n\t");
                input->SetMaxPopulation(atoi(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "INIT_WINDOW") == 0) {
                strtok(nullptr, " \n\t");
                input->SetInitWindow(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "GENE_WINDOW") == 0) {
                strtok(nullptr, " \n\t");
                input->SetGeneWindow(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "BEST_WINDOW") == 0) {
                strtok(nullptr, " \n\t");
                input->SetBestWindow(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "RANDOM_GEN") == 0) {
                strtok(nullptr, " \n\t");
                input->SetRandomGen(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "CROSSOVER") == 0) {
                strtok(nullptr, " \n\t");
                input->SetCrossover(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "PERMUTATION") == 0) {
                strtok(nullptr, " \n\t");
                input->SetPermutation(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "LATTICE_MUT") == 0) {
                strtok(nullptr, " \n\t");
                input->SetLatticeMut(atof(strtok(nullptr, "\n")));
            } else if (strcmp(ptr, "CONSTRAINT") == 0) {
                strtok(nullptr, " \n\t");
                ptr = strtok(nullptr, " \n\t");
                vector<double> constraint;
                while (ptr != nullptr) {
                    constraint.push_back(atof(ptr));
                    ptr = strtok(nullptr, " \n");
                }
                input->SetConstraint(constraint);
            } else if (strcmp(ptr, "NPAR") == 0) {
                strtok(nullptr, " \n\t");
                input->SetNpar(atoi(strtok(nullptr, "\n")));
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
