#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
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
    string line;
    string word;

    ifstream in(filename);
    if (in.is_open()) {
        while (in) {
            getline(in, line); 
            if (line.compare("") == 0) {
                continue;
            }
            /* replace = with whitespace */
            replace(line.begin(), line.end(), '=', ' ');
            istringstream iss;
            iss.str(line);
            /* skip whitespaces */
            getline(iss >> ws, word, ' ');
            if (word.substr(0, 1).compare("#") == 0) {
                continue;
            } else if (word.compare("ELEMENT") == 0) {
                vector<string> element;
                while (getline(iss >> ws, word, ' ')) {
                    element.push_back(word);
                }
                input->SetElement(element);
                input->SetNelement(input->GetElement().size());
            } else if (word.compare("COMPOSITION") == 0) {
                vector<int> composition;
                while (getline(iss >> ws, word, ' ')) {
                    composition.push_back(stoi(word));
                }
                input->SetComposition(composition);
            } else if (word.compare("Z_NUMBER") == 0) {
                getline(iss >> ws, word);
                input->SetZNumber(stoi(word));
            } else if (word.compare("VOLUME") == 0) {
                getline(iss >> ws, word);
                input->SetVolume(stof(word));
            } else if (word.compare("PAIR_STYLE") == 0) {
                getline(iss >> ws, word);
                input->SetPairStyle(word);
            } else if (word.compare("PAIR_COEFF") == 0) {
                getline(iss >> ws, word);
                input->SetPairCoeff(word);
            } else if (word.compare("MAX_FORCE") == 0) {
                getline(iss >> ws, word);
                input->SetMaxForce(stof(word));
            } else if (word.compare("RELAX_ITER") == 0) {
                getline(iss >> ws, word);
                input->SetRelaxIteration(stoi(word));
            } else if (word.compare("GENERATION") == 0) {
                getline(iss >> ws, word);
                input->SetGeneration(stoi(word));
            } else if (word.compare("POPULATION") == 0) {
                getline(iss >> ws, word);
                input->SetPopulation(stoi(word));
            } else if (word.compare("MAX_POPULATION") == 0) {
                getline(iss >> ws, word);
                input->SetMaxPopulation(stoi(word));
            } else if (word.compare("INIT_WINDOW") == 0) {
                getline(iss >> ws, word);
                input->SetInitWindow(stof(word));
            } else if (word.compare("GENE_WINDOW") == 0) {
                getline(iss >> ws, word);
                input->SetGeneWindow(stof(word));
            } else if (word.compare("BEST_WINDOW") == 0) {
                getline(iss >> ws, word);
                input->SetBestWindow(stof(word));
            } else if (word.compare("RANDOM_GEN") == 0) {
                getline(iss >> ws, word);
                input->SetRandomGen(stof(word));
            } else if (word.compare("CROSSOVER") == 0) {
                getline(iss >> ws, word);
                input->SetCrossover(stof(word));
            } else if (word.compare("PERMUTATION") == 0) {
                getline(iss >> ws, word);
                input->SetPermutation(stof(word));
            } else if (word.compare("LATTICE_MUT") == 0) {
                getline(iss >> ws, word);
                input->SetLatticeMut(stof(word));
            } else if (word.compare("CONSTRAINT") == 0) {
                vector<double> constraint;
                while (getline(iss >> ws, word, ' ')) {
                    constraint.push_back(stof(word));
                }
            } else if (word.compare("NPAR") == 0) {
                getline(iss >> ws, word);
                input->SetNpar(stoi(word));
            } else if (word.compare("RANDOM_SEED") == 0) {
                getline(iss >> ws, word);
                input->SetRandomSeed(stoi(word));
            } else {
                cout << "Check the input tag!" << endl;
            }
        }
    } else {
        cout << "Check the filename!" << endl;
    }

    /* mass */
    vector<double> mass;
    vector<string> element = input->GetElement();
    for (int i = 0; i < input->GetNelement(); ++i) {
        mass.push_back(GetMassFromSymbol(element[i]));
    }
    input->SetMass(mass);

    return input;
}
