#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "input.h"
#include "utils.h"

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
