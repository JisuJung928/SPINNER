#include <iostream>
#include <fstream>
#include <cstring>
#include "input.h"

#define MAXLINE 256
using namespace std;
Input ReadInput(string filename)
{
    Input input;
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
                input.SetElement(element);
            } else if (strcmp(ptr, "COMPOSITION") == 0) {
                strtok(nullptr, " \n\t");
                ptr = strtok(nullptr, " \n\t");
                vector<uint> composition;
                while (ptr != nullptr) {
                    composition.push_back(atoi(ptr));
                    ptr = strtok(nullptr, " \n");
                }
                input.SetComposition(composition);
            } else if (strcmp(ptr, "Z_NUMBER") == 0) {
                strtok(nullptr, " \n\t");
                input.SetZnumber(atoi(strtok(nullptr, " \n")));
            } else if (strcmp(ptr, "VOLUME") == 0) {
                strtok(nullptr, " \n\t");
                input.SetVolume(atof(strtok(nullptr, " \n")));
            } else if (strcmp(ptr, "GENERATION") == 0) {
                strtok(nullptr, " \n\t");
                input.SetGeneration(atoi(strtok(nullptr, " \n")));
            } else if (strcmp(ptr, "POPULATION") == 0) {
                strtok(nullptr, " \n\t");
                input.SetPopulation(atoi(strtok(nullptr, " \n")));
            } else if (strcmp(ptr, "POT_PATH") == 0) {
                strtok(nullptr, " \n\t");
                input.SetPotentialPath(string(strtok(nullptr, " \n")));
            } else if (strcmp(ptr, "RANDOM_SEED") == 0) {
                strtok(nullptr, " \n\t");
                input.SetRandomSeed(atoi(strtok(nullptr, " \n")));
            } else {
                cout << "Check the input tag!" << endl;
            }
        }
    } else {
        cout << "Check the filename!" << endl;
    }
    fp.close();

    return input;
}
