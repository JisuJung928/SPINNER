#include <iostream>
#include <fstream>
#include <string>
#include "input.h"


#define MAXLINE 256
using namespace std;
void Input::read(string filename)
{
    char line[MAXLINE];

    ifstream fp;
    fp.open(filename);
    if (fp.is_open()) {
        while (!fp.eof()) {
            fp.getline(line, MAXLINE); 
            cout << line << endl;
        }
    } else {
        cout << "no" << endl;
    }
    fp.close();
}
