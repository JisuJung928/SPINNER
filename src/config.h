#ifndef __CONFIG_H__
#define __CONFIG_H__
#include "elemInfo.h"
#include "randSpg.h"

#include "input.h"

using namespace std;
vector<Crystal> RandomGeneration(Input *, int);
vector<Crystal> LatticeMutation(Input *, int);
vector<Crystal> Permutation(Input *, int);
vector<Crystal> Crossover(Input *, int);
void SortCrystal(Crystal *, vector<unsigned int>);
// TODO: WriteConfig
#endif
