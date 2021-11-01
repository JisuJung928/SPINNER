#ifndef __CONFIG_H__
#define __CONFIG_H__
#include "input.h"

using namespace std;
void RandomGeneration(Input *, Crystal *, int, int);
void LatticeMutation(Input *, Crystal *, int, int);
void Permutation(Input *, Crystal *, int, int);
void Crossover(Input *, Crystal *, int, int);
void SortCrystal(Crystal *, vector<unsigned int>);
// TODO: WriteConfig
#endif
