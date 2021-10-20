#ifndef __CONFIG_H__
#define __CONFIG_H__
#include "elemInfo.h"
#include "randSpg.h"

#include "input.h"

using namespace std;
vector<Crystal> GenerateCrystal(Input *);
void SortCrystal(Crystal *, vector<unsigned int>);
// TODO: WriteConfig
#endif
