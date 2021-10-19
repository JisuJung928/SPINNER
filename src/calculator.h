#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include "randSpg.h"

#include "input.h"

void *LammpsInit(Input *, Crystal *, int, char **);
void Relax(Input *, Crystal *);
#endif
