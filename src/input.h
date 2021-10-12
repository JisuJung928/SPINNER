#ifndef __INPUT_H__
#define __INPUT_H__
#include <string>

using namespace std;
class Input
{
    private:
        int generation;
        int population;
    public:
        void init()
        {
            generation = 10;
            population = 10;
        }
        void read(string);
};


#endif
