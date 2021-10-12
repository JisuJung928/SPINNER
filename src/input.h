#ifndef __INPUT_H__
#define __INPUT_H__
#include <string>

using namespace std;
class Input
{
    private:
        int generation;
        int population;
        string pot_path;
    public:
        void init()
        {
            generation = 1;
            population = 1;
        }
        void read(string);
        /* for debug
        int get_generation()
        {
            return generation;
        }
        int get_population()
        {
            return population;
        }
        string get_pot_path()
        {
            return pot_path;
        } */
};


#endif
