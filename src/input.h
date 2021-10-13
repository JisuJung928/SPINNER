#ifndef __INPUT_H__
#define __INPUT_H__
#include <string>
#include <vector>

using namespace std;
class Input
{
    private:
        vector<string> element;
        vector<int> composition;
        int z_number;

        int generation;
        int population;

        string pot_path;

    public:
        void init()
        {
            generation = 1;
            population = 1;
            pot_path = "./potential_saved";
        }
        void read(string);
        /* for debug */
        vector<string> get_element()
        {
            return element;
        }
        vector<int> get_composition()
        {
            return composition;
        }
        int get_z_number()
        {
            return z_number;
        }
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
        }
};


#endif
