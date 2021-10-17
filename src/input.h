#ifndef __INPUT_H__
#define __INPUT_H__
#include <string>
#include <vector>

using namespace std;
class Input
{
    private:
        vector<string> element;
        vector<uint> composition;
        int z_number;
        double volume;

        int generation;
        int population;

        string pot_path;

        int random_seed;

    public:
        /* setter */
        void set_element(const vector<string> &v)
        {
            element = v;
        }
        void set_composition(const vector<uint> &v)
        {
            composition = v;
        }
        void set_z_number(int i)
        {
            z_number = i;
        }
        void set_volume(double d)
        {
            volume = d;
        }
        void set_generation(int i)
        {
            generation = i;
        }
        void set_population(int i)
        {
            population = i;
        }
        void set_pot_path(const string &s)
        {
            pot_path = s;
        }
        void set_random_seed(int i)
        {
            random_seed = i;
        }
        /* getter */
        vector<string> get_element() const
        {
            return element;
        }
        vector<uint> get_composition() const
        {
            return composition;
        }
        int get_z_number() const
        {
            return z_number;
        }
        double get_volume() const
        {
            return volume;
        }
        int get_generation() const
        {
            return generation;
        }
        int get_population() const
        {
            return population;
        }
        string get_pot_path() const
        {
            return pot_path;
        }
        int get_random_seed() const
        {
            return random_seed;
        }
};
Input read_input(string);

#endif
