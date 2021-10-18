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

        string potential_path;

        int random_seed;

    public:
        /* setter */
        void SetElement(const vector<string> &v)
        {
            element = v;
        }
        void SetComposition(const vector<uint> &v)
        {
            composition = v;
        }
        void SetZnumber(int i)
        {
            z_number = i;
        }
        void SetVolume(double d)
        {
            volume = d;
        }
        void SetGeneration(int i)
        {
            generation = i;
        }
        void SetPopulation(int i)
        {
            population = i;
        }
        void SetPotentialPath(const string &s)
        {
            potential_path = s;
        }
        void SetRandomSeed(int i)
        {
            random_seed = i;
        }
        /* getter */
        vector<string> GetElement() const
        {
            return element;
        }
        vector<uint> GetComposition() const
        {
            return composition;
        }
        int GetZnumber() const
        {
            return z_number;
        }
        double GetVolume() const
        {
            return volume;
        }
        int GetGeneration() const
        {
            return generation;
        }
        int GetPopulation() const
        {
            return population;
        }
        string GetPotentialPath() const
        {
            return potential_path;
        }
        int GetRandomSeed() const
        {
            return random_seed;
        }
};
Input ReadInput(string);

#endif
