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
        int nelement;
        double volume;
        vector<double> mass;

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
        void SetComposition(const vector<int> &v)
        {
            composition = v;
        }
        void SetZNumber(int i)
        {
            z_number = i;
        }
        void SetNelement(int i)
        {
            nelement = i;
        }
        void SetVolume(double d)
        {
            volume = d;
        }
        void SetMass(const vector<double> &v)
        {
            mass = v;
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
        vector<int> GetComposition() const
        {
            return composition;
        }
        int GetZNumber() const
        {
            return z_number;
        }
        int GetNelement() const
        {
            return nelement;
        }
        double GetVolume() const
        {
            return volume;
        }
        vector<double> GetMass() const
        {
            return mass;
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
double GetMassFromSymbol(string);
Input ReadInput(string);

#endif
