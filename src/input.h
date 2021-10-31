#ifndef __INPUT_H__
#define __INPUT_H__
#include <cstring>
#include <string>
#include <vector>

using namespace std;
class Input
{
    private:
        /* structure */
        vector<string> element;
        vector<int> composition;
        int z_number;
        int nelement;
        double volume;
        vector<double> mass;

        /* potential */
        char *pair_style;
        char *pair_coeff;

        /* calculation */
        double max_force;

        /* evolution */
        int generation;
        int population;

        /* operator */
        double random_gen;
        double crossover;
        double permutation;
        double lattice_mut;

        /* constraint */
        vector<double> constraint;

        /* parallelism */
        int npar;

        /* random */
        int random_seed;

    public:
        Input()
        {
            pair_style = new char[128];
            pair_coeff = new char[128];
            random_seed = -1;
        }
        ~Input(){
            delete []pair_style;
            delete []pair_coeff;
        }

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
        void SetPairStyle(char *c)
        {
            strcpy(pair_style, c);
        }
        void SetPairCoeff(char *c)
        {
            strcpy(pair_coeff, c);
        }
        void SetMaxForce(double d)
        {
            max_force = d;
        }
        void SetGeneration(int i)
        {
            generation = i;
        }
        void SetPopulation(int i)
        {
            population = i;
        }
        void SetRandomGen(double d)
        {
            random_gen = d;
        }
        void SetCrossover(double d)
        {
            crossover = d;
        }
        void SetPermutation(double d)
        {
            permutation = d;
        }
        void SetLatticeMut(double d)
        {
            lattice_mut = d;
        }
        void SetConstraint(const vector<double> &v)
        {
            constraint = v;
        }
        void SetNpar(int i)
        {
            npar = i;
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
        double GetMaxForce() const
        {
            return max_force;
        }
        char *GetPairStyle() const
        {
            return pair_style;
        }
        char *GetPairCoeff() const
        {
            return pair_coeff;
        }
        int GetGeneration() const
        {
            return generation;
        }
        int GetPopulation() const
        {
            return population;
        }
        double GetRandomGen() const
        {
            return random_gen;
        }
        double GetCrossover() const
        {
            return crossover;
        }
        double GetPermutation() const
        {
            return permutation;
        }
        double GetLatticeMut() const
        {
            return lattice_mut;
        }
        vector<double> GetConstraint() const
        {
            return constraint;
        }
        int GetNpar() const
        {
            return npar;
        }
        int GetRandomSeed() const
        {
            return random_seed;
        }
};
double GetMassFromSymbol(string);
Input *ReadInput(string);

#endif
