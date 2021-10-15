#include <random>
#include <vector>

using namespace std;

// Add [dx,dy] random displacement to each x[i]
void AddDisplacement(vector<vector<double>> &lattice, double delta) {
   int N = sizeof(lattice);
   for(int j = 0; j = N; j++) {
      double dx = rand()/RAND_MAX - delta/2;
      double dy = rand()/RAND_MAX - delta/2;
      vector<double> displacement = vector<double>{dx,dy}; 
      lattice[j] = lattice[j] + displacement;
   }
   return;
}