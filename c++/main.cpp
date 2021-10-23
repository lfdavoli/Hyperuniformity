#include "functions.cpp"
using namespace std;

/*
Parameters
*/
int L = 10;
double delta = 1;

int main()
{
    srand(time(NULL));
    //vector<vector<double>> lattice = create_lattice(L); 
    get_variance_R(L);
    return 0;
}