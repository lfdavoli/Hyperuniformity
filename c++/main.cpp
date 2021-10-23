#include "functions.cpp"
using namespace std;

/*
Parameters
*/
int L = 50;
double delta = 2;

int main()
{
    srand(time(NULL));
    auto start = chrono::high_resolution_clock::now();
    double* lattice = create_lattice(L);
    AddDisplacement(lattice,delta,L*L);
    get_variance_R(lattice,L);
    //get_variance_x0(L);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    cout<<"Total duration: "<<duration.count()<<endl;
    return 0;
}