#include "functions.cpp"
using namespace std;

/*
Parameters
*/
int L = 200;
double delta = 1;

int main()
{
    srand(time(NULL));
    auto start = chrono::high_resolution_clock::now();
    double* lattice = create_lattice(lattice_size);
    
    get_variance_R(L);
    //get_variance_x0(L);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    cout<<"Total duration: "<<duration.count()<<endl;
    return 0;
}