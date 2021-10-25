#include "functions.cpp"
using namespace std;

/*
Parameters
*/
int L = 200;
double delta = 2;

int main()
{
    srand(time(NULL));
    auto start = chrono::high_resolution_clock::now();
    string output_file = to_string(L)+".csv";
    double* lattice = create_lattice(L);
    if( delta != -1)
    {
        output_file = to_string(L)+"_"+to_string(delta)+".csv";
        AddDisplacement(lattice,delta,L*L);
    }
    get_variance_R(lattice,L,output_file);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    cout<<"Total execution time(ms): "<<duration.count()<<endl;
    return 0;
}