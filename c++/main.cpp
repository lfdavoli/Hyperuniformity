#include "functions.cpp"
using namespace std;

/*
Parameters
*/
int L = 200;            // lattice size
int delta = 2;          // displacement
int nsamples = 10000;   // samples number
bool optimize = false;   // enable optimized computation if "true"

int main()
{
    srand(time(NULL));
    auto start = chrono::high_resolution_clock::now();
    string output_file = to_string(L)+".csv";
    double* lattice = CreateLattice(L);
    
    if( delta != -1)
    {
        output_file = to_string(L)+"_"+to_string(delta)+".csv";
        AddDisplacement(lattice,delta,L);
    }

    if( optimize == false){
        GetVariance(lattice,L,output_file,nsamples);
    }
    else{
        GetVarianceOpt(lattice,L,"opt"+output_file,nsamples);
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    cout<<"Total execution time(ms): "<<duration.count()<<endl;
    return 0;
}