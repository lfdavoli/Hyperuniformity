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
    auto start = chrono::high_resolution_clock::now();
    get_variance_x0(L);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop-start);
    cout<<"Total duration: "<<duration.count()<<endl;
    return 0;
}