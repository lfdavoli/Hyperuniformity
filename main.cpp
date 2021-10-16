#include <fstream>
#include "functions.cpp"
using namespace std;

int main()
{
    srand(time(NULL));
    
    for (size_t i = 0; i < 2000000; i++)
    {
        rand();
    } 
}