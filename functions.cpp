#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <time.h> 

using namespace std;

/*
Input-> L=size of lattice
*/
vector<vector<double>> create_lattice(int L)
{
    vector<vector<double>> points;
    for (size_t i = 0; i < L*L; i++)
    {
        vector<double> coords;
        coords.push_back(i%L);
        coords.push_back(floor(i/L));
        points.push_back(coords);
    }
    return points;
}

/*
Compute distance between two points on lattice
*/
double dist(double X_x0, double Y_x0, int iter, vector<vector<double>> lattice)
{
    int L = sqrt(lattice.size());

    double X_x = lattice[iter][0];
    double Y_x = lattice[iter][1];
    double distx, disty;
    
    if(abs(X_x - X_x0) < L/2)
    {
        distx = abs(X_x - X_x0);
    }
    else
    {
        distx = L - abs(X_x - X_x0);
    }

    if(abs(Y_x - Y_x0) < L/2)
    {
        disty = abs(Y_x - Y_x0);
    }
    else
    {
        disty = L - abs(Y_x - Y_x0);
    }

    return sqrt(distx*distx + disty*disty);
}


/*
Input-> Radius and Lattice with point positions
Output-> Variance

compute number of points inside given circle centered in x0
sample 10000 x0
*/
double variance(double R, vector<vector<double>> lattice)
{
    double sigma_sq;
    int L = lattice.size();

    int iterations = 10000;
    int N = 0;
    int N_sq = 0;

    for (size_t i = 0; i < iterations; i++)
    {
        double X_x0 = ((double) rand() / (RAND_MAX)) * sqrt(L);
        double Y_x0 = ((double) rand() / (RAND_MAX)) * sqrt(L);

        int points_inside = 0;
        for (size_t j = 0; i < L; i++)
        {
            if(dist(X_x0, Y_x0, j, lattice) < R)
            {
                points_inside ++;
            }
        }
        N += points_inside;
        N_sq += points_inside*points_inside;
    }
    sigma_sq = N_sq/iterations - (N/iterations)*(N/iterations);
    
    return sigma_sq;
}