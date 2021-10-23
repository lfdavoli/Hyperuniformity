#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <time.h> 
#include <string>
#include <chrono>
#include <algorithm>

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
        coords.push_back(i/L);
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
double compute_variance_R(double R, vector<vector<double>> lattice)
{
    double sigma_sq;
    int L = sqrt(lattice.size());

    int iterations = 10000;
    int N = 0;
    int N_sq = 0;

    for (size_t i = 0; i < iterations; i++)
    {
        double X_x0 = ((double) rand() / (RAND_MAX)) * L;
        double Y_x0 = ((double) rand() / (RAND_MAX)) * L;

        int points_inside = 0;
        //cout<<"Variance "<<i<<endl;
        for (size_t j = 0; j < L*L; j++)
        {
            //cout<<"Variance dentro "<<j<<endl;
            if(dist(X_x0, Y_x0, j, lattice) < R)
            {
                points_inside ++;
            }
        }
        //cout<<"Variance uscito "<<i<<endl;
        N += points_inside;
        N_sq += points_inside*points_inside;
    }
    sigma_sq = N_sq/iterations - (N/iterations)*(N/iterations);
    
    return sigma_sq;
}

/*
Input-> x0 and Lattice with point positions
Output-> Variance

Given x0 it computes N(R) and N(R)^2 for all radii
*/
vector<int> compute_N_r_x0(double X_x0, double Y_x0, vector<vector<double>> lattice)
{
    double sigma_sq;
    int L = sqrt(lattice.size());
    vector<int> N_x0;
    vector<double> distances;

    double radii[200];
    radii[0] = 0.01;
    radii[199] = L/2;
    double c = exp(log(radii[199]/radii[0])/(200-1));
    for (size_t i = 1; i < 200-1; i++)
    {
        radii[i] = radii[i-1]*c;
    }

    for (size_t j = 0; j < L*L; j++)
    {
        distances.push_back(dist(X_x0, Y_x0, j, lattice));

    }

    sort(distances.begin(), distances.end());

    int new_first = 0;
    int counter = 0;
    int points_inside;

    for (auto &&r : radii)
    {
        points_inside = 0;
        for (int j = new_first; j < distances.size(); j++)
        {
            if(distances[j] > r)
            {
                points_inside = max(j - 1, 0);
                new_first = j;
                break;
            }
        }
        N_x0.push_back(points_inside);
    }
    
    return N_x0;
}


void get_variance_R(int lattice_size)
{
    auto start = chrono::high_resolution_clock::now();

    vector<vector<double>> lattice = create_lattice(lattice_size);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    //cout<<"Create_lattice: "<<duration.count()<<endl;
    double radii[200];
    radii[0] = 0.01;
    radii[199] = lattice_size/2;
    double c = exp(log(radii[199]/radii[0])/(200-1));
    for (size_t i = 1; i < 200-1; i++)
    {
        radii[i] = radii[i-1]*c;
    }

    string file_name = "dati.csv"; 
    ofstream output(file_name);

    for (auto &&i : radii)
    {
        double sigma_sq;
        auto start = chrono::high_resolution_clock::now();
        
        sigma_sq = compute_variance_R(i, lattice)/(i*i);

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::microseconds>(stop-start);
        //cout<<"compute_variance: "<<duration.count()<<endl;

        output<<i<<" "<<sigma_sq<<endl;
    }
}

void get_variance_x0(int lattice_size)
{
    vector<vector<double>> lattice = create_lattice(lattice_size);
    vector<int> N_x0;
    vector<double> N(200, 0);
    vector<double> N_squared(200, 0);
    double sigma;
    int N_ITER = 10000;
    char file_name = (char)lattice_size;
    ofstream output(file_name + ".csv");
    

    for (int i = 0; i < N_ITER; i++)
    {
        double X_x0 = ((double) rand() / (RAND_MAX)) * lattice_size;
        double Y_x0 = ((double) rand() / (RAND_MAX)) * lattice_size;
        N_x0 = compute_N_r_x0(X_x0, Y_x0, lattice);
        //cout<<N_x0.size()<<endl;
        for (size_t i = 0; i < N_x0.size(); i++)
        {
            N[i] = (N[i] + N_x0[i]);
            N_squared[i] = (N_squared[i] + N_x0[i]*N_x0[i]);
        }
    }

    double radii[200];
    radii[0] = 0.01;
    radii[199] = lattice_size/2;
    double c = exp(log(radii[199]/radii[0])/(200-1));
    for (size_t i = 1; i < 200-1; i++)
    {
        radii[i] = radii[i-1]*c;
    }
    //cout<<N.size()<<endl;
    for(int i = 0; i<N.size(); i++)
    {
        //cout<<radii[i]*radii[i]<<endl;
        sigma = ((N_squared[i]/N_ITER) - (N[i]/N_ITER)*(N[i]/N_ITER))/(radii[i]*radii[i]);
        output<<radii[i]<<" "<<sigma<<endl;
    }
    
}


/*
Add [dx,dy] random displacement to each lattice site
*/
void AddDisplacement(vector<vector<double>> &lattice, double delta) {
   int N = lattice.size();
   for(int j = 0; j != N-1; j++) {
      double dx = ( (double)rand() )/RAND_MAX  - delta/2;
      double dy = ( (double)rand() )/RAND_MAX  - delta/2;
      lattice[j] = vector<double>{dx+lattice[j][0],dy+lattice[j][1]};
   }
   return;
}