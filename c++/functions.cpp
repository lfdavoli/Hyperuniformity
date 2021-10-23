#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <time.h> 
#include <string>
#include <chrono>
#include <algorithm>
#include <omp.h>

using namespace std;

/*
Input-> L=size of lattice
*/
double* create_lattice(int L)
{
    //vector<vector<double>> points;
    //for (size_t i = 0; i < L*L; i++)
    //{
    //    vector<double> coords;
    //    coords.push_back(i%L);
    //    coords.push_back(i/L);
    //    points.push_back(coords);
    //}

    double *points = new double[2*L*L];
    for (size_t i = 0; i < L*L; i++)
    {
        points[i*2] = i%L;
        points[i*2+1] = i/L;
        //cout<<i%L<<" "<<i/L<<endl;
    }
    
    return points;
}


/*
Compute distance between two points on lattice
*/
double dist(double X_x0, double Y_x0, int iter, double* lattice, int L)
{
    //int L = sqrt(sizeof(lattice)/sizeof(lattice[0]));

    double X_x = lattice[iter*2];
    double Y_x = lattice[iter*2 + 1];
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
double compute_variance_R(double R, double* lattice, int L)
{
    double sigma_sq;
    //int L = sqrt(sizeof(lattice)/sizeof(lattice[0]));

    int iterations = 10000;
    int N = 0;
    int N_sq = 0;

    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < iterations; i++)
    {
        double X_x0 = ((double) rand() / (RAND_MAX)) * L;
        double Y_x0 = ((double) rand() / (RAND_MAX)) * L;

        int points_inside = 0;
        //cout<<"Variance "<<i<<endl;
        for (size_t j = 0; j < L*L; j++)
        {
            //cout<<"Variance dentro "<<j<<endl;
            if(dist(X_x0, Y_x0, j, lattice, L) < R)
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
int* compute_N_r_x0(double X_x0, double Y_x0, double* lattice, double radii[200], int L)
{
    double sigma_sq;
    //int L = sqrt((sizeof(lattice)/sizeof(double))/2);
    //cout<<L<<endl;
    int *N_x0 = new int[200];
    double distances[L*L];

    for (size_t j = 0; j < L*L; j++)
    {
        distances[j] = (dist(X_x0, Y_x0, j, lattice, L));
        
    }
    //int size_dist = sizeof(distances)/sizeof(distances[0]);

    sort(distances, distances + L*L);

    int new_first = 0;
    int points_inside;

    for (int i = 0; i<200; i++)
    {
        points_inside = 0;
        for (int j = new_first; j < L*L; j++)
        {
            if(distances[j] > radii[i])
            {
                points_inside = max(j - 1, 0);
                new_first = j;
                break;
            }
        }
        //cout<<points_inside<<endl;
        N_x0[i]= points_inside;
    }
    
    return N_x0;
}


void get_variance_R(int lattice_size)
{
    auto start = chrono::high_resolution_clock::now();

    double* lattice = create_lattice(lattice_size);

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

    ofstream output(to_string(lattice_size) + ".csv");

    for (auto &&i : radii)
    {
        double sigma_sq;
        auto start = chrono::high_resolution_clock::now();
        
        sigma_sq = compute_variance_R(i, lattice, lattice_size)/(i*i);

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::microseconds>(stop-start);
        //cout<<"compute_variance: "<<duration.count()<<endl;

        output<<i<<" "<<sigma_sq<<endl;
    }
}

void get_variance_x0(int lattice_size)
{
    double* lattice = create_lattice(lattice_size);
    int *N_x0;
    double N[200] = {0};
    double N_squared[200] = {0};
    double sigma;
    int N_ITER = 10000;
    ofstream output(to_string(lattice_size) + ".csv");

    double radii[200];
    radii[0] = 0.01;
    radii[199] = lattice_size/2;
    double c = exp(log(radii[199]/radii[0])/(200-1));
    for (size_t i = 1; i < 200-1; i++)
    {
        radii[i] = radii[i-1]*c;
    }

    #pragma omp parallel for default(shared)

    for (int i = 0; i < N_ITER; i++)
    {
        double X_x0 = ((double) rand() / (RAND_MAX)) * lattice_size;
        double Y_x0 = ((double) rand() / (RAND_MAX)) * lattice_size;
        N_x0 = compute_N_r_x0(X_x0, Y_x0, lattice, radii, lattice_size);
        
        
        for (size_t j = 0; j < 200; j++)
        {
            N[j] = (N[j] + N_x0[j]);
            N_squared[j] = (N_squared[j] + N_x0[j]*N_x0[j]);   
        }
        
    }

    
    //cout<<N.size()<<endl;
    for(int i = 0; i<200; i++)
    {
        //cout<<radii[i]*radii[i]<<endl;
        sigma = ((N_squared[i]/N_ITER) - (N[i]/N_ITER)*(N[i]/N_ITER))/(radii[i]*radii[i]);
        output<<radii[i]<<" "<<sigma<<endl;
    }
    output.close();
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