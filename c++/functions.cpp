#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <time.h> 
#include <string>
#include <chrono>
#include <algorithm>
#include <omp.h>

using namespace std;

//--------------------------------- Basic functions ----------------------------
/*
Create a square lattice of particles of size L.
Input: L = size of lattice
Output: points = 2L-long array containing particles positions.
*/
double* CreateLattice(int L){
    double *points = new double[2*L*L];
    for (size_t i = 0; i < L*L; i++){
        points[i*2] = i%L;
        points[i*2+1] = i/L;
    }
    return points;
}


/*
Compute distance between two points on lattice.
Input: (X_x0,Y_x0) = position of center of window of observation, index of the observed particle, lattice, lattice size.
Output: distance between the observation point and the particle.
*/
double GetDist_x0_i(double X_x0, double Y_x0, int index, double* lattice, int L){
    double X_x = lattice[index*2];
    double Y_x = lattice[index*2 + 1];
    double distx, disty;
    
    if(abs(X_x - X_x0) < L/2){
        distx = abs(X_x - X_x0);
    }
    else{
        distx = L - abs(X_x - X_x0);
    }
    if(abs(Y_x - Y_x0) < L/2){
        disty = abs(Y_x - Y_x0);
    }
    else{
        disty = L - abs(Y_x - Y_x0);
    }
    return sqrt(distx*distx + disty*disty);
}

/*
Add [dx,dy] uniform random displacement between -delta/2 and delta/2 to each lattice particle.
Input: lattice, displacement, lattice size.
Output: displaced-particles lattice
*/
void AddDisplacement(double* &lattice, double delta, int L) {
    int N = L*L;
    for(int j = 0; j != N-1; j++) {
        double dx = ( (double)rand() )/RAND_MAX  - delta/2;
        double dy = ( (double)rand() )/RAND_MAX  - delta/2;
        lattice[j*2] = lattice[j*2] + dx;
        lattice[j*2 + 1] = lattice[j*2 + 1] + dy;
        // Apply PBCs
        if(lattice[j*2] < 0){
            lattice[j*2] += L;
        }
        else if(lattice[j*2] > L){
            lattice[j*2] -= L;
        }
        if(lattice[j*2+1] < 0){
            lattice[j*2+1] += L;
        }
        else if(lattice[j*2+1] > L){
            lattice[j*2+1] -= L;
        }
   }
   return;
}

//--------------------------------- Traditional approach ----------------------------

/*
Compute the variance of the distances between x0 and the particles inside a windows of observation (center in x0, radius R).
Samples = 10000 observations.
Input: Radius of the windows of observation, lattice of particles, lattice size, number of samples.
Output: Variance.
*/
double GetVariance_R(double R, double* lattice, int L, int nsamples){
    double sigma_sq=0;
    double N[nsamples] = {0};
    double N_ = 0;
    double N_sq_ = 0;
    double N_sq[nsamples] = {0};

    // Parallelisation over samples, saving useful information in N and N_sq. 
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < nsamples; i++){
        double X_x0 = ((double) rand() / (RAND_MAX)) * L;
        double Y_x0 = ((double) rand() / (RAND_MAX)) * L;
        
        int points_inside = 0;
        for (size_t j = 0; j < L*L; j++){
            if(GetDist_x0_i(X_x0, Y_x0, j, lattice, L) < R){
                points_inside ++;
            }
        }
        N[i] = ((double) points_inside)/nsamples;
        N_sq[i] = ((double) points_inside)*((double) points_inside)/nsamples;
    }
    for (size_t i = 0; i < nsamples; i++){
        N_ += N[i];
        N_sq_ += N_sq[i];
    }
    sigma_sq = N_sq_ - N_*N_;

    // Consistency check
    if(sigma_sq<0){
        cout<<"sigma<0: "<<R<<endl<<" N: "<<N<<endl<<" N_sq: "<<N_sq<<endl;
    }
    return sigma_sq;
}

/*
Compute the variances of the position of particles in the lattice w.r.t. observation windows of varying radii.
Input: lattice, lattice size, output file name, number of samples
Output: variance vs radius dataset at fixed lattice size
*/
void GetVariance(double* lattice,int lattice_size, string output_file,int nsamples){
    // Choosine radii between 0.01 and half the lattice size according to an exp distribution.
    double radii[200];
    radii[0] = 0.01;
    radii[199] = lattice_size/2;
    double c = exp(log(radii[199]/radii[0])/(200-1));
    for (size_t i = 1; i < 200-1; i++){
        radii[i] = radii[i-1]*c;
    }

    ofstream output(output_file);
    for (auto &&i : radii){
        double sigma_sq;        
        sigma_sq = GetVariance_R(i, lattice, lattice_size, nsamples)/(i*i);
        output<<i<<" "<<sigma_sq<<endl;
    }
    output.close();
}

//--------------------------------- Optimised approach ----------------------------

/*
Compute the number of particle inside the observation window keeping the center x0 fixed and varying the radius R.
Input: x0 position, lattice, values of the radius, lattice size.
Output: number of point around x0 as a function of the radius of the observation window.
*/
int* GetN_x0(double X_x0, double Y_x0, double* lattice, double* radii, int L){
    double sigma_sq;
    int *N_x0 = new int[200];
    double distances[L*L];
    // Compute distances bewtween the particles and x0
    for (size_t j = 0; j < L*L; j++){
        distances[j] = (GetDist_x0_i(X_x0, Y_x0, j, lattice, L));
    }
    sort(distances, distances + L*L);
    int new_first = 0;
    int points_inside;
    // Count the number of particles inside each observation window, depending on its radius.
    for (int i = 0; i<200; i++){
        points_inside = 0;
        for (int j = new_first; j < L*L; j++){
            if(distances[j] > radii[i]){
                points_inside = max(j - 1, 0);
                new_first = j;
                break;
            }
        }
        N_x0[i] = points_inside;
    }
    return N_x0;
}

/*
Compute the variances of the position of particles in the lattice w.r.t. observation windows of varying radii.
Input: lattice, lattice size, output file name
Output: variance vs radius dataset at fixed lattice size.
*/
void GetVarianceOpt(double* lattice,int lattice_size, string output_file, int nsamples){
    int **N_x0 = new int*[nsamples];
    for(int i = 0; i<nsamples;i++){
        N_x0[i] = new int[200];
    }
    double N[200] = {0};
    double N_sq[200] = {0};
    double sigma;
    ofstream output(output_file);
    // Choosine radii between 0.01 and half the lattice size according to an exp distribution.
    double radii[200];
    radii[0] = 0.01;
    radii[199] = lattice_size/2;
    double c = exp(log(radii[199]/radii[0])/(200-1));
    for (size_t i = 1; i < 200-1; i++){
        radii[i] = radii[i-1]*c;
    }

    // Parallelisation over samples, saving useful information in N and N_sq 
    #pragma omp parallel for default(shared)
    for (int i = 0; i < nsamples; i++){
        double X_x0 = ((double) rand() / (RAND_MAX)) * lattice_size;
        double Y_x0 = ((double) rand() / (RAND_MAX)) * lattice_size;
        N_x0[i] = GetN_x0(X_x0, Y_x0, lattice, radii, lattice_size);
    }
    // Averages and variance computation    
    for(int i = 0; i<200; i++){
        for(int j = 0; j<nsamples; j++){
            N[i] += ((double) N_x0[j][i])/nsamples;
            N_sq[i] += ((double) N_x0[j][i])*((double) N_x0[j][i])/nsamples;
        }
        sigma = (N_sq[i] - (N[i])*(N[i]))/(radii[i]*radii[i]);
        output<<radii[i]<<" "<<sigma<<endl;
    }
    output.close();
}


