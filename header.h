#pragma once
#include <iostream>
#include <vector>
#include <armadillo>
#include <stdlib.h>
#include <omp.h>
#include <iomanip>
#include <execution>
#include <chrono>
#include <complex>
#include <algorithm>

constexpr double ACCURACY = 1.0e-2; //for comparing
constexpr std::complex < double> imaginary(0.0, 1.0);

//------------------------------------------------------------HAMILTONIAN--------------------------------------------------------------------------
arma::mat hamStandWay(int L, double delta); //create whole Hamiltonian
arma::mat hamBlocks(int L, double delta, int which); //create block matrix for given total magnetization  block
arma::mat hamBlocksPBC(int L, double delta, double J, int which);//create Hamiltonian for given magnetization block
arma::mat hamBlocksPBC_KuboLinear(int L, double delta, double J, int which, double h,double q);
std::vector<arma::mat> hamBlocks(int L, double delta); //create vector of dense matrices that create Hamiltonian for each totat spin
arma::mat addKuboLinearToBlockHamiltonian(int L, double h, double q, int which, const arma::mat& Hamiltonian);//adding hcos(qi)*Siz \sum to hamiltonian -> Kubo!
//--------------------------------------------------------for blocks--------------------------------------------------------------------------
std::vector<std::vector<int>> mapTotalSpins(int L); //maps all total spins and puts vector of vectors with positions from 0 to 2^L
std::vector<int> mapTotalSpins(int L,int which); //maps only user choosen total spin, so we can just get one matrix
//---------------------------------------------for base creation-------------------------------------------------------------------------------
int binSearch(std::vector<int> vec, int l, int r, int num);
int binSearch(std::vector<int> vec, int l, int r, double num);
arma::vec stateToBinary(int state, int L);
int binaryToState(arma::vec state, int L);
//----------------------------------tests------------------------------------------------------------------------------------------------------------------
bool testEigen(const arma::mat& matrix, const arma::vec& vector, double eigen); //tests if we get true eigenvalues from the eigen alghorithm
bool* isOrthonormal(const arma::vec &vec1, const arma::vec& vec2);//check orthonormality
arma::mat normalize(int L, const arma::mat &eigvec);//normalize vectors
arma::vec constructStateFromBase(int L, int whichBlock, int state); // returns state from base vectors |spinsOnSite> to work with Block Hamiltonian
//----------------------------------------Spin maps-------------------------------------------------------------------------------------
arma::cx_vec findMeanS_zMap(int L, int whichBlock, const arma::cx_vec& state);//which is for blocks
arma::vec findSiz(int L, int whichBlock, double i); //find Siz on ith lattice site
//=========for SpinStructureFactor ==============
arma::cx_vec findSq_m(int L, double q, int whichBlock, const arma::vec& state); //find sum over L e^iqi Siz<n| we want to create it for conjugation as Siz is hermitan and it will speed up the calculations
arma::cx_vec findSq(int L, double q, int whichBlock); //find sum over L e^iqi Siz<n| we want to create it for conjugation as Siz is hermitan and it will speed up the calculations
double* S_q_omega(int L, double T1, double T2, double T3, const arma::vec& eigenVal, const arma::mat& eigenVec, double q, double omega,double Z11, double Z22, double Z33);
std::vector<std::vector<double>> S_q_omega_exact(int L,double q, double T2, const arma::vec& eigenVal, const arma::mat& eigenVec, double dw, double Z22, double Z33);//use this rather than S_q_omega
arma::cx_double SqFromFourierForGivenQ_timeInstance(int L, double q, double dt, const arma::cx_vec meanSz); //creating a map for given q for Spin Structure Factor obtained from double Fourier transform after evolution

//-------------------------------------------time evolution---------------------------------------------------------------------------
arma::cx_vec timeEvo_H_not_depend_forManyTest(int L, double delta, double step, double timeMax, double Tinit, const arma::vec& eigenVal, const arma::vec& coeff, const arma::mat& eigenVec, const arma::cx_vec phi_T);//we send eigenvalues, psi, and coefficients <n|psi(0)> for each n, phi_T is for checking infinite time
arma::cx_vec timeEvolutionHNotDependent(int L, double q, double delta, double step, double timeMax, const arma::vec& eigenVal, const arma::vec& coeff, const arma::mat& eigenVec);//shorter version of above because we won't make termal function
//-------------------------------------statistics----------------------------------------------------------------------------
double getGap(int L, const arma::vec& eigval);//printers
void heatCap(int L, const arma::vec& eigval, const char* filename);//heat capacity
void meanEnergy(int L, const arma::vec& eigval, std::ofstream& file);
double meanEnergy(int L, const arma::vec& eigval, std::ostream& file, double EInit);//return temperature for approximately Einit
double find_stat_sum(const arma::vec& eigenVal, double T);

