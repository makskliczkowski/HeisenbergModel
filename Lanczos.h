#pragma once
#include "header.h"

void LanczosDiag(int SystemSize, int stepNumber, const arma::mat& Hamiltonian, arma::vec & eigvalLan, arma::mat & eigvecLan);
void LanczosDiagWithVector(int SystemSize, int stepNumber, const arma::mat& Hamiltonian, arma::vec& eigvalLan, arma::mat& lanczosHam, const arma::cx_vec& vector); //starting from a given vector
arma::mat LanczosDiagGiveKrylov(int systemSize, int lanSteps, const arma::mat& Hamiltonian, arma::vec& eigvalLan, arma::mat& eigvecLan); //giving Krylov space vectors
//--------------------STATISTICS--------------------------------
void heatCapLanczos(int L, int r, int lanSize,int systemSize, const arma::mat& Hamiltonian);
void spinSpectralResponseLanczos(int L, int whichBlock, int systemSize, int stepNumber, const arma::mat& Hamiltonian);