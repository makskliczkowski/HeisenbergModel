#include "Lanczos.h"
#include <time.h>  

void LanczosDiag(int SystemSize, int stepNumber, const arma::mat& Hamiltonian, arma::vec& eigvalLan, arma::mat& eigvecLan)
{

	//srand(time(NULL));
	arma::arma_rng::set_seed_random();  // set the seed to a random value

	if (stepNumber < 2) return;
	arma::vec Psi0(SystemSize,arma::fill::randu);
	arma::mat pseudoMatrix(stepNumber, stepNumber, arma::fill::zeros);
	//Psi0.print("Psi0=");
	for (int i = 0; i < SystemSize; i++)Psi0(i) = static_cast<double>(Psi0(i)) - 0.5;//create random vector first
	Psi0 = Psi0 / sqrt(cdot(Psi0,Psi0));//normalise
	//Psi0.print("Psi0-0.5I=");
	//BEGINNING---------FIRST ELEMENTS ABOVE LOOP ------------------------------------------------------------
	arma::vec tmp = Hamiltonian * Psi0;//whole vector(can be devided into perpendicular and parallel parts)
	double ai = cdot(Psi0, tmp);
	double bi =0.0;
	arma::vec tmp2 = tmp - ai * Psi0;
	double bip1 = sqrt(cdot(tmp2, tmp2));
	pseudoMatrix(0,0) = ai;
	pseudoMatrix(0, 1) = bip1;
	//--------------------------LOOP-----------------------------------
	for (int i = 1; i < stepNumber-1; i++) {
		//create new ith vector----------
		arma::vec Psii = tmp2 / bip1;
		//put on hamiltonian--------------
		tmp = Hamiltonian * Psii;
		//--------H|phii> = bi|phi_i-1>+ai|phi_i>*bip1|phi_i+1>
		//and next goes getting all coefficients
		ai = cdot(Psii, tmp);
		bi = bip1;// bi is bip1 from last step so we don't have to calculate it again
		tmp2 = tmp - (ai * Psii) - (bi * Psi0);//it will be for new Psii
		bip1 = sqrt(cdot(tmp2, tmp2));
		pseudoMatrix(i, i-1) = bi;
		pseudoMatrix(i, i) = ai;
		pseudoMatrix(i, i+1) = bip1;
		//after writing to matrix we need to change psi_0(which is psi_i-1) to psii
		Psi0 = Psii;
	}
	//ENDING WE PUT LAST CONDITIONS------------------------------------------
	arma::vec Psii = tmp2 / bip1;
	//Psii.print();
	tmp = Hamiltonian * Psii;
	//tmp.print("tmp = ");
	ai = cdot(Psii, tmp);
	pseudoMatrix(stepNumber - 1, stepNumber-2) = bip1;
	pseudoMatrix(stepNumber - 1, stepNumber - 1) = ai;
	//pseudoMatrix.print("Mat=");
	//-----------------------------------------------
	arma::eig_sym(eigvalLan, eigvecLan,pseudoMatrix);
	return;
}

void LanczosDiagWithVector(int SystemSize, int stepNumber, const arma::mat& Hamiltonian, arma::vec& eigvalLan, arma::mat& lanczosHam, const arma::cx_vec & vector)
{

	srand(time(NULL));

	if (stepNumber < 2) return;
	arma::cx_vec Psi0 = arma::normalise(vector); //starting from that
	//imaginary parts are almost zero so for the sake of speed let's make Hamiltonian real
	//Psi0.print("Psi0");
	//BEGINNING---------FIRST ELEMENTS ABOVE LOOP ------------------------------------------------------------
	arma::cx_vec tmp = Hamiltonian * Psi0;//whole vector(can be devided into perpendicular and parallel parts)
	arma::cx_double ai = cdot(Psi0, tmp);
	arma::cx_double bi = 0.0;
	arma::cx_vec tmp2 = tmp - ai * Psi0;
	arma::cx_double bip1 = sqrt(cdot(tmp2, tmp2));
	lanczosHam(0, 0) = ai.real();
	lanczosHam(0, 1) = bip1.real();
	//--------------------------LOOP-----------------------------------
	for (int i = 1; i < stepNumber - 1; i++) {
		//create new ith vector----------
		arma::cx_vec Psii = tmp2 / bip1;
		//put on hamiltonian--------------
		tmp = Hamiltonian * Psii;
		//--------H|phii> = bi|phi_i-1>+ai|phi_i>*bip1|phi_i+1>
		//and next goes getting all coefficients
		ai = cdot(Psii, tmp);
		bi = bip1;// bi is bip1 from last step so we don't have to calculate it again
		tmp2 = tmp - (ai * Psii) - (bi * Psi0);//it will be for new Psii
		bip1 = sqrt(cdot(tmp2, tmp2));
		lanczosHam(i, i - 1) = bi.real();
		lanczosHam(i, i) = ai.real();
		lanczosHam(i, i + 1) = bip1.real();
		//after writing to matrix we need to change psi_0(which is psi_i-1) to psii
		Psi0 = Psii;
	}
	//ENDING WE PUT LAST CONDITIONS------------------------------------------
	arma::cx_vec Psii = tmp2 / bip1;
	//Psii.print();
	tmp = Hamiltonian * Psii;
	//tmp.print("tmp = ");
	ai = cdot(Psii, tmp);
	lanczosHam(stepNumber - 1, stepNumber - 2) = bip1.real();
	lanczosHam(stepNumber - 1, stepNumber - 1) = ai.real();
//	lanczosHam.print("Mat=");
	//-----------------------------------------------
	arma::eig_sym(eigvalLan, lanczosHam);
	return;



}

arma::mat LanczosDiagGiveKrylov(int systemSize, int lanSteps, const arma::mat& Hamiltonian, arma::vec& eigvalLan, arma::mat& eigvecLan) {

	arma::arma_rng::set_seed_random();  // set the seed to a random value
	if (lanSteps < 2 || lanSteps > systemSize) return arma::mat(); //check boundaries
	arma::mat krilovVec(systemSize, lanSteps); //create matrix of vectors to return
	arma::vec Psi0(systemSize, arma::fill::randu);//random vector
	arma::mat pseudoMatrix(lanSteps, lanSteps, arma::fill::zeros); //Lanczos matrix


	for (int i = 0; i < systemSize; i++) Psi0(i) = static_cast<double>(Psi0(i)) - 0.5;//create random vector first
	Psi0 = Psi0 / sqrt(cdot(Psi0, Psi0));//normalise
	krilovVec.col(0) = Psi0; //add vector to matrix zeroth

	//Krilov space {Psi0, H|Psi0> , H H |psi > etc. so we will need change something during calculations, not to calculate this again
	//BEGINNING---------FIRST ELEMENTS ABOVE LOOP ------------------------------------------------------------
	arma::vec tmp = Hamiltonian * Psi0;//whole vector(can be devided into perpendicular and parallel parts)
	double ai = cdot(Psi0, tmp);
	double bi = 0.0;
	arma::vec tmp2 = tmp - ai * Psi0;
	double bip1 = sqrt(cdot(tmp2, tmp2));
	pseudoMatrix(0, 0) = ai;
	pseudoMatrix(0, 1) = bip1;
	//--------------------------LOOP-----------------------------------
	for (int i = 1; i < lanSteps - 1; i++) {
		//create new ith vector----------
		arma::vec Psii = tmp2 / bip1;
		krilovVec.col(i ) = Psii; //add vector to matrix |Psi1>, |Psi2> etc
		//put on hamiltonian--------------
		tmp = Hamiltonian * Psii;
		//--------H|phii> = bi|phi_i-1>+ai|phi_i>*bip1|phi_i+1>
		//and next goes getting all coefficients
		ai = cdot(Psii, tmp);
		bi = bip1;// bi is bip1 from last step so we don't have to calculate it again
		tmp2 = tmp - (ai * Psii) - (bi * Psi0);//it will be for new Psii
		bip1 = sqrt(cdot(tmp2, tmp2));
		pseudoMatrix(i, i - 1) = bi;
		pseudoMatrix(i, i) = ai;
		pseudoMatrix(i, i + 1) = bip1;
		//after writing to matrix we need to change psi_0(which is psi_i-1) to psii
		Psi0 = Psii;
	}
	//ENDING WE PUT LAST CONDITIONS------------------------------------------
	arma::vec Psii = tmp2 / bip1;
	krilovVec.col(lanSteps-1) = Psii; //add vector to matrix first
	//Psii.print();
	tmp = Hamiltonian * Psii;
	//tmp.print("tmp = ");
	ai = cdot(Psii, tmp);
	pseudoMatrix(lanSteps - 1, lanSteps - 2) = bip1;
	pseudoMatrix(lanSteps - 1, lanSteps - 1) = ai;
	//pseudoMatrix.print("Mat=");
	//-----------------------------------------------
	arma::eig_sym(eigvalLan, eigvecLan, pseudoMatrix);
	return krilovVec;


}
//--------------------------------------STATISTICS------------------------------------
void heatCapLanczos(int L,int r,int lanSize,int systemSize, const arma::mat & Hamiltonian)
{
	const char* filename = ("hCL" + std::to_string(L) + "r" + std::to_string(r) + "lanSize" + std::to_string(lanSize) + ".dat").c_str();
	std::ofstream file;
	file.open(filename, std::ios::out | std::ios::trunc);
	if (!file.is_open()) {
		throw "can't open a file\n";
	}
	arma::vec eigvalLan;
	arma::mat eigvecLan;

	double C = 0;
	int steps = 1000;
	double T_max = 3;
	double step = T_max / steps;
	arma::mat overlaps(r,lanSize, arma::fill::zeros);
	for (int n = 0; n < r; n++) {//loop over distinct Lanczos procedures
		LanczosDiag(systemSize, lanSize, Hamiltonian, eigvalLan, eigvecLan);
		#pragma omp parallel for shared(overlaps,eigvecLan) num_threads(8) 
		for (int i = 0; i < lanSize; i++) {
			//eigvalLan.print("Lanczos=");
			overlaps(n, i) = abs(eigvecLan.col(i)[0]) * abs(eigvecLan.col(i)[0]);
		}
	}//first elements that come from the overlapping, as the eigvecs are constructed in Krimov base
		for (int j = 1; j < steps; j++) {//different times
			double T = 1.0*j * step;
			double beta = 1.0 / T;
			double Z = 0;
			double meanE = 0;
			double meanE2 = 0;
			#pragma omp parallel for shared(beta,meanE,meanE2,eigvalLan,Z) num_threads(8) 
			for (int n = 0; n < r; n++) { //loop over distinct Lanczos procedures
				for (int i = 0; i < lanSize; i++) {
					double tempo = -beta * eigvalLan[i];
					double piWithoutZ = exp(tempo) * overlaps(n,i); //we take it as a starting vector in creating the Krilov space, which is also orthogonal
					//, we don't need to calculate the overlap as it's just the first term in this base, because we constructed the hamiltonian from such a base!;
#pragma omp atomic
					Z += piWithoutZ; //adding probablities to get Z
#pragma omp atomic
					meanE = meanE + eigvalLan[i] * piWithoutZ;
#pragma omp atomic
					meanE2 = meanE2 + eigvalLan[i] * eigvalLan[i] * piWithoutZ;
				}
			}
			meanE = meanE / Z;//divide for pi
			meanE2 = meanE2 / Z; 
			C = (meanE2 - (meanE * meanE)) * beta * beta/(L); //divide by L
			if (C < 1e-12 || C==INFINITY) C = 0;
			file << T << " " << C << std::endl;
		}

	file.close();
}

void spinSpectralResponseLanczos(int L,int whichBlock,int systemSize,int stepNumber, const arma::mat& Hamiltonian)
{
	//looking for a Green function representation of S(q,w) for a system with a Lanczos method
	//P. PRELOVŠEK, J. BONÈA, ARXIV: 1111.5931 (2011)
	//E. DAGOTTO, REV. MOD. PHYS. 66, 763 (1994)
	const char* filename = ("SqwLanL" + std::to_string(L) + "lanSize" + std::to_string(stepNumber) + ".dat").c_str();
	std::ofstream file;
	file.open(filename, std::ios::out | std::ios::trunc);
	if (!file.is_open()) {
		throw "can't open a file\n";
	}
	arma::vec eigvalLan;
	arma::mat krilovVec;
	arma::mat eigvecLan;
	arma::vec groundStateLan;
	arma::vec groundState(systemSize,arma::fill::zeros);
	double groundEn;
	krilovVec = LanczosDiagGiveKrylov(systemSize, stepNumber, Hamiltonian, eigvalLan, eigvecLan);
	groundStateLan = eigvecLan.col(0);
	eigvecLan.clear();//clearing for memory seve as we want to calculate only groundstate
	groundEn = eigvalLan(0);
	eigvalLan.clear();//clearing for memory save as we want to calculate only groundstate

	for (int i = 0; i < systemSize; i++) {
		//loop to calculate the ground state from Krilov vectors
		#pragma omp parallel for shared(groundState,groundStateLan,krilovVec) num_threads(8) 
		for (int j = 0; j < stepNumber ; j++) {
			groundState(i) = groundState(i) + groundStateLan(j) * krilovVec.col(j)(i);
		}
	}
	//groundState.print("Numerical ground=");
	arma::mat lanHam(stepNumber, stepNumber, arma::fill::zeros); //create Lanczos matrix
	const double omegaStep = 0.1;
	const double minusOneOverPi = -(1.0 / std::_Pi);
	const int howManyOmegas = static_cast<int>((std::_Pi) / omegaStep);
	for (int i = 1; i <= static_cast<int>(L*1.0 / 2.0); i++) {
		double q = i * 2 * std::_Pi / (1.0 * L);
		arma::cx_vec Sq(systemSize);
		Sq = findSq(L, q, whichBlock);
		//Sq.print("Sq=");
		arma::cx_vec phi0(Sq.n_elem);
		#pragma omp parallel for shared(phi0,groundState,Sq) num_threads(8) 
		for (int j = 0; j < Sq.n_elem; j++) {
			phi0(j) = groundState(j) * Sq(j); //calculating new vector for start(Sq is diagonal)
		}
		//phi0.print("phi0");
		//we need to starn new Lanczos Procces now, starting from phi0
		
		LanczosDiagWithVector(systemSize, stepNumber, Hamiltonian, eigvalLan, lanHam, phi0);
		//eigvalLan.print("eigValLAn");
		arma::cx_double numerator = (cdot(phi0, phi0));
		//loop over the frequencies
		arma::vec Sqw(howManyOmegas);
		#pragma omp parallel for shared(Sqw,lanHam,groundEn,numerator) num_threads(8) 
		for (int j = 1; j <= howManyOmegas; j++)
		{
			double omega = j * omegaStep;
			double eta = ACCURACY;
			arma::cx_double z = omega + imaginary * eta + groundEn;
			arma::cx_double denominator = z - lanHam(stepNumber -1, stepNumber - 1);//last denominator//b_m+1 =0
			for (int k = stepNumber - 1; k > 0; k--) {
				arma::cx_double temp = denominator;
				denominator = z - lanHam(k-1, k-1) - (lanHam.row(k)(k-1) * lanHam.row(k)(k - 1) /temp);
			}//finding the denominator 
			Sqw(j-1) =  minusOneOverPi* ( numerator / denominator).imag();

		}
		for(int j=1;j<= howManyOmegas;j++)	file << q / std::_Pi << " " << j * omegaStep / std::_Pi << " " << Sqw(j-1) << std::endl; //maybe will be faster
		file << std::endl;
	}
	file.close();
}
