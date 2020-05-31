#include "testingMains.h"
//----------------------------------------------------------LANCZOSES---------------------------------------------

void testingLanczos(int argc, char** argv)
{
	std::ofstream LanczosSteps, LanczosDivergence;
	LanczosSteps.open("lansteps.dat", std::ios::out);
	LanczosDivergence.open("landiv.dat", std::ios::out);
	if (!LanczosSteps.is_open()) throw"couldn't open LanczosSteps";
	if (!LanczosDivergence.is_open()) throw"couldn't open LanczosDivergence";
	//-----------arguments-------------------
	//1 - L
	//2 - delta
	//3 - J
	if (argc == 4) {
		int L = std::stoi(argv[1], NULL);
		double delta = std::stod(argv[2], NULL);
		double J = std::stod(argv[3], NULL);
		int whichBlock = static_cast<int>((1.0 * L) / 2);
		arma::mat Hamiltonian = hamBlocks(L, delta,whichBlock);
		arma::vec eigValuesStandard;
		arma::vec eigValuesLanczos;
		arma::mat eigVectorsStandard;
		arma::mat eigVectorsLanczos;
		//arma::mat eigvectors;
		arma::eig_sym(eigValuesStandard, eigVectorsStandard, Hamiltonian);
		//eigValuesStandard.print("eigenvalues_REAL");
		int systemSize = Hamiltonian.n_cols;
		std::cout << systemSize << std::endl;
		for (int steps = 3; steps < systemSize; steps++) {
			LanczosDiag(systemSize, steps, Hamiltonian,eigValuesLanczos,eigVectorsLanczos);
			LanczosDivergence << steps << "\t\t\t" << ((eigValuesLanczos[0] - eigValuesStandard[0]) / eigValuesStandard[0]) << std::endl;
			if (steps < 35) {
				LanczosSteps << steps;
				for (int i = 0; i < eigValuesLanczos.n_elem; i++) {
					LanczosSteps << " " << eigValuesLanczos[i];
					
				}
				LanczosSteps << std::endl;
			}
			eigValuesLanczos.clear();
		}

	}
	LanczosSteps.close();
	LanczosDivergence.close();

} //convergence

void testHeatCapacityLanczos(int argc, char** argv) {
	//-----------arguments-------------------
	//1 - L
	//2 - delta
	//3 - J
	//4 - Lanczos size
	//5 - lanczon evaluation steps for capacity 
	if (argc == 6) {
		int L = std::stoi(argv[1], NULL);
		double delta = std::stod(argv[2], NULL);
		double J = std::stod(argv[3], NULL);
		int lanSize = std::stoi(argv[4], NULL);
		int r = std::stoi(argv[5], NULL);//how many Lanczos procedures
		int whichBlock = static_cast<int>((1.0 * L) / 2);
		//-----file-------------
		const char* filename = ("Tru_hCL" + std::to_string(L) + "r" + std::to_string(r) + "lanSize" + std::to_string(lanSize) + ".dat").c_str();
		//-----------------------
		arma::mat Hamiltonian = hamBlocks(L, delta, whichBlock);
		arma::vec eigValuesStandard;
		arma::mat eigVectorsStandard;
		arma::eig_sym(eigValuesStandard, Hamiltonian);
		int systemSize = Hamiltonian.n_cols;
		std::cout << systemSize << std::endl;
		heatCapLanczos(L, r, lanSize, systemSize, Hamiltonian);
		heatCap(L, eigValuesStandard, filename);
		//eigValuesStandard.print("Standard=");
	}

}

void testSpinStructureLanczos(int argc, char** argv) {
	//-----------arguments-------------------
	//1 - L
	//2 - delta
	//3 - J
	//4 - Lanczos size
	if (argc == 5) {
		int L = std::stoi(argv[1], NULL);
		double delta = std::stod(argv[2], NULL);
		double J = std::stod(argv[3], NULL);
		int lanSize = std::stoi(argv[4], NULL);
		int whichBlock = static_cast<int>((1.0 * L) / 2);
		arma::mat Hamiltonian = hamBlocksPBC(L, delta,J, whichBlock);//with PBC?
		//arma::vec eigValuesStandard;
		//arma::mat eigVectorsStandard;
		//arma::eig_sym(eigValuesStandard,eigVectorsStandard, Hamiltonian);
		int systemSize = Hamiltonian.n_cols;
		std::cout << systemSize << std::endl;
		//eigVectorsStandard.col(0).print("real groundVec = ");
		spinSpectralResponseLanczos(L, whichBlock, systemSize, lanSize, Hamiltonian);

	}
	return;


}
//========================================================================MAINs FOR TESTING FUNCTIONS AND PLOTTING===========================================================================
void findS_q_omega_map(int L, double delta, double J)
{
	//--------------------------------files------------------------------------------
	std::ofstream MapT1, MapT2, MapT3;
	const char* Map = ("Tt0_L" + std::to_string(L) + "_J" + std::to_string(int(J)) + "_d" + std::to_string(int(delta)) + ".dat").c_str();
	MapT1.open(Map, std::ios::out);//T=0
	Map = ("Tt1_L" + std::to_string(L) + "_J" + std::to_string(int(J)) + "_d" + std::to_string(int(delta)) + ".dat").c_str();
	MapT2.open(Map, std::ios::out);//T=1
	Map = ("T00_L" + std::to_string(L) + "_J" + std::to_string(int(J)) + "_d" + std::to_string(int(delta)) + ".dat").c_str();
	MapT3.open(Map, std::ios::out);//T=infinite
	if (!MapT1.is_open()) throw"couldn't open MapT1";
	if (!MapT2.is_open()) throw"couldn't open MapT2";
	if (!MapT3.is_open()) throw"couldn't open MapT3";
	//MapT1 << "q" << "\t\t" << "w[in pi units]" << "\t\t" << "S(q,w)" << std::endl;
	//MapT2 << "q" << "\t\t" << "w[in pi units]" << "\t\t" << "S(q,w)" << std::endl;
	//MapT3 << "q" << "\t\t" << "w[in pi units]" << "\t\t" << "S(q,w)" << std::endl;
	//-------------------------------------------------------------------------------
	int whichBlock = static_cast<int>((1.0 * L) / 2);
	//-------------------------------Hamiltonian-------------------------------------
	arma::mat Hd = hamBlocksPBC(L, delta, J, whichBlock);
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, Hd);
	//eigval.print("eigval=");
	Hd.clear();
	//double step = 0.01;
	//int stepNumber = static_cast<int>(1 * std::_Pi / step);//for omega steps
	//std::cout << stepNumber << "\n";
	//-------------------------------loops-------------------------------------------
	//double T1 = 0.1;
	double T2 = 1.0;
	double T3 = DBL_MAX;
	//double Z1 = find_stat_sum(eigval, T1);
	double Z2 = find_stat_sum(eigval, T2);
	double Z3 = find_stat_sum(eigval, T3);
	double dw = 0.05;
	int omeganum = static_cast<int>((eigval[eigval.n_elem - 1] - eigval[0]) / dw);
	/*for (int j = 0; j < stepNumber; j++) {
		double omega = j * step;
		bool shouldPrintEndl=0;
		//auto get_omega = [&](double m) {return static_cast<double>(m - n); };//lambda function for a loop

		if (shouldPrintEndl) {
			MapT1 << std::endl;
			MapT2 << std::endl;
			MapT3 << std::endl;
		}
	}*/
	std::vector<std::vector<double>> omegas;
	for (int i = 0; i <= L; i++) {
		double q = i * 2 * std::_Pi / (1.0 * L);
		omegas = S_q_omega_exact(L, q, 1.0, eigval, eigvec, dw, Z2, Z3);
		for (int j = 0; j < omeganum; j++) {
			double omega = j * dw;
			MapT1 << q / std::_Pi << "\t\t" << omega / std::_Pi << "\t\t" << omegas[0][j] << std::endl;
			MapT2 << q / std::_Pi << "\t\t" << omega / std::_Pi << "\t\t" << omegas[1][j] << std::endl;
			MapT3 << q / std::_Pi << "\t\t" << omega / std::_Pi << "\t\t" << omegas[2][j] << std::endl;
		}
		MapT1 << std::endl;
		MapT2 << std::endl;
		MapT3 << std::endl;
	}

	//omegas.clear();
	//eigvec.clear();
	//eigval.clear();
	MapT1.close(); MapT2.close(); MapT3.close();
}

void findSqw_mapWithTimeEvo(int L, double delta, double J, double h)
{
	//checking relaxation after turning off the perturbation
	const int whichBlock = static_cast<int>((1.0 * L) / 2);
	double timeStep = 0.01;
	double maxTime = 5 * L;
	//file--------------------------------------------------------------------
	const char* fileForSpinsMapFourier_name = (std::to_string(L) + "_d" + std::to_string(int(delta)) + "_spinMapFourier.dat").c_str();//For second method - Fourier transform (double)
	std::ofstream fileForSpinsMapFourier(fileForSpinsMapFourier_name, std::ios::out);
	if (!fileForSpinsMapFourier.is_open()) throw "not opened\n";
	//------------------------------------------------------------------------
	arma::mat non_perturbedH = hamBlocksPBC(L, delta, J, whichBlock);
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, non_perturbedH);
	//counting mean difference between levels
	double sum = 0;
	for (int i = 0; i < eigval.n_elem - 2; i++) {
		sum = sum + (eigval[i + 1] - eigval[i]);
	}
	sum = sum / eigval.n_elem;
	std::cout << sum << std::endl;

	for (int m = 0; m <= L / 2; m++) {
		double q = 2 * std::_Pi * m / L;
		arma::mat perturbedH = addKuboLinearToBlockHamiltonian(L, h, q, whichBlock, non_perturbedH);
		//std::cout << "q=" << q/std::_Pi << std::endl;
		//perturbedH.print("perturbedH");
		arma::vec eigvalPert;
		arma::mat eigvecPert;
		arma::eig_sym(eigvalPert, eigvecPert, perturbedH);
		arma::vec psi_0 = eigvecPert.col(0);//saving the ground state -> turning off the perturbation!
		//cleaning-------------------------------
		eigvalPert.clear();
		eigvecPert.clear();
		perturbedH.clear();
		//overlapping----------------------------
		arma::vec coeff = arma::vec(eigval.n_elem, arma::fill::zeros);
#pragma omp parallel for shared(coeff,eigvec,psi_0) num_threads(8)
		for (int j = 0; j < eigval.n_elem; j++) {//overlapping of psi(0) and eigenfunctions of H
			coeff[j] = cdot(eigvec.col(j), psi_0);
		}
		//time evolution--------------------------
		arma::cx_vec Sq_t_map = timeEvolutionHNotDependent(L, q, delta, timeStep, maxTime, eigval, coeff, eigvec);
		//Sq_t_map.print("a");
		//finding map-----------------------------
		int omegasteps = std::_Pi / ACCURACY;
		for (int j = 0; j < omegasteps; j++) {
			double omega = j * ACCURACY;
			arma::cx_double Sqw = 0.0 + imaginary * 0.0;
			//#pragma omp parallel for shared(Sqw,omega,Sq_t_map, timeStep) num_threads(8)
			for (int k = 0; k < Sq_t_map.n_elem; k++) {
				double t = k * timeStep;
				//#pragma omp atomic
				Sqw = Sqw + (timeStep * exp(imaginary * omega * t)) * Sq_t_map[k];
			}
			fileForSpinsMapFourier << q / std::_Pi << "\t\t" << omega / std::_Pi << "\t\t" << sqrt(Sqw * conj(Sqw)).real() << std::endl;
			omega += ACCURACY;
		}
		fileForSpinsMapFourier << std::endl;
	}
	fileForSpinsMapFourier.close();
}
//for finding map without perturturbation for T=0,1,inf
