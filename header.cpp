#include "header.h"

//=====================================================HAMILTONIANS CONSTRUCTION=====================================================
arma::mat hamStandWay(int L, double delta)
{

	//J=1
	long numberOfStates = pow(2, L);
	arma::mat H(numberOfStates, numberOfStates, arma::fill::zeros);
	for (int i = 0; i < numberOfStates; i++) { //loop over states that we are working on
		arma::vec vector = stateToBinary(i, L);
		for (int j = 0; j < L -1; j++) { //loop over latice poitns -> spin sides
			double s_j =   static_cast<bool>(vector[j]) ? 0.5 : -0.5; //jth spin z value
			double s_jp1 = static_cast<bool>(vector[j+1]) ? 0.5 : -0.5;; // j+1 spin z value

			//Hp - diagonal part
			H(i, i) = H(i, i) + delta * s_j * s_jp1;
			//Hk Si+Sip1-
			if (vector[j] == 0 && vector[j + 1] == 1) {
				arma::vec changedVec = vector;
				int temp = vector[j];
				changedVec[j] = vector[j + 1];
				changedVec[j + 1] = temp;  //change their positions
				int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
				H(changedState,i) += 0.5;
			}
			//Sip1+Si-
			else if (vector[j] == 1 && vector[j + 1] == 0) {
				arma::vec changedVec = vector;
				int temp = vector[j];
				changedVec[j] = vector[j + 1];
				changedVec[j + 1] = temp;  //change their positions
				int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
				H(changedState, i) += 0.5;
			}
		}	
	}
	return H;
}

arma::mat hamBlocks(int L, double delta, int which)
{
	//we will create a vector of matrices that account for the whole matrix
	long numberOfStates = pow(2, L);
	std::vector<int> map = mapTotalSpins(L,which); //sorted vectors of spins

	int size = map.size();
	arma::mat blockMatRef(size, size, arma::fill::zeros);
	//arma::mat* blockMat = new arma::mat(size, size, arma::fill::zeros);//define each block matrix
	//arma::mat& blockMatRef = *blockMat; //must create reference to acces pointer to matrix, probably
	for (int j = 0; j < size; j++) { //loop over each subgroup spins
		arma::vec vector = stateToBinary(map[j], L); // temp vector for each number
		//vector.print("v" + std::to_string(j)+"=");
		for (int k = 0; k < L - 1; k++) { //go over a vector
			double s_j = static_cast<bool>(vector[k]) ? 0.5 : -0.5; //jth spin z value
			double s_jp1 = static_cast<bool>(vector[k + 1]) ? 0.5 : -0.5;; // j+1 spin z value
			//Hp
			blockMatRef(j, j) = (blockMatRef(j, j)) + delta * s_j * s_jp1;
			//Hk
			if (vector[k] == 0 && vector[k + 1] == 1) {
				arma::vec changedVec = vector;
				int temp = vector[k];
				changedVec[k] = vector[k + 1];
				changedVec[k + 1] = temp;  //change their positions
				int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
				int positionOfChanged = binSearch(map, 0, size - 1, changedState); //we find it via binary search and it outputs the position in our sourted array vecForStot
				blockMatRef(positionOfChanged, j) += 0.5; //get it inside
			}
			//Sip1+Si-
			else if (vector[k] == 1 && vector[k + 1] == 0) {
				arma::vec changedVec = vector;
				int temp = vector[k];
				changedVec[k] = vector[k + 1];
				changedVec[k + 1] = temp;  //change their positions
				int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
				int positionOfChanged = binSearch(map, 0, size - 1, changedState); //we find it via binary search and it outputs the position in our sourted array vecForStot
				blockMatRef(positionOfChanged, j) += 0.5; //get it inside
			}

			}
		}
	//blockMatRef.print("MAT= ");
	return blockMatRef;
}

arma::mat hamBlocksPBC(int L, double delta, double J, int which)
{
	//we will create a vector of matrices that account for the whole matrix
	long numberOfStates = pow(2, L);
	std::vector<int> map = mapTotalSpins(L, which); //sorted vectors of spins

	int size = map.size();
	int PBC = 0;
	arma::mat blockMatRef(size, size, arma::fill::zeros);
	//arma::mat* blockMat = new arma::mat(size, size, arma::fill::zeros);//define each block matrix
	//arma::mat& blockMatRef = *blockMat; //must create reference to acces pointer to matrix, probably
	for (int j = 0; j < size; j++) { //loop over each subgroup spins
		arma::vec vector = stateToBinary(map[j], L); // temp vector for each number
		//vector.print("v" + std::to_string(j)+"=");
		for (int k = 0; k < L; k++) { //go over a vector
			if (k == L - 1) PBC = 0;
			else PBC = k + 1;
			double s_j = static_cast<bool>(vector[k]) ? 0.5 : -0.5; //jth spin z value
			double s_jp1 = static_cast<bool>(vector[PBC]) ? 0.5 : -0.5;; // j+1 spin z value
			//Hp
			blockMatRef(j, j) = (blockMatRef(j, j)) +J* delta * s_j * s_jp1;
			//Hk
			if (vector[k] == 0 && vector[PBC] == 1) {
				arma::vec changedVec = vector;
				int temp = vector[k];
				changedVec[k] = vector[PBC];
				changedVec[PBC] = temp;  //change their positions
				int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
				int positionOfChanged = binSearch(map, 0, size - 1, changedState); //we find it via binary search and it outputs the position in our sourted array vecForStot
				blockMatRef(positionOfChanged, j) += J* 0.5; //get it inside
			}
			//Sip1+Si-
			else if (vector[k] == 1 && vector[PBC] == 0) {
				arma::vec changedVec = vector;
				int temp = vector[k];
				changedVec[k] = vector[PBC];
				changedVec[PBC] = temp;  //change their positions
				int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
				int positionOfChanged = binSearch(map, 0, size - 1, changedState); //we find it via binary search and it outputs the position in our sourted array vecForStot
				blockMatRef(positionOfChanged, j) += J*0.5; //get it inside
			}

		}
	}
	//blockMatRef.print("MAT= ");
	return blockMatRef;
}

arma::mat hamBlocksPBC_KuboLinear(int L, double delta, double J, int which, double h,double q)
{
	arma::mat blockMatRef = hamBlocksPBC(L, delta, J, which);//non perturbed Hamiltonian
	int matrixSize = blockMatRef.n_cols;
	for (int m = 0; m < L; m++){
		arma::vec dH = h * cos(q*m) * findSiz(L, which, m);//finding Siz for it will change diagonal elements
		for (int j = 0; j < matrixSize; j++) {
			blockMatRef(j, j) = blockMatRef(j, j) + dH(j);
		}
	}
	return blockMatRef;
}

arma::mat addKuboLinearToBlockHamiltonian(int L, double h,double q, int which, const arma::mat & Hamiltonian)
{
	arma::mat blockMatRef = Hamiltonian;
	int matrixSize = blockMatRef.n_cols;
	for (int m = 0; m < L; m++) {
		arma::vec dH = h * cos(q * (m*1.0)) * findSiz(L, which, m);//finding Siz for it will change diagonal elements
		for (int j = 0; j < matrixSize; j++) {
			blockMatRef(j, j) = blockMatRef(j, j) + dH(j);
		}
	}
	return blockMatRef;
}

std::vector<arma::mat> hamBlocks(int L, double delta)
{

	//we will create a vector of matrices that account for the whole matrix
	long numberOfStates = pow(2, L);
	std::vector<arma::mat> vecOfMat;
	std::vector<std::vector<int>> map = mapTotalSpins(L); //sorted vectors of spins


	for (int i = 0; i < L+1; i++) { // loop over each Stot
		std::vector<int> vecForStot = map[i];
		int size = vecForStot.size();
		arma::mat blockMatRef(size, size, arma::fill::zeros);
		//arma::mat* blockMat = new arma::mat(size, size, arma::fill::zeros);//define each block matrix
		//arma::mat& blockMatRef = *blockMat; //must create reference to acces pointer to matrix, probably
		for (int j = 0; j < size ; j++) { //loop over each subgroup spins
			arma::vec vector = stateToBinary(vecForStot[j], L); // temp vector for each number
			for (int k = 0; k < L - 1; k++) { //go over a vector
				double s_j = static_cast<bool>(vector[k]) ? 0.5 : -0.5; //jth spin z value
				double s_jp1 = static_cast<bool>(vector[k + 1]) ? 0.5 : -0.5;; // j+1 spin z value
				//Hp
				blockMatRef(j, j) = (blockMatRef(j, j)) + delta * s_j * s_jp1;
				//Hk
				if (vector[k] == 0 && vector[k + 1] == 1) {
					arma::vec changedVec = vector;
					int temp = vector[k];
					changedVec[k] = vector[k + 1];
					changedVec[k + 1] = temp;  //change their positions
					int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
					int positionOfChanged = binSearch(vecForStot, 0, size - 1, changedState); //we find it via binary search and it outputs the position in our sourted array vecForStot
					blockMatRef(positionOfChanged, j) += 0.5; //get it inside
				}
				//Sip1+Si-
				else if (vector[k] == 1 && vector[k + 1] == 0) {
					arma::vec changedVec = vector;
					int temp = vector[k];
					changedVec[k] = vector[k + 1];
					changedVec[k + 1] = temp;  //change their positions
					int changedState = binaryToState(changedVec, L); //binary number accoutning for the changed state
					int positionOfChanged = binSearch(vecForStot, 0, size - 1, changedState); //we find it via binary search and it outputs the position in our sourted array vecForStot
					blockMatRef(positionOfChanged, j) += 0.5; //get it inside
				}

			}
		}
		blockMatRef.print("MAT= ");
		vecOfMat.push_back(blockMatRef);
	}
	return vecOfMat;
}

std::vector<std::vector<int>> mapTotalSpins(int L)//this is needed for creating block Hamiltonian
{
	
	std::vector<std::vector<int>> map;
	double ss = 0.0;
	std::vector<int> states;
	for (int i = 0; i < pow(2,L); i++) {//fill vector with numbers from 0 to L-1 to later erease them O(n)
		states.push_back(i);
	}

	for (int j = 0; j < L + 1; j++) { //loop over total spin O(n)
		std::vector<int> tempSubvector;
		double stot = -(1.0*L) / 2 + j;
		int idx = 0;
		for (int i = 0; i < states.size(); i++) { //we make it smaller by erasing elements each time not to check all again
			
			arma::vec baseVec = stateToBinary(states[i], L);
			ss = 0;
			for (int k = 0; k < L; k++) {
				double sz = static_cast<bool>(baseVec[k]) ? 0.5 : -0.5; //check each spin
				ss += sz;
			}
			if (ss == stot) {
				tempSubvector.push_back(states[i]);//if good add to segment
				states.erase(states.begin() + i);
				i = i - 1;
				idx++;
			}

		}
		map.push_back(tempSubvector);
		tempSubvector.clear();
		//std::cout << idx << std::endl;
	}
	return map;
}
	/*--------------------------------------------------------------------------------------------------NOTE------------------------------  ------------------------------------------------------------------       */

	//we could also think of knowing number of Ones each time we make the binary change, that would allow us to map each element to segment again much quicker - > number of ones would be our table element number. 
	//also then getting back to Hamiltonian we would again use binary search to create each submatrix easily, probably...
	
	/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

std::vector<int> mapTotalSpins(int L, int which) { //when we want just one matrix!


	std::vector<int> map;
	double stot;
	if (which < 0 || which > L) {
		if (L % 2 == 0) stot = 0;
		else stot = 0.5;
		 //returns 0 or 1/2 whem someome gives crazy spin number (from 0 to L) ->L+1 totoal spins!
	}
	else {
		stot = -(1.0 * L) / 2 + which; // we choose which submatrix we want
	}
	int statesNum = pow(2, L);
	int idx = 0;
	for (int i = 0; i < statesNum; i++) { //we make it smaller by erasing elements each time not to check all again
		arma::vec baseVec = stateToBinary(i, L);
		double ss = 0;
		for (int k = 0; k < L; k++) {
			double sz = static_cast<bool>(baseVec[k]) ? 0.5 : -0.5; //check each spin
			ss += sz;
		}
		if (ss == stot) {
			map.push_back(i);//if good add to segment
			idx++;
		}

	}
	//	std::cout << idx << std::endl;
	return map;
}
//functions for base creation-------------------------------------------------------------------------
int binSearch(std::vector<int> vec, int l, int r, int num)
{//o(sqrt(N))
	if (l <= r) {
		int mid = (l + r) / 2;
		if (vec[mid] == num)
			return mid;
		if (vec[mid] > num)
			return binSearch(vec, l, mid - 1, num);
		if (vec[mid] < num)
			return binSearch(vec, mid + 1, r, num);
	}
	return -1;
}

int binSearch(std::vector<double> vec, int l, int r, double num)
{
	if (vec.size() == 0) return -1;
	if (l <= r) {
		int mid = (l + r) / 2;
		if (vec[mid] == num)
			return mid;
		if (vec[mid] > num)
			return binSearch(vec, l, mid - 1, num);
		if (vec[mid] < num)
			return binSearch(vec, mid + 1, r, num);
	}
	return -1;
}

arma::vec stateToBinary(int state, int L)
{
	arma::vec temp(L, arma::fill::zeros);
	int rest = state;
	for (int i = 0; i < L; i++) {
		temp[i] = rest % 2;
		rest = static_cast<int>(rest / 2);
	}
	return temp;
}

int binaryToState(arma::vec state, int L)
{
	int temp = 0;
	for (int i = 0; i < L; i++) {
		temp = temp + state[i] * pow(2, i);
	}
	return temp;
}
//===========================================tests==============================================
bool testEigen(const arma::mat& matrix, const arma::vec& vector,double eigen)
{
	using namespace arma;
	vec phi = matrix * vector; 
	double test = cdot(vector,phi);
	if (test > eigen -0.05 && test < eigen + 0.05) return true;
	else return false;
}

bool* isOrthonormal(const arma::vec& vec1, const arma::vec& vec2)
{
	//first place is orthogonality, second is normality

	double a = cdot(vec1, vec2);
	double b = cdot(vec1, vec1);
	bool* temp = new bool[2];
	if (a == 0 && b < 1.01 && b > 0.98)
	{
		temp[0] = true;
		temp[1] = true;
	}
	else if (a == 0 && b > 1.01 || b < 0.98) {
		temp[0] = false;
		temp[1] = true;
	}
	else if (a != 0 && b < 1.01 && b > 0.98) {
		temp[0] = true;
		temp[1] = false;
	}
	else {
		temp[0] = false;
		temp[1] = false;
	}

	return temp;
}

arma::mat normalize(int L, const arma::mat& eigvec)
{
	arma::mat tempMat;
	arma::vec tempVec;
	int size = eigvec.n_elem;
	for (int i = 0; i < size; i++) {
		tempVec = eigvec.col(i);
		double temp = cdot(tempVec , tempVec); 
		tempVec = tempVec / sqrt(temp);
		tempMat.insert_cols(i, tempVec);
	}


	return tempMat;
}

arma::vec constructStateFromBase(int L, int whichBlock, int state)//searches given state in base if there is one
{
	std::vector<int> map = mapTotalSpins(L, whichBlock);
	int size = map.size();
	arma::vec returnState(size, arma::fill::zeros);
	int number = binSearch(map, 0, size - 1, state);
	if (number == -1) return arma::vec();
	else {
		returnState[number] = 1;
		return returnState;
	}
}

//=====================================================================================
double getGap(int L, const arma::vec& eigval) //getting the bandgap
{
	int size = eigval.n_elem;
	double first = std::numeric_limits<double>::max();//find the two smallest elements
	double second = first;
	for (int i = 0; i < size; i++) {
		if (eigval[i] < first) {
			second = first;
			first = eigval[i];
		}
		else if (eigval[i] < second && eigval[i] != first) {
			second = eigval[i];
		}
	}
	return (second - first)/L;//divide by L
}
//==================================Statistical functions===============================
void heatCap(int L, const arma::vec& eigval, const char* filename)
{
	double C = 0;
	int size = eigval.n_elem;
	std::ofstream file;
	file.open(filename, std::ios::out | std::ios::trunc);
	if (!file.is_open()) {
		throw "can't open a file\n";
	}
	int steps = 1000;
	double T_max = 3;
	double step = T_max / steps;
	for (int j = 1; j < steps; j++) {
		double T = 1.0*j *step ;
		double beta = 1.0 / T;
		double Z = 0;
		double meanE = 0;
		double meanE2 = 0;
		for (int i = 0; i < size; i++) {
			double tempo = -beta * eigval[i];
			double piWithoutZ = exp(tempo);
			Z += piWithoutZ;
			meanE += eigval[i] * piWithoutZ;
			meanE2 += eigval[i] * eigval[i] * piWithoutZ;
		}
		meanE = meanE / Z;//divide for pi
		meanE2 = meanE2 / Z; 
		C = (meanE2 - (meanE * meanE)) * beta * beta/L; //divide by L
		file << T << " " << C << std::endl;
	}
	file.close();
}

void meanEnergy(int L, const arma::vec& eigval, std::ofstream &file)
{

	int size = eigval.n_elem;
	int steps = 1000;
	double T_max = 50;
	double step = T_max / steps;
	double T = 0.05;
	for (int i = 1; i < steps; i++) {
		double beta = 1.0 / T;
		double Z = 0;
		double meanE = 0;
		for (int i = 0; i < size; i++) {
			double temp = -beta * eigval[i];
			double piWithoutZ = exp(temp);
			Z += piWithoutZ;
			meanE += eigval[i] * piWithoutZ;
		}
		meanE = meanE / Z /L;
		file << T << "\t\t\t" << meanE << std::endl;
		T = T + step;
	}
}

double meanEnergy(int L, const arma::vec& eigval, std::ostream& file, double EInit) //returns Tinit for a given Einit of thermal state
{
	double TInit = -1.0;//if -1 couldn't find
	int size = eigval.n_elem;
	int steps = 3000;
	double T_max = 100;
	double step = T_max / steps;
	double T = 0.05;
#pragma loop(hint_parallel(8))
	for (int i = 1; i < steps; i++) {
		double beta = 1.0 / T;
		double Z = 0;
		double meanE = 0;
		for (int i = 0; i < size; i++) {
			double temp = -beta * eigval[i];
			double piWithoutZ = exp(temp);
			Z += piWithoutZ;
			meanE += eigval[i] * piWithoutZ;
		}
		meanE = meanE / Z / L;
		if (meanE > EInit - ACCURACY && meanE < EInit + ACCURACY) TInit = T;
		file << T << "\t\t\t" << meanE << std::endl;
		T = T + step;
	}
	return TInit;
}

double find_stat_sum(const arma::vec& eigenVal, double T)
{
	double Z = 0.0;
	double beta = 1.0 / T;
	//if (T > 1e20) return 1.0 / eigenVal.n_elem;
#pragma loop(hint_parallel(8))
	std::for_each(eigenVal.begin(), eigenVal.end(), [&](double n) {Z += exp(-beta * n); });//lambda function
	if (Z == INFINITY) return DBL_MAX;
	if (isnand(Z)) return 1;
	else return Z;
}

//=========================Spin Mapping====================================================================
arma::cx_vec findMeanS_zMap(int L, int whichBlock, const arma::cx_vec& state){
	std::vector<int> states = mapTotalSpins(L, whichBlock);//we need to have our base vectors again for the eigenvector consists of them
	int vecNumber = states.size();
	arma::cx_vec mean(L, arma::fill::zeros);

	for (int j = 0; j < L; j++) {//go over each site, when we are on each, just take the ith element of the multiplication beacause its our base, check coefficient for j'th spin and go
		//#pragma omp parallel for shared(states,mean,state) num_threads(8) 
		for (int i = 0; i < vecNumber; i++) {//over states
			arma::vec vector = stateToBinary(states[i], L);//our ith state in changed base. We changed it to smaller, which represents numbers from 0 to (L|L/2) when creating blocks
			double coeff = 0;
			if (vector[j] == 1.0) { coeff = 0.5; }
			else { coeff = -0.5; }
			//#pragma omp atomic
			mean[j] += coeff * state[i] * conj(state[i]);
		}
	}
	return mean;
	//make a sum over base vectors in form (100...), (01000... ) etc. and they are represented inside by spin vectors from map, then we just check spin on each site from each base vector and add
}

arma::cx_vec findSq_m(int L, double q, int whichBlock, const arma::vec& state)
{ //a function that corresponds to calculating <m|Sq or for a real vectors Sq|m>
	std::vector<int> states = mapTotalSpins(L, whichBlock);//we need to have our base vectors again for the eigenvector consists of them
	int vecNumber = states.size();
	arma::cx_vec returner = arma::cx_vec(state.n_elem, arma::fill::zeros);
	//---------------------------------------------------------------------------try to change to for_each-------------------------------------------------------------------------------
	arma::cx_double exponent = (0, 0);
	arma::cx_vec temp = arma::cx_vec(state.n_elem, arma::fill::zeros);
	for (int j = 0; j <= L; j++) {
		exponent = exp(imaginary * static_cast<double>(j) * q);
		//bad for memory but we won't deal with that big L because of calculation time
#pragma omp parallel for shared(states,exponent,returner,state,q) num_threads(8) 
		for (int i = 0; i < vecNumber; i++) {//over states
			arma::vec vector = stateToBinary(states[i], L);//our ith state in changed base. We changed it to smaller, which represents numbers from 0 to (L|L/2) when creating blocks
			//vector.clear();
			double coeff = 0;
			if (vector[j] == 1.0) { coeff = 0.5; }
			else { coeff = -0.5; }
			temp[i] = exponent * coeff * state[i];//should be conjugate, but the vector is real for that case now
		}
		returner += temp;
	}
	return returner;
}

arma::cx_vec findSq(int L, double q, int whichBlock)
{
	std::vector<int> states = mapTotalSpins(L, whichBlock);//we need to have our base vectors again for the eigenvector consists of them
	int vecNumber = states.size();
	arma::cx_vec returner = arma::cx_vec(states.size(), arma::fill::zeros);
	arma::cx_double exponent = (0, 0);
#pragma omp single
	for (int j = 0; j <= L; j++) {
		exponent = exp(imaginary * static_cast<double>(j) * q);
#pragma omp parallel for shared(states,exponent,returner,q) num_threads(8) 
		for (int i = 0; i < vecNumber; i++) {//over states
			arma::vec vector = stateToBinary(states[i], L);//our ith state in changed base. We changed it to smaller, which represents numbers from 0 to (L|L/2) when creating blocks
			//vector.clear();
			arma::cx_double coeff = 0.0+imaginary*0.0;
			if (vector[j] == 1.0) { coeff = 0.5 + 0.0 * imaginary; }
			else { coeff = -0.5 +0.0*imaginary; }
			returner[i] = returner[i]+ exponent * coeff;
		}
		//returner.print("returner");
	}
	return returner;
}
//creating matrix(diagonal so it's a vector of Sq) -> it is used in S_q_omega_exact(worse for memory but much better speed when keeping it before calculating <m|Sq|n>)
arma::vec findSiz(int L, int whichBlock, double i)//we want to have a matrix(diagonal so vector) corresponding to our base vectors
{
	std::vector<int> states = mapTotalSpins(L, whichBlock);//we need to have our base vectors again for the eigenvector consists of them
	int vecNumber = states.size();
	arma::vec returner = arma::vec(states.size(), arma::fill::zeros);
	#pragma omp parallel for shared(states,returner) num_threads(8) 
	for (int j = 0; j < vecNumber; j++) {
		arma::vec vector = stateToBinary(states[j], L);//our ith state in changed base
		double coeff = 0;
		if (vector[i] == 1.0) { coeff = 0.5; }
		else { coeff = -0.5; }
		returner[j] = coeff;
	}
	return returner;
}

double* S_q_omega(int L, double T1, double T2, double T3, const arma::vec& eigenVal, const arma::mat& eigenVec, double q, double omega, double Z11, double Z22, double Z33)
{	//T1,T2,T3 for different betas in one loop because of calculation time
	int size = eigenVal.n_elem;
	double beta1 = -1.0 / T1;
	double beta2 = -1.0 / T2;
	double beta3 = -1.0 / T3;
	double S_q_om1, S_q_om2, S_q_om3 = 0;
	//#pragma omp flush [omega,q,eigenVec,eigenVal,beta1,beta2,beta3, L,size]
#pragma omp parallel for shared(L,Z11,Z22,Z33,beta1,beta2,beta3,eigenVec,eigenVal,q,omega) num_threads(8) reduction(+: S_q_om1) reduction(+: S_q_om2)  reduction(+: S_q_om3) 
	for (int n = 0; n < size; n++) {
		arma::cx_vec psi_n = findSq_m(L, q, static_cast<int>((1.0 * L) / 2), eigenVec.col(n));
		double En = eigenVal[n];
		double exponent1 = exp(beta1 * En);
		//if (exponent1 == INFINITY || exponent1 > DBL_MAX) exponent1 = DBL_MAX;
		//if (isnand(exponent1)) exponent1 = 1;
		double exponent2 = exp(beta2 * En);
		double exponent3 = exp(beta3 * En);
		for (int m = n; m < size; m++) {//loop over |m>
			double Em = eigenVal[m];
			double tempOmega = -(En - Em);
			if (omega > (tempOmega - ACCURACY / 2) && omega < (tempOmega + ACCURACY / 2)) {//Dirac delta with a certain accuracy-------------
				double omegaReal = -(En - Em);//We are certain that omega is in there!
				//---------------------------------------------------------------------------try to change to for_each-------------------------------------------------------------------------------
				/*arma::cx_vec::iterator psi_n_Iterator = psi_n.begin();
				auto multiplyToGetModule[&](double m_i) {
					double result = m_i * m_i * (*psi_n_Iterator);
				};
				std::for_each(eigenVec.col(m).begin(),eigenVec.col(m).end(),
				 //lambda function for multiplication*/
				arma::cx_double module = cdot(arma::cx_vec(eigenVec.col(m), arma::vec(eigenVec.col(m).n_elem, arma::fill::zeros)), psi_n);
#pragma omp atomic
				S_q_om1 += ((exponent1 * module * conj(module)) / Z11).real();
#pragma omp atomic
				S_q_om2 += ((exponent2 * module * conj(module)) / Z22).real();
#pragma omp atomic
				S_q_om3 += ((exponent3 * module * conj(module)) / Z33).real();
			}
		}
		//S_q_om1 =  S_q_om1;
		/*if (S_q_om1 == INFINITY || S_q_om1 > 1e30)
		{
			S_q_om1 = DBL_MAX;
		}
		//S_q_om2 = ;
		if (S_q_om2 == INFINITY || S_q_om2 > 1e30)
		{
			S_q_om2 = DBL_MAX;
		}
		//S_q_om3 = S_q_om3 * exponent3;
		if (S_q_om3 == INFINITY || S_q_om3 > 1e30)
		{
			S_q_om3 = DBL_MAX;
		}*/
	}

	static double* returner = new double[3];
	returner[0] = S_q_om1;
	returner[1] = S_q_om2;
	returner[2] = S_q_om3;
	return returner;
}
//rather not use, better version bellow
std::vector<std::vector<double>> S_q_omega_exact(int L, double q, double T2, const arma::vec& eigenVal, const arma::mat& eigenVec, double dw, double Z22, double Z33) //fastest working version
{
	int size = eigenVal.n_elem;
	double beta2 = -1.0 / T2;
	int omeganumber = static_cast<int>(((eigenVal[size - 1] - eigenVal[0]) / dw)+dw/2);
	std::vector<double> S_q_om1, S_q_om2, S_q_om3;
	for (int a = 0; a < omeganumber; a++) { //clear for new q;
		S_q_om1.push_back(0.0);
		S_q_om2.push_back(0.0);
		S_q_om3.push_back(0.0);
	}
	arma::cx_vec Sq = arma::cx_vec(eigenVec.n_elem, arma::fill::zeros);
	Sq = findSq(L, q, static_cast<int>((1.0 * L) / 2)); // we leave it like that with no enourmous calculation speed down but some memory loss prevention
#pragma omp parallel for shared(Sq,S_q_om1,S_q_om2,S_q_om3,Z22,Z33,beta2,eigenVec,eigenVal,q,omeganumber,size) num_threads(8)
	for (int n = 0; n < size; n++) {
		double En = eigenVal[n];
		//double exponent1 = exp(beta1 * En); no exponent for T=0;
		double exponent2 = exp(beta2 * En);
		//double exponent3 = 1; // for infinite temp
		//for T=0 we can go outside, as statistic coefficien is just a delta function with E_0-E_n;
		arma::cx_double module0 = (0, 0);
		for (int a = 0; a < eigenVec.col(n).n_elem; a++) {
			module0 = module0 + (eigenVec.col(0)[a]) * (eigenVec.col(n)[a]) * Sq[a];
		}
		int position0 = static_cast<int>(((En - eigenVal[0]) / dw)+dw/2);
		//if (position0 == 0) std::cout << module0 << "\n";
#pragma omp atomic
		S_q_om1[position0] = S_q_om1[position0]+((module0 * conj(module0))).real();
		//--------------------------------------------------------------------------------------
		for (int m = n; m < size; m++) {//loop over |m>
			double Em = eigenVal[m];
			double tempOmega = -(En - Em);
			int position = (static_cast<int>((tempOmega / dw)));
			arma::cx_double module = (0, 0);
			for (int a = 0; a < eigenVec.col(m).n_elem; a++) {//<m|Sq|n>
				module = module + (eigenVec.col(m)[a]) * (eigenVec.col(n)[a]) * Sq[a];
			}
#pragma omp atomic
			S_q_om2[position] = S_q_om2[position]+((exponent2 * module * conj(module)) /( Z22)).real() ;
#pragma omp atomic
			S_q_om3[position] = S_q_om3[position]+((module * conj(module)) /( Z33)).real() ;
		}
	}
//#pragma barrier
	std::vector<std::vector<double>> returner(3, std::vector<double>(omeganumber));
	returner[0] = (S_q_om1);
	returner[1] = (S_q_om2);
	returner[2] = (S_q_om3);
	return returner;
}

arma::cx_double SqFromFourierForGivenQ_timeInstance(int L, double q, double dt, const arma::cx_vec meanSz)//giving \sum_l exp(iql)<Slz(t)>
{
	arma::cx_double returner = 0.0 + imaginary*0.0;
	for (int i = 0; i < L; i++) {
		returner += exp(imaginary * q * (i * 1.0)) * meanSz[i];
	}
	return returner;
}

//==========================Time evolution=================================================================
//method that assumes that Hamiltonian is not time dependent, also it calculates thermal state so better use the 2nd one "timeEvolutionHNotDependent", which assumes nothing.
arma::cx_vec timeEvo_H_not_depend_forManyTest(int L, double delta, double step, double timeMax, double Tinit ,const arma::vec& eigenVal, const arma::vec & coeff, const arma::mat& eigenVec, const arma::cx_vec phi_T)
{
	const char* Spinsname = (std::to_string(L) + "_d" + std::to_string(int(delta)) + "_spinMap.dat").c_str();
	const char* TimeEvoName = (std::to_string(L) + "_d" + std::to_string(int(delta)) + "_timeEvo.dat").c_str();
	std::ofstream Spins(Spinsname, std::ios::out);
	if (!Spins.is_open()) throw "not opened\n";
	std::ofstream TimeEvo(TimeEvoName, std::ios::out);
	if (!TimeEvo.is_open()) throw "not opened\n";

	const int size = coeff.n_elem;

	arma::cx_vec Psi_t(size, arma::fill::zeros);//for time dependent function psi(t)

	TimeEvo << "t" << "\t\t\t" << "Re(<psi(t)|phiT>)" << "\t\t\t" << "Im(<psi(t)|phiT>)" << std::endl;
	int stepNum = static_cast<int>(timeMax / step);
	double t = 0.0;
	for (int i = 0; i < stepNum; i++) {
#pragma loop(hint_parallel(8))
//#pragma omp parallel shared(Psi_t)
		for (int j = 0; j < size; j++) {
			//------------------make function-----------------------
			arma::cx_double exponential;
			exponential = exp(-imaginary * t * eigenVal[j]);//arma::cx_double(cos(-t * eigenVal[j]), sin(-t * eigenVal[j]));
//#pragma omp critical
			Psi_t = Psi_t + (coeff[j] * exponential * eigenVec.col(j));//exponent changed to complex

		}
		//Psi_t.print("psi(t)");
		//overlapping----------------------------------------
		if (delta == 0) {
			arma::cx_double overlapping(0, 0);
#pragma loop(hint_parallel(8))
			overlapping = cdot(Psi_t, phi_T);
			/*for (int j = 0; j < size; j++) {
				arma::cx_double exponential = arma::cx_double(cos(t * eigenVal[j]), sin(t * eigenVal[j]));
				overlapping = overlapping + exponential * coeff[j];
				}*/
			TimeEvo << t << "\t\t\t" << overlapping.real() << "\t\t\t" << overlapping.imag() << "\t\t\t" << abs(overlapping) << std::endl;
		}
		//---------------map spins--------------------------
		arma::cx_vec SpinMap = findMeanS_zMap(L, static_cast<int>((1.0 * L) / 2), Psi_t);
		for (int k = 0; k < L; k++) {//can be put into spinMap function, but it's O(L) so not a big problem
			Spins << k << "\t\t\t" << t << "\t\t\t" << SpinMap[k].real() << std::endl;
		}
		Spins << std::endl;
		//---------------------------------------------------
		t += step;
		Psi_t=arma::cx_vec(size, arma::fill::zeros);
	}
	Spins.close();
	TimeEvo.close();
	return Psi_t;

}

arma::cx_vec timeEvolutionHNotDependent(int L ,double q,double delta, double step, double timeMax, const arma::vec& eigenVal, const arma::vec& coeff, const arma::mat& eigenVec)
{
	//coeff is a vector that contains of overlapps of a given wavefunction with Hamiltonian eigenvectors
	std::ostringstream forq;
	//forq << std::fixed << std::setprecision(1) << q/std::_Pi;
	//const char* fileForSpinsMapLattice_name = (std::to_string(L) + "_d" + std::to_string(int(delta)) +"_q"+ forq.str() + "_spinMapLattice.dat").c_str();//For first method that has been used before
	//std::ofstream fileForSpinsMapLattice(fileForSpinsMapLattice_name, std::ios::out);
	//if (!fileForSpinsMapLattice.is_open()) throw "not opened\n";

	const int sizeOfCoeff = coeff.n_elem;
	arma::cx_vec Psi_t(sizeOfCoeff, arma::fill::zeros);//for time dependent function psi(t)
	int stepNum = static_cast<int>(timeMax / step);
	arma::cx_vec Sq_t(stepNum, arma::fill::zeros);
	double t = 0.0;
	for (int i = 0; i < stepNum; i++) {
		//#pragma omp parallel for shared(Psi_t,eigenVal)
		for (int j = 0; j < sizeOfCoeff; j++) {
			arma::cx_double exponential = exp(-1.0*imaginary * t * eigenVal[j]);
			//------------------make function-----------------------
			//#pragma omp atomic
			Psi_t = Psi_t + (coeff[j] * exponential * eigenVec.col(j));
		}
		//Psi_t.print("Psi_t=");
		//---------------map spins----------------------------------------------------------------------------------
		arma::cx_vec SpinMap = findMeanS_zMap(L, static_cast<int>((1.0 * L) / 2), Psi_t);//we find map for each psi_t
		//for (int k = 0; k < L; k++) {//can be put into spinMap function, but it's O(L) so not a big problem
		//	fileForSpinsMapLattice << k << "\t\t\t" << t << "\t\t\t" << SpinMap[k].real() << std::endl;
		//}
		//fileForSpinsMapLattice << std::endl;
		Sq_t[i] = SqFromFourierForGivenQ_timeInstance(L, q, step, SpinMap);
		//------------------------------------------------------------------------------------------------------------
		t += step;
		Psi_t = arma::cx_vec(sizeOfCoeff, arma::fill::zeros);
	}
	//fileForSpinsMapLattice.close();
	Psi_t.clear();
	return Sq_t;
}


