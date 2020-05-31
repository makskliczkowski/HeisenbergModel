// Hamiltonian.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
// Eigenvalue calculations for Heisenberg model Hamiltonian


#include "header.h"
#include "Lanczos.h"
#include "testingMains.h"

int main(int argc, char** argv)
{
    try {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //for (int L = 6; L <=12; L = L + 2) {
            int L = std::stoi(argv[1], NULL);
            double delta = std::stod(argv[2], NULL);
            double J = std::stod(argv[3], NULL);
            /*std::ofstream file;
            file.open("sinusy.dat", std::ios::out | std::ios::trunc);
            if (!file.is_open()) {
                throw "can't open a file\n";
            }
            for (int i = 1; i <= static_cast<int>(L * 1.0 / 2.0); i++) {
                double q = i * 2 * std::_Pi / (1.0 * L);
                file << q / std::_Pi << " " << 0.5 * sin(q * 1.0) << " " << sin(q / 2) << std::endl;
            }*/
            //testHeatCapacityLanczos(argc, argv);
            testSpinStructureLanczos(argc, argv);
            //findS_q_omega_map(L, delta, J);
           // findSqw_mapWithTimeEvo(L, delta, J, h);

            /*double delta =0.0;
            const char* meanEnName = (std::to_string(L) + "_d"+ std::to_string(int(delta)) + "_meanEnergy.dat").c_str();
            meanEn.open(meanEnName, std::ios::out);
            //-------define outside scope for parallelism-----------
            arma::mat Hd;
            arma::vec eigval;
            arma::mat eigvec;
            arma::cx_vec psi,phi_T;
            arma::vec psi_0,coeff;
            double Einit = 0.0;
            double Tinit = 0.0;
            int whichBlock = static_cast<int>((1.0 * L) / 2);
            //-------------------------------------------------------
            //-----creating psi_0 as base vector for spins(half up half down)---------------
            arma::vec psi_0_inBase(L,arma::fill::zeros);
            for (int i = 0; i < L; i++) {
                if (i < L / 2) psi_0_inBase[i] = 1;
                else psi_0_inBase[i] = 0;
            }
            psi_0 = constructStateFromBase(L, whichBlock, binaryToState(psi_0_inBase,L));
                psi_0_inBase.clear();//clear memory
            //--------------HAMILTONIAN-------------------------------------------------------
            Hd = hamBlocks(L, delta, whichBlock);
            arma::eig_sym(eigval, eigvec, Hd);
           // eigval.print("eigval");
                Hd.clear();//clear memory
            //-------Einit-------------------------------------------------------------------
            coeff = arma::vec(eigval.n_elem, arma::fill::zeros);
            for (int j = 0; j < eigval.n_elem; j++) {//overlapping of psi(0) and eigenfunctions of H
                coeff[j] = cdot(eigvec.col(j), psi_0);
                Einit = Einit + (eigval[j] *coeff[j]*coeff[j]);
                //Einit = cdot(psi_0, Hd * psi_0);
            }
            Einit = Einit / L;
            Tinit = meanEnergy(L, eigval, meanEn,Einit);
            //std::cout << "Einit=" << Einit << "\t" << "Tinit= " << Tinit << std::endl;
                psi_0.clear();//clear memory
            //----------construct |Phi_T> = e^(-beta_init *En)|n>---------------------------
            phi_T = arma::cx_vec(eigval.n_elem, arma::fill::zeros);
            for (int j = 0; j < eigval.n_elem; j++) {
                phi_T = phi_T + exp(-eigval[j]/DBL_MAX)*eigvec.col(j);
            }
            phi_T = phi_T / eigval.n_elem;
            //phi_T.print("phi_T");
            //arma::cx_vec eig0(eigval.n_elem, arma::fill::zeros);
            //eig0.set_real(eigvec.col(0));
            //arma::cx_vec Szi0 = findMeanS_zMap(L, whichBlock,eig0 );
            //Szi0.print("For eigvec with smallest energy spins on each site:");
            //----------time dep -----------------------------------------------------------
            psi = timeEvo_H_not_depend(L, delta, 0.05, 2*L*1.0, Tinit, eigval, coeff, eigvec, phi_T);
            //-----------------------------------------------------------------------
                eigvec.clear();// clear memory from not used vectors and matrices for next step
                eigval.clear();
                psi.clear();
                phi_T.clear();
                coeff.clear();

            meanEn.close();
        }*/
        //-------------timer---------------------------------------------------------------------
            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            std::cout << "For L = "<< L << " it took: " << time_span.count() << "s" << std::endl;
        //----------------------------------------------------------------------------------------- 

        //}
        return 0;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
    catch (...) {
        std::cout << "something\n";
    }
}
