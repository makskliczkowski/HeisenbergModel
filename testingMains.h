#pragma once
#include "Lanczos.h"
#include "header.h"

//------------------------------------------------------------LANCZOSES----------------------------------------
void testingLanczos(int argc, char** argv);
void testHeatCapacityLanczos(int argc, char** argv);
void testSpinStructureLanczos(int argc, char** argv);

//========================================================================MAINs FOR TESTING FUNCTIONS AND PLOTTING===========================================================================
void findS_q_omega_map(int L, double delta, double J);//for testing SpinStructureFactor without time evolution
void findSqw_mapWithTimeEvo(int L, double delta, double J, double h);//testing linear response Kubo
