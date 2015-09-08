#define _USE_MATH_DEFINES
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include "mkl.h"
#include "mkl_dfti.h"
#include<float.h>
#include<fstream>
#include "fitsio.h"
#include<iomanip>



void createKineticEnergy(MKL_Complex16 *kinEnergy, struct simPars pars, bool imProp);
void createKinEnergyY(MKL_Complex16 *kinEnergyY, struct simPars pars);
void createNonlinearEnergy(MKL_Complex16 *posPot, MKL_Complex16 *psi, struct simPars pars, bool imProp, bool hPotOn);


