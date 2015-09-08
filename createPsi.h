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


void init(MKL_Complex16 *psi, struct simPars pars);

void createSuperposition(MKL_Complex16 *psi, struct simPars pars);
