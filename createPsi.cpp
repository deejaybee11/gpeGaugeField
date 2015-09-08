#include "main.h"
#include "createPsi.h"
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
//#include "mkl_vml_functions.h"
#include "fitsio.h"
#include<iomanip>

#if !defined(MKL_ILP64)
	#define LI "%li"
	#else
	#define LI "%lli"
#endif

#define HBAR 1.054e-34



/*Function to initialize values in psi array*/
void init(MKL_Complex16 *psi, struct simPars pars)
{
	int index;
	printf("Creating a %s for a %s evolution.\n", (pars.imProp == true) ? "circle" : "Gaussian", (pars.imProp == true) ? "imaginary time" : "real time");
	/*If imProp is true, creates a sphere which is minimised to take on the shape of the trap when evolved in imaginary time. If imProp is false, creates a 3D Gaussian. The Gaussian is a superposition of two condensates, one with "zero" momentum and one with "2*hbar*k" momentum. k is taken as the wavenumber of 780nm light.*/
	for (int i = 0; i < pars.nX; ++i)
	{
		for (int j = 0; j < pars.nY; ++j)
		{
			int index = i*pars.nY + j;
			if(!pars.imProp){
			
				psi[index].real = exp(-pow(pars.x[i],2.0)/pow(pars.sigmaX,2.0) - pow(pars.y[j],2.0)/pow(pars.sigmaY,2.0));
			
				psi[index].imag = 0;
			}
			else
			{
				if(pow((pow(pars.x[i],2.0) + pow(pars.y[j],2.0)),0.5) <= 8e-6)
				{
					psi[index].real = 1;
				}
				else
				{
					psi[index].real = 0;
				}
				psi[index].imag = 0;
			}
		}
	}
}
