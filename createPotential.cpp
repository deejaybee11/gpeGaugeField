/*The MIT License(MIT)
	
	Copyright(c)[2015][Dylan Brown]
	
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files(the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions :
	
	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "main.h"
#include "createPotential.h"
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

#if !defined(MKL_ILP64)
	#define LI "%li"
	#else
	#define LI "%lli"
#endif

#define HBAR 1.054e-34

 
void createKineticEnergy(MKL_Complex16 *kinEnergy, struct simPars pars, bool imProp)
{
	int index;
	double hPotVal;
	if(imProp==false)
	{
		for (int i = 0; i < pars.nX; ++i)
		{
			for (int j = 0; j < pars.nY; ++j)
			{	
				index = i*pars.nY + j;
				kinEnergy[index].real = cos((pars.dt / HBAR) * (pow(HBAR,2.0)*(pow(pars.kX[i],2.0) + pow(pars.kY[j],2.0))) / (2.0 * pars.mass));
				kinEnergy[index].imag = -1*sin((pars.dt / HBAR) * (pow(HBAR,2.0)*(pow(pars.kX[i],2.0) + pow(pars.kY[j],2.0))) / (2.0 * pars.mass));
			}
		}
	}
	else
	{
		for (int i = 0; i < pars.nX; ++i)
		{
			for (int j = 0; j < pars.nY; ++j)
			{
				index = i*pars.nY + j;
				kinEnergy[index].real = exp(-1.0*(pars.dt / HBAR) * (pow(HBAR,2.0)*(pow(pars.kX[i],2.0) + pow(pars.kY[j],2.0))) / (2.0 * pars.mass));
				kinEnergy[index].imag = 0;
			}
		}
	}
}

void createKinEnergyY(MKL_Complex16 *kinEnergyY, struct simPars pars)
{
	int index;
	double *saveE = 0;
	saveE = (double*)mkl_malloc(pars.N*sizeof(double),64);
        for (int i = 0; i < pars.nX; ++i)
        {
                for (int j = 0; j < pars.nY; ++j)
                {
                        index = i*pars.nY + j;
                        kinEnergyY[index].real = cos((pars.dt / HBAR) * (pow(HBAR,2.0)*(pow(pars.kY[i],2.0))) / (2.0 * pars.mass));
                        kinEnergyY[index].imag = -1*sin((pars.dt / HBAR) *(pow(HBAR,2.0)*(pow(pars.kY[i],2.0))) / (2.0 * pars.mass));
                }
        }
}


void createNonlinearEnergy(MKL_Complex16 *posPot, MKL_Complex16 *psi, struct simPars pars, bool imProp, bool hPotOn)
{	
	double hPotVal;
	int index;
	if(imProp==true && (hPotOn==true || hPotOn==false))
	{
		for(int ii = 0; ii < pars.nX; ++ii)
		{
			for(int jj = 0; jj < pars.nY; ++jj)
			{
				index = ii*pars.nY + jj;
				hPotVal = (pars.mass/2.0) * (pow(pars.omegaX,2.0)*pow(pars.x[ii],2.0) + pow(pars.omegaY,2.0)*pow(pars.y[jj],2.0));
				posPot[index].real = exp(-1.0 * (pars.intPot * (psi[index].real*psi[index].real + psi[index].imag*psi[index].imag)/pars.nAtoms + hPotVal)* pars.dt / HBAR);
				posPot[index].imag = 0;
			}
		}
	}
	else if(imProp==false && hPotOn==false)
	{
		for(int ii = 0; ii < pars.nX; ++ii)
		{
			for(int jj = 0; jj < pars.nY; ++jj)
			{
				index = ii*pars.nY + jj;
				hPotVal = 0;
				posPot[index].real = cos((pars.intPot * (pars.intPot * (psi[index].real*psi[index].real + psi[index].imag*psi[index].imag) + hPotVal) * pars.dt / HBAR));
				posPot[index].imag = -1*sin((pars.intPot * (pars.intPot * (psi[index].real*psi[index].real + psi[index].imag*psi[index].imag) + hPotVal) * pars.dt / HBAR));
			}
		}
	}
	else
	{
		for(int ii = 0; ii < pars.nX; ++ii)
		{
			for(int jj = 0; jj < pars.nY; ++jj)
			{
				index = ii*pars.nY + jj;
				hPotVal = ((pars.mass/2.0) * (pow(pars.omegaX,2.0)*pow(pars.x[ii],2.0) + pow(pars.omegaY,2.0)*pow(pars.y[jj],2.0)));
				posPot[index].real = cos((pars.intPot * (pars.intPot * (psi[index].real*psi[index].real + psi[index].imag*psi[index].imag) + hPotVal) * pars.dt / HBAR));
				posPot[index].imag = -1*sin((pars.intPot * (pars.intPot * (psi[index].real*psi[index].real + psi[index].imag*psi[index].imag) + hPotVal) * pars.dt / HBAR));
			}
		}
	}
		
}

static void saveArray(double *absPsi, int N, const char * fitsFileName, struct simPars pars, int posSpace)
{	
	/*Saves data as a 2D fits image by summing along one dimension of the array*/
	double *saveFile = 0;
	saveFile = (double*)mkl_malloc(pars.nX*pars.nY*sizeof(double),64);
	if (saveFile==0) printf("Somethings wronnngng\n");
	fitsfile *fptr;
	int status = 0;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {pars.nY, pars.nX};
	
	for(int i = 0; i < pars.N; ++i)
	{
		saveFile[i] = absPsi[i];
	}

	//printf("Starting save routine\n");
	fits_create_file(&fptr, fitsFileName, &status);
	//printf("Create file %s\n", 0==status ? "PASSED" : "FAILED");
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	//printf("Create image %s\n", 0==status ? "PASSED" : "FAILED");
	fits_update_key(fptr, TINT, "NX", &pars.nX, "Number of X points", &status);
	fits_update_key(fptr, TINT, "NY", &pars.nY, "Number of Y points", &status);
	fits_update_key(fptr, TINT, "POSSPACE", &posSpace, "Image is in Position space?", &status);
	//printf("Update keys %s\n", 0==status ? "PASSED" : "FAILED");
	nelements = naxes[0] * naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, saveFile, &status);
	//printf("Write image %s\n", 0==status ? "PASSED" : "FAILED");
	fits_close_file(fptr, &status);
	//printf("Close file %s\n", 0==status ? "PASSED" : "FAILED");
	fits_report_error(stderr, status);
	//printf("Save complete!\n");
	mkl_free(saveFile);


}
