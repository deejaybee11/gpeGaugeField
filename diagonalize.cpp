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
#include "diagonalize.h"
#define _USE_MATH_DEFINES
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include "mkl.h"
#include "mkl_dfti.h"
#include "mkl_lapacke.h"
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
#define Nm 3


int diagonalize(double *pMinusA, MKL_Complex16 *kinEnergy, int count, int maxCount, struct simPars pars)
{
	double H1, H2, H3;
	MKL_INT info, n = Nm, lda = 3, ldvl = 3, ldvr = 3;
	MKL_Complex16 eigenVals[Nm], vl[Nm*Nm], vr[Nm*Nm];
	double *delta;
	delta = (double*)mkl_malloc(pars.nY*sizeof(double), 64);
	int index;
	double min;
	
	int counter[2] = {0,0};
	for(int i = 0; i < pars.nY; ++i)
	{
		delta[i] = ((double)count/(double)maxCount)*pars.detuningGradient*pars.y[i];
	
	}
	
	#pragma omp parallel for		
	for(int i = 0; i < pars.nY; ++i)
	{
		for(int j = 0; j < pars.nX; ++j)
		{
			index = i*pars.nY + j;
			
			H1 = (pow(HBAR,2))*pow((pars.kX[j] + 2*pars.kRaman),2)/(2*pars.mass) - HBAR*delta[i];
			H2 = (pow(HBAR,2))*pow(pars.kX[j],2)/(2*pars.mass);
			H3 = (pow(HBAR,2))*pow((pars.kX[j] - 2*pars.kRaman),2)/(2*pars.mass) + HBAR*delta[i];
			    
			MKL_Complex16 Ht[9] = {{H1,0}, {pars.omegaR/2.0,0}, {0,0}, {pars.omegaR/2.0,0}, {H2,0}, {pars.omegaR/2.0,0}, {0,0}, {pars.omegaR/2.0,0}, {H3,0}};       
				
			info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'N', n, Ht, lda, eigenVals, vl, ldvl, vr, ldvr);
			
			if(eigenVals[0].real < eigenVals[1].real)
			{
				min = std::min(eigenVals[0].real, eigenVals[2].real);
			}
			else
			{
				min = std::min(eigenVals[1].real, eigenVals[2].real);
			}
			
			pMinusA[index] = min;		
			
		} 
	}
	for (int i = 0; i < pars.N; ++i)
	{;
		kinEnergy[i].real = cos((0.5*pars.dt / HBAR) * pMinusA[i]);
		kinEnergy[i].imag = -1*sin((0.5*pars.dt / HBAR) * pMinusA[i]);
	}

	goto cleanup;
		
	cleanup:
		mkl_free(delta);
		if(info != 0) exit(EXIT_FAILURE);
		return 0;
	
}

