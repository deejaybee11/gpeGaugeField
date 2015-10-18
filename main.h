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


#include<math.h>
#include<stdio.h>
#include<iostream> 
#include "mkl.h"
#include "mkl_dfti.h"
#include<float.h>

static void saveArray(double *absPsi, int N, const char * fitsFileName, struct simPars pars, int posSpace);
void savefile(MKL_Complex16 *psi, const char *filename, struct simPars pars);
void loadfile(MKL_Complex16 *psi, int n, const char * filename);
void savefile(double *array, const char *filename, int size);

struct simPars
{
	int nX;							//length of x array
	int nY;							//length of y array
	int N;							//Total number of points
	double sigmaX;						//Width of initial Gaussian in x
	double sigmaY;						//Width of initial Gaussian in y
	double xGridLength;					//Physical length of x array in meters
	double yGridLength;					//Physical length of y array in meters
	double dx;						//Mesh size
	double dy;						//Mesh size
	double dt;						//Time discretization
	int nSteps;						//Number of real iterations of algorithm
	int iSteps;						//Number of imaginary iterations of algorithm
	double aScatt;						//Rb-87 s wave scattering length
	double mass;						//Mass of Rubidium-87
	double omegaX;						//Harmonic Trap Frequency in x
	double omegaY;						//Harmonic Trap Frequency in y
	double *x;						//Array holding discrete x position values
	double *y;						//Array holding discrete y position values
	double *kX;						//Array holding discrete x wavenumber values
	double *kY;						//Array holding discrete y wavenumber values
	double intPot;						//Strength of nonlinearity
	int nAtoms;						//Number of atoms in condensate
	bool imProp;
	double kRaman;						//Recoil momentum of two-photon transition
	double eRec;						//recoil energy
	double detuningGradient;				//Gradient of detuning
	double omegaR;						//Coupling strength
	
};
