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
#include "createPsi.h"
#include "createPotential.h"
#include "diagonalize.h"
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

/*
	MKL_Complex16 is a data structure that contains two 64bit doubles, one called real, one called imag. The real and imaginary components need to be evaluated separately by calling
	psi.real and psi.imag respectively. 
	
	Creating an array of MKL_Complex16, eg: MKL_Complex16 *psi[10], creates an array of length 10, each element contains the data structure MKL_Complex16.
	This allows the passing of *psi to a function as a pointer, accessing the elements of psi (eg: psi[4].real) and modifying them in the original memory location without the need for referencing (&)
	or assigning the result to a new array allowing for more efficient memory management (See function "static void init(MKL_Complex16 *psi, struct simPars pars)" for an example).
	
	To view data, run analyze3dCdata.py located in the source folder. This will plot the input and output arrays used in the code.
	
	2D Arrays index in order column-row, eg: arr[col][row]. This translates to a 1D representation where elements are obtained by i*nCols + j where i,j belong to the set [0,nCols), [0,nRows) respectively. 
*/

int main(){
	printf("\033[0m");
	putenv("KMP_BLOCKTIME=infinite");
	putenv("KMP_AFFINITY=verbose,granularity=fine,compact,norespect");
	mkl_set_num_threads(mkl_get_max_threads());
	system("echo KMP_BLOCKTIME = $KMP_BLOCKTIME");
	system("echo KMP_AFFINITY = $KMP_AFFINITY");
	mkl_disable_fast_mm();
	
	/*Structure that holds all simulation parameters*/
	simPars pars;
	int index;
	
	/*Set the number of elements in each dimension*/
	pars.nX = 128;
	pars.nY = 128;
	pars.N = pars.nX*pars.nY;
	
	/*Width of Gaussian in each dimension*/
	pars.sigmaX = 2.5e-6;
	pars.sigmaY = 2.5e-6;
	
	/*Length of each dimension in metres*/
	pars.xGridLength = 20e-6;
	pars.yGridLength = 20e-6;
	
	/*Trap frequencies in each dimension*/
	pars.omegaX = 2*M_PI*150;
	pars.omegaY = 2*M_PI*150;
	
	/*Mass of atomic species*/
	pars.mass = 87 * 1.667e-27;
	
	/*Iteration parameters*/
	pars.nSteps = 1000000000;
	pars.iSteps = 1000000;//500000;
	int maxCount = 200000000; //Number of steps until field is at max strength
	pars.imProp = true;
	pars.dt = 1e-9;
	
	/*Gauge potential parameters*/
	pars.kRaman = 2*M_PI/780e-9;
	pars.eRec = pow(HBAR,2)*pow(pars.kRaman,2)/(2*pars.mass);
	pars.detuningGradient = -0.04*pars.eRec*pars.kRaman/HBAR;
	pars.omegaR = 8.2*pars.eRec;
	
	/*Interaction Potential*/
	pars.nAtoms = 30000;
	pars.aScatt = 5.186e-9;
	double U = (4*M_PI*pow(HBAR,2.0)*pars.aScatt)/pars.mass;
	double lx = pow((HBAR/(pars.omegaX*pars.mass)),0.5);
	pars.intPot = U/lx;
	std::cout << "Intpot = " << pars.intPot << std::endl;
	
	/*Clears all files from fits folder. fitsio throws a tantrum if the file already exists.*/
	system("exec rm -r fits/*.fits");
	
	/*Pointer to psi data*/
	MKL_Complex16 *psi = 0;
	double *pMinusA = 0;
	pMinusA = (double*)mkl_malloc(pars.N*sizeof(double),64);
	
	/*Pointer to |psi| data*/
	double *absPsi = 0;
	MKL_Complex16 *posPot = 0;
	
	/*Int to determine execution status*/
	MKL_LONG status = 0;
	
	/*DFTI descriptor handle as described in the Intel MKL Reference Manual*/
	DFTI_DESCRIPTOR_HANDLE handle = 0, descx = 0, descy = 0;
	
	/*Assign memory locations to x, y arrays*/
	pars.x = 0;
	pars.x = (double*)mkl_malloc(pars.nX*sizeof(double),64);
	if (pars.x == 0) goto failed;
	pars.y = 0;
	pars.y = (double*)mkl_malloc(pars.nY*sizeof(double),64);
	if (pars.y == 0) goto failed;
	
	/*Fill elements of x, y arrays*/
	for (int i = 0; i < pars.nX; ++i)
	{
		pars.x[i] = -1.0 * pars.xGridLength / 2.0 + i*pars.xGridLength / ((double)pars.nX - 1.0);
	}
	for (int i = 0; i < pars.nY; ++i)
	{
		pars.y[i] = -1.0 * pars.yGridLength / 2.0 + i*pars.yGridLength / ((double)pars.nY - 1.0);
	}
	printf("x, y arrays allocated\n");
	
	pars.dx = pars.x[2] - pars.x[1];
	pars.dy = pars.y[2] - pars.y[1];
	
	/*Assign memory locations for wavenumber arrays*/
	pars.kX = 0;
	pars.kX = (double*)mkl_malloc(pars.nX*sizeof(double),64);
	if (pars.kX == 0) goto failed;
	pars.kY = 0;
	pars.kY = (double*)mkl_malloc(pars.nY*sizeof(double),64);
	if (pars.kY == 0) goto failed;
	
	/*Fill elements of kx, ky arrays*/
	double aX, aY;
	aX = -pars.nX / 2.0;
	aY = -pars.nY / 2.0;
	double bX, bY;
	bX = pars.nX / 2.0 - 1;
	bY = pars.nY / 2.0 - 1;
	double stepX, stepY;
	stepX = (2*M_PI / pars.xGridLength) * ((bX - aX) / (pars.nX - 1.0));
	stepY = (2*M_PI / pars.yGridLength) * ((bY - aY) / (pars.nY - 1.0));
	
	for (int i = 0; i < pars.nX; ++i)
	{
		pars.kX[i] = (2 * M_PI / pars.xGridLength) * aX + i*stepX;
	}
	for (int i = 0; i < pars.nY; ++i)
	{
		pars.kY[i] = (2 * M_PI / pars.yGridLength) * aY + i*stepY;
	}
	
	printf("kx, ky arrays allocated\n");

	/*Swap left and right hand sides of wavenumber arrays so they are inline with the DFT*/
	double tmp;
	int n2X, n2Y;
	n2X= pars.nX / 2.0;
	n2Y = pars.nY / 2.0;
	for (int i = 0; i < n2X; ++i)
	{
		tmp = pars.kX[i];
		pars.kX[i] = pars.kX[i + n2X];
		pars.kX[i+n2X] = tmp;
	}
	for (int i = 0; i < n2Y; ++i)
	{
		tmp = pars.kY[i];
		pars.kY[i] = pars.kY[i + n2Y];
		pars.kY[i+n2Y] = tmp;
	}
	printf("FFT shift of k arrays complete\n");
	
	/*Initialize psi as a 1D array accessed as psi[i*nRows + j]*/
	psi = (MKL_Complex16*)mkl_malloc(pars.N*sizeof(MKL_Complex16),64);
	init(psi,pars);
	if (psi == 0) goto failed;
	printf("Psi initialized successfully\n");
	
	/*Allocate memory for absPsi*/
	absPsi = (double*)mkl_malloc(pars.N*sizeof(double),64);
	if (absPsi==0) goto failed;
	printf("Kinetic energy initialized successfully\n");
	
	/*Construct Kinetic Energy Array*/
	MKL_Complex16 *kinEnergy = 0;
	kinEnergy = (MKL_Complex16*)mkl_malloc(pars.N*sizeof(MKL_Complex16),64);
	createKineticEnergy(kinEnergy, pars, true);
	MKL_Complex16 *kinEnergyY = 0;
	kinEnergyY = (MKL_Complex16*)mkl_malloc(pars.N*sizeof(MKL_Complex16),64);
	createKinEnergyY(kinEnergyY,pars);
	posPot = (MKL_Complex16*)mkl_malloc(pars.N*sizeof(MKL_Complex16),64);
	
	//////////////////////////////////// PREPARE THE DISCRETE FOURIER TRANSFORM ROUTINES /////////////////////////////////////////
	/*Check version of DFTI package*/
	char version[DFTI_VERSION_LENGTH];
	DftiGetValue(0, DFTI_VERSION, version);	
	printf("%s\n", version);
	/*Prints details about the DFTI descriptor*/
	printf(" DFTI_PRECISION      = DFTI_DOUBLE\n");
	printf(" DFTI_FORWARD_DOMAIN = DFTI_COMPLEX\n");
	printf(" DFTI_DIMENSION      = 2 one dimensional transforms\n");
	printf(" DFTI_LENGTHS        = {%i, %i}\n", pars.nX,pars.nY);
	/*Initialize the DFTI descriptor*/
	printf("Creating DFTI descriptor\n");
////////////////////////////////////////////////////////////////////////////////////
	/*Performs a Fourier transform along one axis at a time*/	
	MKL_LONG Nx = pars.nX, Ny = pars.nY;
	DftiCreateDescriptor(&descx, DFTI_DOUBLE, DFTI_COMPLEX, 1, Nx);
	DftiSetValue(descx, DFTI_NUMBER_OF_TRANSFORMS, Ny);
	DftiSetValue(descx, DFTI_INPUT_DISTANCE, Nx);
	DftiSetValue(descx, DFTI_BACKWARD_SCALE, 1.0/Nx);
	DftiCommitDescriptor(descx);
	
	MKL_LONG strides[] = {0, Nx};
        DftiCreateDescriptor(&descy, DFTI_DOUBLE, DFTI_COMPLEX, 1, Ny); 
        DftiSetValue(descy, DFTI_NUMBER_OF_TRANSFORMS, Nx);
        DftiSetValue(descy, DFTI_INPUT_DISTANCE, 1);
	DftiSetValue(descy, DFTI_BACKWARD_SCALE, 1.0/Ny);
	DftiSetValue(descy, DFTI_INPUT_STRIDES,  strides);
        DftiCommitDescriptor(descy);
/////////////////////////////////////////////////////////////////////////////////
	MKL_LONG N[2]; N[0] = pars.nX; N[1] = pars.nY;
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
	if (0 != status) goto failed;
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, (1.0/(N[0]*N[1])));
	if (0 != status) goto failed;
	printf("Commit descriptor\n");
	status = DftiCommitDescriptor(handle);
	if (status != 0) goto failed;
	printf("Setting DFTI scale for inverse transform\n");
	
	printf("Beginning imaginary time iteration\n");
	/*Loop which performs the TSSP algorithm*/
	int modFac = 500;
	
	double tempR = 0;
	double tempI = 0;
	double psiSum = 0;
	
	/*Normalize Psi*/
	/*Compute absolute value of psi*/
	vzAbs(pars.N,psi,absPsi);
	/*Square absPsi*/
	vdMul(pars.N,absPsi,absPsi,absPsi);
	
	for(int ii = 0; ii < pars.N; ++ii)
	{
		psiSum += absPsi[ii];		
	}

	for(int ii = 0; ii < pars.N; ++ii)
	{	
		tempR = psi[ii].real * sqrt(pars.nAtoms/(psiSum * pars.dx * pars.dy));
		tempI = psi[ii].imag * sqrt(pars.nAtoms/(psiSum * pars.dx * pars.dy));
		psi[ii].real = tempR;
		psi[ii].imag = tempI;
	}

	////////////////////////////////////////////////////
	psiSum = 0;
	remove("initPsi.fits");
	saveArray(absPsi, pars.N, "initPsi.fits", pars, 0);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*Calculate ground state of trap*/
	/*IMAGINARY LOOP BEGINS HERE*/
	if(pars.imProp){
		printf("IMAGINARY TIME PROGRESS:\n");
		for (int i = 0; i < pars.iSteps; ++i)
		{
		
			if(i % modFac == 0)
			{
				printf("Imaginary time step %d out of %d\n", i, pars.iSteps);
			}
					
			/*Forward FFT*/
			status = DftiComputeForward(handle, psi, psi);	
			if (0 != status) goto failed;
		
			/*Multiplies the momentum space wavefunction with kinEnergy and saves result to psi. Is an in-place routine.*/
			vzMul(pars.N,psi,kinEnergy,psi);

			/*Backward FFT*/		
			status = DftiComputeBackward(handle, psi, psi);
			if (0 != status) goto failed;
			
			/*Calculates the nonlinear potential exponential*/	
			createNonlinearEnergy(posPot, psi, pars, true, true);
		
			/*Multiplies psi with posPot and saves result to psi. Is an in-place routine.*/		
			vzMul(pars.N, psi, posPot, psi);
		
			/*Compute absolute value of psi*/
			vzAbs(pars.N,psi,absPsi);
			/*Square absPsi*/
			vdMul(pars.N,absPsi,absPsi,absPsi);
		
			/*Calculates the sum of |psi|^2*/
			for(int ii = 0; ii < pars.N; ++ii)
			{
				psiSum += absPsi[ii];	
			}
		
			/*Normalizes wavefunction so integral is = 1*/
			for(int ii = 0; ii < pars.N; ++ii)
			{	
				tempR = psi[ii].real * sqrt(pars.nAtoms/(psiSum * pars.dx * pars.dy));
				tempI = psi[ii].imag * sqrt(pars.nAtoms/(psiSum * pars.dx * pars.dy));
				psi[ii].real = tempR;
				psi[ii].imag = tempI;	
			}
			psiSum = 0;
		}	
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
		printf("\033[0m");
		printf("\033[300C\n");
		printf("Ground state computation complete\n");
		
		//Compute absolute value of psi
		vzAbs(pars.N,psi,absPsi);
		//Square absPsi
		vdMul(pars.N,absPsi,absPsi,absPsi);
		
		//Save to text file initialPsi.txt
		remove("groundState.fits");
		saveArray(absPsi, pars.N, "groundState.fits", pars, 0);
		//savefile(psi, "initPsi.txt", pars);
		vzAbs(pars.N,psi,absPsi);
		for(int i = 0; i < pars.N; ++i)
		{
			psi[i].real = absPsi[i];
			psi[i].imag = 0;
		}
	}
	else
	{
		//loadfile(psi, pars.N, "initPsi.txt");
		//printf("Ground state loaded\n");
		init(psi,pars);
	}
	
	/*Reassign kinetic energy and harmonic trap arrays for real time evolution*/

	//createKineticEnergy(kinEnergy, pars, false);
	//printf("Ek and V arrays reassigned for real time evolution.\n");
	
	
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////	
	
	/*REAL LOOP BEGINS HERE*/
	printf("REAL TIME PROGRESS:\n");
	int count = 0;
	for (int i = 0; i < pars.nSteps; ++i)
	{

		diagonalize(pMinusA, kinEnergy, maxCount, maxCount, pars);
		if (i % modFac == 0)
		{ 
			printf("Real time step %d out of %d\n", i, pars.nSteps);
			
			//Compute absolute value of psi
			vzAbs(pars.N,psi,absPsi);
			//Square absPsi
			vdMul(pars.N,absPsi,absPsi,absPsi);
			//Save to text file initialPsi.txt
			std::string currFileName = "fits/psi" + std::to_string(i/modFac) + ".fits";
			saveArray(absPsi, pars.N, currFileName.c_str(), pars, 0);
		}

//		remove("energyX.fits");
//		saveArray(pMinusA, pars.N, "energyX.fits", pars, 0);
		
		/*Forward FFT*/
		status = DftiComputeForward(descx, psi, psi);
		if (0 != status) goto failed;
		vzMul(pars.N,psi,kinEnergy,psi);
		status = DftiComputeBackward(descx, psi, psi);
		if (0 != status) goto failed;
		status = DftiComputeForward(descy, psi, psi);
		if (0 != status) goto failed;
		vzMul(pars.N, psi, kinEnergyY, psi);
		status = DftiComputeBackward(descy, psi, psi);
		if (0 != status) goto failed;

		createNonlinearEnergy(posPot, psi, pars, false, true);	
		vzMul(pars.N, psi, posPot, psi);
		
		status = DftiComputeForward(descy, psi, psi);
		if (0 != status) goto failed;
		vzMul(pars.N, psi, kinEnergyY, psi);		
		status = DftiComputeBackward(descy, psi, psi);
		if (0 != status) goto failed;
		status = DftiComputeForward(descx, psi, psi);
		if (0 != status) goto failed;
		vzMul(pars.N,psi,kinEnergy,psi);
		status = DftiComputeBackward(descx, psi, psi);
		
		/*Compute absolute value of psi*/
		vzAbs(pars.N,psi,absPsi);
		/*Square absPsi*/
		vdMul(pars.N,absPsi,absPsi,absPsi);
		
		/*Calculates the sum of |psi|^2*/
		for(int ii = 0; ii < pars.N; ++ii)
		{
			psiSum += absPsi[ii];	
		}
		
		/*Normalizes wavefunction so integral is = 1*/
		for(int ii = 0; ii < pars.N; ++ii)
		{	
			tempR = psi[ii].real * sqrt(1.0/(psiSum * pars.dx * pars.dy));
			tempI = psi[ii].imag * sqrt(1.0/(psiSum * pars.dx * pars.dy));
			psi[ii].real = tempR;
			psi[ii].imag = tempI;
		}
		psiSum = 0;
		count++;
	}
	printf("\033[0m");
	printf("\033[300C\n");
	printf("Real time computation complete.\n");
	goto cleanup;
	
	cleanup:
		
		printf("Free DFTI Descriptor\n");
		DftiFreeDescriptor(&handle);
		printf("Free x, y arrays\n");
		mkl_free(pars.x);
		mkl_free(pars.y);
		mkl_free(pars.kX);
		mkl_free(pars.kY);
		mkl_free(posPot);
		mkl_free(psi);
		mkl_free(absPsi);
		mkl_free(kinEnergy);
		printf("TSSP Routine %s\n", 0==status ? "PASSED" : "FAILED");
		return status;
	
	failed:
		printf(DftiErrorMessage(status));
		printf(" ERROR, status = "LI"\n", status);
		status = 1;
		goto cleanup;

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
	fits_update_key(fptr, TDOUBLE, "LX", &pars.xGridLength, "Length in X", &status);
	fits_update_key(fptr, TDOUBLE, "LY", &pars.yGridLength, "Length in Y", &status);
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

void savefile(MKL_Complex16 *psi, const char *filename, struct simPars pars)
{
	remove(filename);
	std::ofstream output(filename);
	output << pars.nX << " " << pars.nY << "\n\n";
	for (int i = 0; i < pars.N; ++i)
	{
		output << psi[i].real << " " << psi[i].imag << "\n";
    	}
    	output.close();

}

void loadfile(MKL_Complex16 *psi, int n, const char * filename)
{
	int nx;
	int ny;
	
	std::ifstream file(filename);
	
	file >> nx;
	file >> ny;
	
	for (int i = 0; i < nx*ny; i++) {
		
		file >> psi[i].real;
		file >> psi[i].imag;
	}
	file.close();

}
