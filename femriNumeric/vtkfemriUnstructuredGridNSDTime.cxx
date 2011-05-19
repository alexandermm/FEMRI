/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridNSDTime.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/04 15:46:07 $
  Version:   $Revision: 1.13 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriUnstructuredGridNSDTime.h"
#include "vtkfemriNSDQuadrature.h"
#include "vtkfemriGaussArbQuadrature.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkQuadraticTriangle.h"

//alexmbcm
#include <ctime>  //For timing algorithm
#include <cfloat> //For definition of double type epsilon 
//#include <iostream>
//using namespace std;

vtkStandardNewMacro(vtkfemriUnstructuredGridNSDTime);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridNSDTime, "$Revision: 1.13 $");

vtkfemriUnstructuredGridNSDTime::vtkfemriUnstructuredGridNSDTime()
{
  this->MagnetizationValue = 1.0;
  this->SetNumberOfInputPorts(1);
  
  this->qLab = NULL;
  this->qHab = NULL;
  this->qGab = NULL;
  this->gaussQuadrature = NULL;
}

vtkfemriUnstructuredGridNSDTime::~vtkfemriUnstructuredGridNSDTime()
{
	if (this->qLab)
    {
		this->qLab->Delete(); this->qLab = NULL;
    }
	if (this->qHab)
    {
		this->qHab->Delete(); this->qHab = NULL;
    }
	if (this->qGab)
    {
		this->qGab->Delete(); this->qGab = NULL;
    }
	if (this->gaussQuadrature)
    {
		this->gaussQuadrature->Delete(); this->gaussQuadrature = NULL;
    }

}

int vtkfemriUnstructuredGridNSDTime::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

//This function is called from parent class in femriCommon/vtkfemriKSpaceGenerator.h
void vtkfemriUnstructuredGridNSDTime::Initialize()
{
	//Declare and initialize quadrature objects
	if (this->qLab)
    {
		this->qLab->Delete(); this->qLab = NULL;
    }
	this->qLab = vtkfemriNSDQuadrature::New();
	this->qLab->SetOrder(this->QuadratureOrder);
   	this->qLab->Initialize(GAUSS_LAGUERRE, GAUSS_LAGUERRE);
	
	if (this->qHab)
    {
		this->qHab->Delete(); this->qHab = NULL;
    }
	this->qHab = vtkfemriNSDQuadrature::New();
	this->qHab->SetOrder(this->QuadratureOrder);
   	this->qHab->Initialize(GAUSS_HALF_HERMITE, GAUSS_LAGUERRE);
	
	if (this->qGab)
    {
		this->qGab->Delete(); this->qGab = NULL;
    }
	this->qGab = vtkfemriNSDQuadrature::New();
	this->qGab->SetOrder(this->QuadratureOrder);
   	this->qGab->Initialize(GAUSS_LEGENDRE, GAUSS_LAGUERRE);
	
	//Declare and initialize normal quadrature object
	if (this->gaussQuadrature)
    {
		this->gaussQuadrature->Delete(); this->gaussQuadrature = NULL;
    }
	this->gaussQuadrature = vtkfemriGaussArbQuadrature::New();
	this->gaussQuadrature->SetOrder(this->QuadratureOrder);
    this->gaussQuadrature->Initialize(VTK_QUADRATIC_TRIANGLE);
}

/*inline double fx(double xi, double eta)
{
	return 8.03*xi*xi + 5.08*eta*eta + 1.42*xi*eta + 0.78*xi - 0.47*eta + 0.19;
}*/	
	
//Evaluate fourier transform at given k space values
void vtkfemriUnstructuredGridNSDTime::EvaluateFourierFunction(double frequency[3], 
double value[2], ofstream& writer)
{
  //Get information for all vtk cells (elements)
  vtkUnstructuredGrid* input = vtkUnstructuredGrid::SafeDownCast(this->GetInput());
  
 	  
  //Complex signal real and imaginary values at the given kspace value
  value[0] = 0.0;
  value[1] = 0.0;
  double cellValue[2];
 
  //Loop to get signal for each element 
  int numberOfCells = input->GetNumberOfCells();
  
  int i;
  //cout << endl << " Changing k-space value: " << endl;
  for (i=0; i<numberOfCells; i++)
    {
	
	//alexmbcm: start timer
	clock_t startTime = clock();
	
	
    vtkCell* cell = input->GetCell(i);
    
	cellValue[0] = 0.0;
	cellValue[1] = 0.0;
	
	//If element is not tri6, get out of loop
	if (cell->GetCellType() != VTK_QUADRATIC_TRIANGLE)
	{
		vtkErrorMacro("femri Error: unsupported cell type.");
		break;
	}
	
	    
    int numberOfCellPoints = cell->GetNumberOfPoints();
 
 
	//Get oscillator g
	int q, j, k;
	
	double xs[6];
	double ys[6];
	double zs[6];
	double p[3];
	double a = 0.0;
	double b = 1.0;
	
	for (q=0; q<numberOfCellPoints; q++)
	{
		cell->GetPoints()->GetPoint(q,p);
		
		xs[q] = p[0];
		ys[q] = p[1];
		zs[q] = p[2];
	}
	
			
	double twoPi = 2.0 * vtkMath::Pi();

	double ga[6];
	ga[0] = -twoPi*((2*xs[0]-4*xs[3]+2*xs[1])*frequency[0] + 
	                (2*ys[0]-4*ys[3]+2*ys[1])*frequency[1] + 
					(2*zs[0]-4*zs[3]+2*zs[1])*frequency[2]);
					
	ga[1] = -twoPi*((2*xs[0]-4*xs[5]+2*xs[2])*frequency[0] + 
					(2*ys[0]-4*ys[5]+2*ys[2])*frequency[1] + 
					(2*zs[0]-4*zs[5]+2*zs[2])*frequency[2]);				
					
	ga[2] = -twoPi*((4*xs[0]-4*xs[3]+4*xs[4]-4*xs[5])*frequency[0] +
	                (4*ys[0]-4*ys[3]+4*ys[4]-4*ys[5])*frequency[1] +
					(4*zs[0]-4*zs[3]+4*zs[4]-4*zs[5])*frequency[2]);
	
	ga[3] = -twoPi*((4*xs[3]-3*xs[0]-xs[1])*frequency[0] +
	                (4*ys[3]-3*ys[0]-ys[1])*frequency[1] +
					(4*zs[3]-3*zs[0]-zs[1])*frequency[2]);				
													
	ga[4] = -twoPi*((4*xs[5]-3*xs[0]-xs[2])*frequency[0] +
	                (4*ys[5]-3*ys[0]-ys[2])*frequency[1] +
					(4*zs[5]-3*zs[0]-zs[2])*frequency[2]);
					
	ga[5] = -twoPi*(xs[0]*frequency[0] + ys[0]*frequency[1] + zs[0]*frequency[2]);

    //cout << " ga values: " << endl << setprecision(16) << 
	//ga[0]<<" "<<ga[1]<<" "<<ga[2]<<" "<<ga[3]<<" "<<ga[4]<<" "<<ga[5]<<" "<< endl;
			
	//Decompose integral
	bool nsdSuccess = 0;
	if ((frequency[0] == 0.0) && (frequency[1] == 0.0) && (frequency[2] == 0.0))
	{
		nsdSuccess = 0;
	}
	else
	{
		nsdSuccess = decompNSD(cell, cellValue, ga, a,b, numberOfCellPoints, frequency);
	}
	
	//cout << " NSD signal: " << setprecision(16) << cellValue[0] << " " << cellValue[1] << endl;
	
	//IF NSD CANT DO IT USE SIMEDREA METHOD
	if (nsdSuccess == 0)
	{
		double twoPi = 2.0 * vtkMath::Pi();
		int subId = 0;
		double localQuadPoint[2], globalQuadPoint[3];
		double quadratureWeight;
		double weights[6];
		double derivs[12];
		int numberOfQuadraturePoints = 0;
	
		//Loop to go though all quadrature points
		numberOfQuadraturePoints = this->gaussQuadrature->GetNumberOfQuadraturePoints();
	
		int q, j, k;
		for (q=0; q<numberOfQuadraturePoints; q++)
		{
			//Get quadrature point and weight
			this->gaussQuadrature->GetQuadraturePoint(q,localQuadPoint);
			quadratureWeight = this->gaussQuadrature->GetQuadratureWeight(q);
	  
			cell->EvaluateLocation(subId,localQuadPoint,globalQuadPoint,weights);
			//double kdotx = vtkMath::Dot(globalQuadPoint,frequency);
			
			//cellValue[0] += quadratureWeight * fx(localQuadPoint[0], localQuadPoint[1]) * cos(twoPi * kdotx);
			//cellValue[1] -= quadratureWeight * fx(localQuadPoint[0], localQuadPoint[1]) * sin(twoPi * kdotx);

				
			//Calculate globalQuadPoint (the quadrature point in global coordinates)
			cell->EvaluateLocation(subId,localQuadPoint,globalQuadPoint,weights);
	  
			//Calculate normal at local quadrature point
			double Repsilon[3] = {0.0, 0.0, 0.0};
			double Reta[3]     = {0.0, 0.0, 0.0};
			double normal[3];
	  
			//Get shape function derivatives at local quadrature point 
			//First part of derivs is delN/delEpsilon, second part delN/delEta
			vtkQuadraticTriangle::SafeDownCast(cell)->InterpolationDerivs(localQuadPoint,derivs);
	  	
			//Compute Repsilon, Reta and the normal
			double x[3];
			for (j=0; j<numberOfCellPoints; j++)
			{
				cell->GetPoints()->GetPoint(j,x);
		  
				for (k=0; k<3; k++)
				{
					Repsilon[k] += x[k] * derivs[j];
					Reta[k]     += x[k] * derivs[numberOfCellPoints+j];
				}
			}
			vtkMath::Cross(Repsilon,Reta,normal);
			vtkMath::Normalize(normal);
	  	  
			//Calculate the coordinate transformation jacobian	
			double jacobian = this->ComputeJacobian(cell,localQuadPoint);
		
			double kdotx = vtkMath::Dot(globalQuadPoint,frequency);
			double kdotn = vtkMath::Dot(normal,frequency);
			double twoPik2 = 
			twoPi*(frequency[0]*frequency[0]+frequency[1]*frequency[1]+frequency[2]*frequency[2]);

			if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
			{
				cellValue[0] += this->MagnetizationValue * jacobian * quadratureWeight * 
				(normal[0]*globalQuadPoint[0]+normal[1]*globalQuadPoint[1]+normal[2]*globalQuadPoint[2])/3.0;
				cellValue[1] += 0.0;
			}
			else
			{
				cellValue[0] += 
				this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * sin(twoPi * kdotx);
				cellValue[1] += 
				this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * cos(twoPi * kdotx);
			}
		}
	}

	//Get time in ms
    double theTime = (double) (clock() - startTime) / CLOCKS_PER_SEC * 1000.0;
	
	int numberOfQuadPoints = this->qLab->GetNumberOfQuadraturePoints();
	
	//Record values for each element (cell)
	if (nsdSuccess)
		{ 
		writer << setprecision(16) << "1  " << frequency[0] << " " << frequency[1] << " " << frequency[2] << "  " << i << "  " << numberOfQuadPoints << "  " << theTime << "  " << cellValue[0] << " " << cellValue[1] << endl;
		}
	else
		{
		writer << "0  " << frequency[0] << " " << frequency[1] << " " << frequency[2] << "  " << i << "  " << numberOfQuadPoints << "  " << theTime << "  " << cellValue[0] << " " << cellValue[1] << endl;
		}
	
	
	//Add cell values to total kspace signal
	value[0] += cellValue[0];
	value[1] += cellValue[1];
	
	//Delete any dynamically allocated arrays 
	//delete[] weights;
	//delete[] derivs;
    }
}


//Function needed for EvaluateFourierFunction
double vtkfemriUnstructuredGridNSDTime::ComputeJacobian(vtkCell* cell, double pcoords[3])
{
  double jacobian = 0.0;
  
  int cellDimension = cell->GetCellDimension();
  
  if (cellDimension != 2)
  {
    return 0.0;
  }

  int numberOfCellPoints = cell->GetNumberOfPoints();
  
  double* derivs = new double[2*numberOfCellPoints];
  
  vtkQuadraticTriangle::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
  
  int i, j;
  double jacobianMatrixTr[2][3];
  for (i=0; i<3; i++)
  {
    jacobianMatrixTr[0][i] = jacobianMatrixTr[1][i] = 0.0;
  }

  double x[3];
  for (j=0; j<numberOfCellPoints; j++)
  {
    cell->GetPoints()->GetPoint(j,x);
    for (i=0; i<3; i++)
    {
      jacobianMatrixTr[0][i] += x[i] * derivs[j];
      jacobianMatrixTr[1][i] += x[i] * derivs[numberOfCellPoints+j];
    }
  }
  delete[] derivs;

  double jacobianMatrixSquared[2][2];
  jacobianMatrixSquared[0][0] = vtkMath::Dot(jacobianMatrixTr[0],jacobianMatrixTr[0]);
  jacobianMatrixSquared[0][1] = vtkMath::Dot(jacobianMatrixTr[0],jacobianMatrixTr[1]);
  jacobianMatrixSquared[1][0] = vtkMath::Dot(jacobianMatrixTr[1],jacobianMatrixTr[0]);
  jacobianMatrixSquared[1][1] = vtkMath::Dot(jacobianMatrixTr[1],jacobianMatrixTr[1]);

  double jacobianSquared = vtkMath::Determinant2x2(jacobianMatrixSquared[0],jacobianMatrixSquared[1]);

  if (jacobianSquared < 0.0)
  {
    jacobianSquared = fabs(jacobianSquared);
  }

  jacobian = sqrt(jacobianSquared);
 
  return jacobian;
}



//NSD FUNCTIONS
bool vtkfemriUnstructuredGridNSDTime::decompNSD(vtkCell* cell, double* cellValue, 
double* ga, double a, double b, int numberOfCellPoints, double* freq)
{
	int dirY;
	bool nsdSuccess = 1;
	
	comp fVal, G1,G4;
	
	//Choose direction in y dimension
	if (ga[4] < 0.0) dirY = -1; 
	else dirY = 1;
	
	//Get x value of stationary point in 'y' boundaries
	double spa = -ga[4]/ga[2];
    double spb =  (2*ga[1]+ga[4])/(2*ga[1]-ga[2]);
	
	//If sp is outside [a,b], don't use it
	if  ( ((spa < a) || (spa > b)) && ((spb < a) || (spb > b)) )
	{
		nsdSuccess = 1;
		G1 = innerDecompA(cell, dirY, ga, a,b, freq);
        G4 = innerDecompB(cell, dirY, ga, a,b, freq);
		
		//cout << " G1: " << G1 << endl;
		//cout << " G2: " << G4 << endl;
		
		fVal = G1 - G4;
	}
	else
	{
		nsdSuccess = 0;
		//Same as fVal = complex(0.0, 0.0);
		fVal = 0.0;
	}
	
	cellValue[0] = fVal.real();
	cellValue[1] = fVal.imag();
	
	return nsdSuccess;
}


//COMPLEX MATH FUNCTIONS BASED ON VTK FUNCTIONS
//xi  -> x
//eta -> y
inline comp complexDot(const comp x[3], const comp   y[3])
{
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}
inline comp complexDot(const comp x[3], const double y[3])
{
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline comp complexDeterminant2x2(const comp c1[2], const comp c2[2]) 
{
	return c1[0]*c2[1] - c2[0]*c1[1];
}
inline void complexCross(const comp x[3], const comp y[3], comp z[3])
{
	z[0] = x[1] * y[2] - x[2] * y[1]; 
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}
inline comp complexNormalize(comp x[3])
{
	comp den; 
	if ( ( den = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) ) != 0.0 )
	{
		for (int i=0; i < 3; i++)
		{
			x[i] /= den;
		}
	}
	return den;
}


//GENERAL FUNCTIONS
inline comp computeComplexJacobian(comp xi, comp eta, vtkCell* cell)
{
	comp jacobian = 0.0;
  
	int cellDimension = cell->GetCellDimension();
  
	if (cellDimension != 2)
	{
		return 0.0;
	}

	int numberOfCellPoints = cell->GetNumberOfPoints();
  
	comp derivs[12];
  
    //delN/delXi
	derivs[0] =  4.0*(xi+eta) - 3.0;
	derivs[1] =  4.0*xi - 1.0; 
	derivs[2] =  0.0;
	derivs[3] =  4.0*(1.0-eta-2.0*xi);
	derivs[4] =  4.0*eta;
	derivs[5] = -4.0*eta;
	//delN/delEta
	derivs[6]  =  4.0*(xi+eta) - 3.0;
    derivs[7]  =  0.0; 
    derivs[8]  =  4.0*eta - 1.0;
    derivs[9]  = -4.0*xi;
    derivs[10] =  4.0*xi;
    derivs[11] =  4.0*(1.0-xi-2.0*eta);
  
	int i, j;
	comp jacobianMatrixTr[2][3];
	for (i=0; i<3; i++)
	{
		jacobianMatrixTr[0][i] = jacobianMatrixTr[1][i] = 0.0;
	}

	double x[3];
	for (j=0; j<numberOfCellPoints; j++)
	{
		cell->GetPoints()->GetPoint(j,x);
		for (i=0; i<3; i++)
		{
			jacobianMatrixTr[0][i] += x[i] * derivs[j];
			jacobianMatrixTr[1][i] += x[i] * derivs[numberOfCellPoints+j];
		}
	}

	comp jacobianMatrixSquared[2][2];
	jacobianMatrixSquared[0][0] = complexDot(jacobianMatrixTr[0],jacobianMatrixTr[0]);
	jacobianMatrixSquared[0][1] = complexDot(jacobianMatrixTr[0],jacobianMatrixTr[1]);
	jacobianMatrixSquared[1][0] = complexDot(jacobianMatrixTr[1],jacobianMatrixTr[0]);
	jacobianMatrixSquared[1][1] = complexDot(jacobianMatrixTr[1],jacobianMatrixTr[1]);

	comp jacobianSquared = complexDeterminant2x2(jacobianMatrixSquared[0],jacobianMatrixSquared[1]);

	return sqrt(jacobianSquared);
}


inline comp fx(comp xi, comp eta, vtkCell* cell, double* frequency)
{
	//return 8.03*xi*xi + 5.08*eta*eta + 1.42*xi*eta + 0.78*xi - 0.47*eta + 0.19;
	
	//Compute k dot normal
	int numberOfCellPoints = cell->GetNumberOfPoints();
	
	comp derivs[12];
  
    //delN/delXi
	derivs[0] =  4.0*(xi+eta) - 3.0;
	derivs[1] =  4.0*xi - 1.0; 
	derivs[2] =  0.0;
	derivs[3] =  4.0*(1.0-eta-2.0*xi);
	derivs[4] =  4.0*eta;
	derivs[5] = -4.0*eta;
	//delN/delEta
	derivs[6]  =  4.0*(xi+eta) - 3.0;
    derivs[7]  =  0.0; 
    derivs[8]  =  4.0*eta - 1.0;
    derivs[9]  = -4.0*xi;
    derivs[10] =  4.0*xi;
    derivs[11] =  4.0*(1.0-xi-2.0*eta);
	  
	//Calculate normal at local quadrature point
	comp Repsilon[3] = {0.0, 0.0, 0.0};
	comp Reta[3]     = {0.0, 0.0, 0.0};
	comp normal[3];
	  	  	
	//Compute Repsilon, Reta and the normal
	int j, k;
	double x[3];
	for (j=0; j<numberOfCellPoints; j++)
	{
		cell->GetPoints()->GetPoint(j,x);
		  
		for (k=0; k<3; k++)
		{
			Repsilon[k] += x[k] * derivs[j];
			Reta[k]     += x[k] * derivs[numberOfCellPoints+j];
		}
	}
	complexCross(Repsilon,Reta,normal);
	complexNormalize(normal);
	
	comp   kdotn   = complexDot(normal, frequency);
	double twoPik2 = 2.0*vtkMath::Pi() *
	(frequency[0]*frequency[0]+frequency[1]*frequency[1]+frequency[2]*frequency[2]);


	return computeComplexJacobian(xi, eta, cell) * kdotn * comp(0.0,1.0) / twoPik2;
}

inline comp gx(comp x, comp y, double* ga)
{
	return ga[0]*x*x + ga[1]*y*y + ga[2]*x*y + ga[3]*x + ga[4]*y + ga[5];
}
inline double polyval(double* gax, double x)
{
	return gax[0]*x*x + gax[1]*x + gax[2];
}


//G FUNCTIONS
inline comp  h_ay(int dir, double p, comp x, double* ga)
{
	return ( -(ga[2]*x+ga[4]) + ((double) dir)*sqrt(pow(ga[2]*x+ga[4], 2.0) + 4.0*ga[1]*comp(0.0,1.0)*p) ) /          
		   (2.0*ga[1]);
}
inline comp dh_ay(int dir, double p, comp x, double* ga)
{
	return ((double) dir)*comp(0.0,1.0) / sqrt(pow(ga[2]*x+ga[4],2.0) + 4.0*ga[1]*comp(0.0,1.0)*p);
}
inline comp  h_by(int dir, double p, comp x, double* ga)
{
	return ( -(ga[2]*x+ga[4]) + ((double) dir)*sqrt(pow((ga[2]-2.0*ga[1])*x + 2.0*ga[1] + ga[4],2.0) + 4.0*ga[1]*comp(0.0,1.0)*p) ) / (2.0*ga[1]);
}
inline comp dh_by(int dir, double p, comp x, double* ga)
{
	return ((double) dir)*comp(0.0,1.0) / sqrt(pow((ga[2]-2.0*ga[1])*x+2.0*ga[1]+ga[4], 2.0) + 4.0*ga[1]*comp(0.0,1.0)*p);
}


//F FUNCTIONS
inline comp  hl(double x, int dir, double p, double a2)
{
	return x + ((double) dir)*comp(0.0,1.0)*p/a2;
}
inline comp dhl(int dir, double a2)
{
	return ((double) dir)*comp(0.0,1.0)/a2;
}

inline comp  hab(int dir, double p, double x, double a1, double a2)
{
	return (-a2 + ((double) dir)*sqrt(pow(2.0*a1*x+a2,2.0) + 4.0*a1*comp(0.0,1.0)*p))/(2.0*a1);
}
inline comp dhab(int dir, double p, double x, double a1, double a2)
{
	return ((double) dir)*comp(0.0,1.0)/sqrt(pow(2.0*a1*x+a2,2.0) + 4.0*a1*comp(0.0,1.0)*p);
}

inline comp  hsp(int dir, double q, double a1, double a2)
{
	return (-a2 + ((double) dir)*q*sqrt(4.0*a1*comp(0.0,1.0)) ) / (2.0*a1);
}
inline comp dhsp(int dir, double a1)
{
	return ((double) dir)*comp(0.0,1.0) / sqrt(a1*comp(0.0,1.0));
}


//INNER DECOMPOSITIONS
comp vtkfemriUnstructuredGridNSDTime::innerDecompA(vtkCell* cell, int dirY, double* ga, double a, double b, double *freq)
{
	double gax[3];
	gax[0] = ga[0];
	gax[1] = ga[3];
	gax[2] = ga[5]; 
	
	double a1 = gax[0];
	double a2 = gax[1];
	double sp = -a2/(2.0*a1); //Stationary point
	int q;
	
	
    //If oscillator can be considered nonoscillatory in interval integrate using standard quadrature
	if ((fabs(a1) < 0.1) && (fabs(a2) < 0.1))
	{
		comp G = 0.0;
		comp f;
		comp g;
		double quadPoint[2];
		double quadWeight;
		
        int numQuadPs = this->qGab->GetNumberOfQuadraturePoints();

		for (q=0; q < numQuadPs; q++)
		{
			this->qGab->GetQuadraturePoint(q,quadPoint);
			quadWeight = this->qGab->GetQuadratureWeight(q);

			f = fx(quadPoint[0], h_ay(dirY, quadPoint[1], quadPoint[0],ga), cell, freq);
			g = gx(quadPoint[0], h_ay(dirY, quadPoint[1], quadPoint[0],ga), ga  );
			
			G = G + quadWeight*f*exp(comp(0.0,1.0)*g) * dh_ay(dirY, quadPoint[1], quadPoint[0], ga);
		}
		return G;
	}

	if (fabs(a1) < 100.0*DBL_EPSILON) //If a1 is 0, oscillator is linear  
    {
		comp F1 = 0.0; 
		comp F2 = 0.0;
		comp u, f, g;
		double quadPoint[2];
		double quadWeight;
		
        int numQuadPs = this->qLab->GetNumberOfQuadraturePoints();
		
		for (q=0; q < numQuadPs; q++)
		{
			this->qLab->GetQuadraturePoint(q,quadPoint);
			quadWeight = this->qLab->GetQuadratureWeight(q);
			
			//CALCULATE F1
			u = hl(a, 1, quadPoint[0], a2);
			f = fx(u, h_ay(dirY, quadPoint[1], u,ga), cell, freq);
			F1 = F1 + quadWeight*f * dh_ay(dirY, quadPoint[1], u,ga)*dhl(1, a2);

			//CALCULATE F2
			u = hl(b, 1, quadPoint[0], a2);
			f = fx(u, h_ay(dirY, quadPoint[1], u,ga), cell, freq);
			F2 = F2 + quadWeight*f * dh_ay(dirY, quadPoint[1], u,ga)*dhl(1, a2);
		}

		return exp(comp(0.0,1.0)*polyval(gax,a))*F1 - exp(comp(0.0,1.0)*polyval(gax,b))*F2;
    }     
	else //If a1 is not 0, its quadratic
	{ 
		//Switch path direction based on direction of h on imaginary axis from real axis
		int dir;
		if (2.0*a*a1+a2 < 0.0) dir = -1; else dir = 1;
    
		if ((sp < a) || (sp > b)) //If sp is outside [a,b], don't use it 
        {
			comp F1 = 0.0; 
			comp F2 = 0.0;
			comp u, f, g;
			double quadPoint[2];
			double quadWeight;
		
			int numQuadPs = this->qLab->GetNumberOfQuadraturePoints();
			
			for (q=0; q < numQuadPs; q++)
			{
				this->qLab->GetQuadraturePoint(q,quadPoint);
				quadWeight = this->qLab->GetQuadratureWeight(q);
				
				//CALCULATE F1
				u = hab( dir, quadPoint[0], a,a1,a2);
				f = fx(u, h_ay(dirY, quadPoint[1], u,ga), cell, freq);
				F1 = F1 + quadWeight*f * dh_ay(dirY, quadPoint[1], u,ga)*dhab( dir, quadPoint[0], a,a1,a2);

				//CALCULATE F2
				u = hab( dir, quadPoint[0], b,a1,a2);
				f = fx(u, h_ay(dirY, quadPoint[1], u,ga), cell, freq);
				F2 = F2 + quadWeight*f * dh_ay(dirY, quadPoint[1], u,ga)*dhab( dir, quadPoint[0], b,a1,a2);
			}

			return exp(comp(0.0,1.0)*polyval(gax,a))*F1 - exp(comp(0.0,1.0)*polyval(gax,b))*F2;
        } 
		else //If sp is inside [a,b], calculate NSD integral with sp
		{  
			comp F1 = 0.0; comp F2 = 0.0;
			comp F3 = 0.0; comp F4 = 0.0;
			comp u, f, g;
			double quadPointL[2];
			double quadWeightL;
			double quadPointH[2];
			double quadWeightH;
		
			int numQuadPs = this->qLab->GetNumberOfQuadraturePoints();

			for (q=0; q < numQuadPs; q++)
			{
				this->qLab->GetQuadraturePoint(q,quadPointL);
				quadWeightL = this->qLab->GetQuadratureWeight(q);
				this->qHab->GetQuadraturePoint(q,quadPointH);
				quadWeightH = this->qHab->GetQuadratureWeight(q);
				
				//CALCULATE F1
				u = hab( dir, quadPointL[0], a,a1,a2);
				f = fx(u, h_ay(dirY, quadPointL[1], u,ga), cell, freq);
				F1 = F1 + quadWeightL*f * dh_ay(dirY, quadPointL[1], u,ga)*dhab( dir, quadPointL[0], a,a1,a2);

				//CALCULATE F2 and F3
				u = hsp( dir, quadPointH[0], a1,a2);
				f = fx(u, h_ay(dirY, quadPointH[1], u,ga), cell, freq);
				F2 = F2 + quadWeightH*f * dh_ay(dirY, quadPointH[1], u,ga)*dhsp( dir,a1);

				u = hsp(-dir, quadPointH[0], a1,a2);
				f = fx(u, h_ay(dirY, quadPointH[1], u,ga), cell, freq);
				F3 = F3 + quadWeightH*f * dh_ay(dirY, quadPointH[1], u,ga)*dhsp(-dir,a1);

				//CALCULATE F4
				u = hab(-dir, quadPointL[0], b,a1,a2);
				f = fx(u, h_ay(dirY, quadPointL[1], u,ga), cell, freq);
				F4 = F4 + quadWeightL*f * dh_ay(dirY, quadPointL[1], u,ga)*dhab(-dir, quadPointL[0], b,a1,a2);
			}

        return exp(comp(0.0,1.0)*polyval(gax,a ))*F1 - exp(comp(0.0,1.0)*polyval(gax,sp))*F2 + 
			   exp(comp(0.0,1.0)*polyval(gax,sp))*F3 - exp(comp(0.0,1.0)*polyval(gax,b ))*F4;
		}
	}
}

comp vtkfemriUnstructuredGridNSDTime::innerDecompB(vtkCell* cell, int dirY, double* ga, double a, double b, double* freq)
{
	double gax[3];
	gax[0] = ga[0] + ga[1] - ga[2];
	gax[1] = ga[2] + ga[3] - ga[4] - 2.0*ga[1];
	gax[2] = ga[1] + ga[4] + ga[5]; 
	
	double a1 = gax[0];
	double a2 = gax[1];
	double sp = -a2/(2.0*a1); //Stationary point

	int q;
	
	
    //If oscillator can be considered nonoscillatory in interval integrate using standard quadrature
	if ((fabs(a1) < 0.1) && (fabs(a2) < 0.1))
	{
		comp G = 0.0;
		comp f;
		comp g;
		double quadPoint[2];
		double quadWeight;
		
        int numQuadPs = this->qGab->GetNumberOfQuadraturePoints();
		
		for (q=0; q < numQuadPs; q++)
		{
			this->qGab->GetQuadraturePoint(q,quadPoint);
			quadWeight = this->qGab->GetQuadratureWeight(q);

			f = fx(quadPoint[0], h_by(dirY, quadPoint[1], quadPoint[0],ga), cell, freq);
			g = gx(quadPoint[0], h_by(dirY, quadPoint[1], quadPoint[0],ga), ga  );
			
			G = G + quadWeight*f*exp(comp(0.0,1.0)*g) * dh_by(dirY, quadPoint[1], quadPoint[0], ga);
		}
		return G;
	}

	if (fabs(a1) < 100.0*DBL_EPSILON) //If a1 is 0, oscillator is linear  
    {
		comp F1 = 0.0; 
		comp F2 = 0.0;
		comp u, f, g;
		double quadPoint[2];
		double quadWeight;
		
        int numQuadPs = this->qLab->GetNumberOfQuadraturePoints();
		
		for (q=0; q < numQuadPs; q++)
		{
			this->qLab->GetQuadraturePoint(q,quadPoint);
			quadWeight = this->qLab->GetQuadratureWeight(q);
			
			//CALCULATE F1
			u = hl(a, 1, quadPoint[0], a2);
			f = fx(u, h_by(dirY, quadPoint[1], u,ga), cell, freq);
			F1 = F1 + quadWeight*f * dh_by(dirY, quadPoint[1], u,ga)*dhl(1, a2);

			//CALCULATE F2
			u = hl(b, 1, quadPoint[0], a2);
			f = fx(u, h_by(dirY, quadPoint[1], u,ga), cell, freq);
			F2 = F2 + quadWeight*f * dh_by(dirY, quadPoint[1], u,ga)*dhl(1, a2);
		}

		return exp(comp(0.0,1.0)*polyval(gax,a))*F1 - exp(comp(0.0,1.0)*polyval(gax,b))*F2;
    }     
	else //If a1 is not 0, its quadratic
	{ 
		//Switch path direction based on direction of h on imaginary axis from real axis
		int dir;
		if (2.0*a*a1+a2 < 0.0) dir = -1; else dir = 1;
    
		if ((sp < a) || (sp > b)) //If sp is outside [a,b], don't use it 
        {
			comp F1 = 0.0; 
			comp F2 = 0.0;
			comp u, f, g;
			double quadPoint[2];
			double quadWeight;
		
			int numQuadPs = this->qLab->GetNumberOfQuadraturePoints();
			
			for (q=0; q < numQuadPs; q++)
			{
				this->qLab->GetQuadraturePoint(q,quadPoint);
				quadWeight = this->qLab->GetQuadratureWeight(q);
				
				//CALCULATE F1
				u = hab( dir, quadPoint[0], a,a1,a2);
				f = fx(u, h_by(dirY, quadPoint[1], u,ga), cell, freq);
				F1 = F1 + quadWeight*f * dh_by(dirY, quadPoint[1], u,ga)*dhab( dir, quadPoint[0], a,a1,a2);

				//CALCULATE F2
				u = hab( dir, quadPoint[0], b,a1,a2);
				f = fx(u, h_by(dirY, quadPoint[1], u,ga), cell, freq);
				F2 = F2 + quadWeight*f * dh_by(dirY, quadPoint[1], u,ga)*dhab( dir, quadPoint[0], b,a1,a2);
			}
			
			return exp(comp(0.0,1.0)*polyval(gax,a))*F1 - exp(comp(0.0,1.0)*polyval(gax,b))*F2;
        } 
		else //If sp is inside [a,b], calculate NSD integral with sp
		{  
			comp F1 = 0.0; comp F2 = 0.0;
			comp F3 = 0.0; comp F4 = 0.0;
			comp u, f, g;
			double quadPointL[2];
			double quadWeightL;
			double quadPointH[2];
			double quadWeightH;
		
			int numQuadPs = this->qLab->GetNumberOfQuadraturePoints();

			for (q=0; q < numQuadPs; q++)
			{
				this->qLab->GetQuadraturePoint(q,quadPointL);
				quadWeightL = this->qLab->GetQuadratureWeight(q);
				this->qHab->GetQuadraturePoint(q,quadPointH);
				quadWeightH = this->qHab->GetQuadratureWeight(q);
				
				//CALCULATE F1
				u = hab( dir, quadPointL[0], a,a1,a2);
				f = fx(u, h_by(dirY, quadPointL[1], u,ga), cell, freq);
				F1 = F1 + quadWeightL*f * dh_by(dirY, quadPointL[1], u,ga)*dhab( dir, quadPointL[0], a,a1,a2);

				//CALCULATE F2 and F3
				u = hsp( dir, quadPointH[0], a1,a2);
				f = fx(u, h_by(dirY, quadPointH[1], u,ga), cell, freq);
				F2 = F2 + quadWeightH*f * dh_by(dirY, quadPointH[1], u,ga)*dhsp( dir,a1);

				u = hsp(-dir, quadPointH[0], a1,a2);
				f = fx(u, h_by(dirY, quadPointH[1], u,ga), cell, freq);
				F3 = F3 + quadWeightH*f * dh_by(dirY, quadPointH[1], u,ga)*dhsp(-dir,a1);

				//CALCULATE F4
				u = hab(-dir, quadPointL[0], b,a1,a2);
				f = fx(u, h_by(dirY, quadPointL[1], u,ga), cell, freq);
				F4 = F4 + quadWeightL*f * dh_by(dirY, quadPointL[1], u,ga)*dhab(-dir, quadPointL[0], b,a1,a2);
			}

        return exp(comp(0.0,1.0)*polyval(gax,a ))*F1 - exp(comp(0.0,1.0)*polyval(gax,sp))*F2 + 
			   exp(comp(0.0,1.0)*polyval(gax,sp))*F3 - exp(comp(0.0,1.0)*polyval(gax,b ))*F4;
		}
	}
}




