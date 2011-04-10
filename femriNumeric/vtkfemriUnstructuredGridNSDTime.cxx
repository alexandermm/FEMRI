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
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkQuadraticTriangle.h"

//alexmbcm
#include <ctime>


vtkStandardNewMacro(vtkfemriUnstructuredGridNSDTime);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridNSDTime, "$Revision: 1.13 $");

vtkfemriUnstructuredGridNSDTime::vtkfemriUnstructuredGridNSDTime()
{
  this->MagnetizationValue = 1.0;
  this->SetNumberOfInputPorts(1);
  
  this->qLab = NULL;
  this->qHab = NULL;
  this->qGab = NULL;
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
}

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

    frequency[0] = 2.0;
	frequency[1] = 3.0;
	frequency[2] = 4.0;
	
	
	double ga[6];
	ga[0] = -twoPi*((2*xs[0]-4*xs[5]+2*xs[2])*frequency[0] + 
					(2*ys[0]-4*ys[5]+2*ys[2])*frequency[1] + 
					(2*zs[0]-4*zs[5]+2*zs[2])*frequency[2]);
					
	ga[1] = -twoPi*((2*xs[0]-4*xs[3]+2*xs[1])*frequency[0] + 
	                (2*ys[0]-4*ys[3]+2*ys[1])*frequency[1] + 
					(2*zs[0]-4*zs[3]+2*zs[1])*frequency[2]);
					
	ga[2] = -twoPi*((4*xs[0]-4*xs[3]+4*xs[4]-4*xs[5])*frequency[0] +
	                (4*ys[0]-4*ys[3]+4*ys[4]-4*ys[5])*frequency[1] +
					(4*zs[0]-4*zs[3]+4*zs[4]-4*zs[5])*frequency[2]);
					
	ga[3] = -twoPi*((4*xs[5]-3*xs[0]-xs[2])*frequency[0] +
	                (4*ys[5]-3*ys[0]-ys[2])*frequency[1] +
					(4*zs[5]-3*zs[0]-zs[2])*frequency[2]);
					
	ga[4] = -twoPi*((4*xs[3]-3*xs[0]-xs[1])*frequency[0] +
	                (4*ys[3]-3*ys[0]-ys[1])*frequency[1] +
					(4*zs[3]-3*zs[0]-zs[1])*frequency[2]);
					
	ga[5] = -twoPi*(xs[0]*frequency[0] + ys[0]*frequency[1] + zs[0]*frequency[2]);

	
	
	//Decompose integral
	bool nsdSuccess;
	if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
	{
		nsdSuccess = 0;
		cellValue[0] = 0.0;
		cellValue[1] = 0.0;
	}
	else
	{
		nsdSuccess = decompNSD(cell, cellValue, ga, a,b, numberOfCellPoints);
	}
	
			
	int numberOfQuadraturePoints = 4;
	
	//Get time
    double theTime = (double) (clock() - startTime) / CLOCKS_PER_SEC * 1000.0;
	
	//Record values for each element (cell)
	if (nsdSuccess)
		{ 
		writer << "0  " << frequency[0] << " " << frequency[1] << " " << frequency[2] << "  " << i << "  " << numberOfQuadraturePoints << "  " << theTime << "  " << cellValue[0] << " " << cellValue[1] << endl;
		}
	else
		{
		writer << "1  " << frequency[0] << " " << frequency[1] << " " << frequency[2] << "  " << i << "  " << numberOfQuadraturePoints << "  " << theTime << "  " << cellValue[0] << " " << cellValue[1] << endl;
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
double* ga, double a, double b, int numberOfCellPoints)
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
		G1 = innerDecompA(cell, dirY, ga, a,b);
        G4 = innerDecompB(cell, dirY, ga, a,b);
		
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

//GENERAL FUNCTIONS
inline comp fx(comp x, comp y)
{
	return 0.8*x*x + 0.5*y*y + 0.1*x*y + 0.7*x - 0.4*y + 0.1;
}
inline comp gx(comp x, comp y, double* ga)
{
	return ga[0]*x*x + ga[1]*y*y + ga[2]*x*y + ga[3]*x + ga[4]*y + ga[5];
}


//G FUNCTIONS
inline comp  h_ay(int dir, double p, double x, double* ga)
{
	return ( -(ga[2]*x+ga[4]) + ((double) dir)*sqrt(pow(ga[2]*x+ga[4], 2.0) + 4.0*ga[1]*comp(0.0,1.0)*p) ) /          
		   (2.0*ga[1]);
}
inline comp dh_ay(int dir, double p, double x, double *ga)
{
	return ((double) dir)*comp(0.0,1.0) / sqrt(pow(ga[2]*x+ga[4], 2.0) + 4.0*ga[1]*comp(0.0,1.0)*p);
}


//INNER DECOMPOSITIONS
comp vtkfemriUnstructuredGridNSDTime::innerDecompA(vtkCell* cell, int dirY, double* ga, double a, double b)
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

			f = fx(quadPoint[0], h_ay(dirY, quadPoint[1], quadPoint[0],ga)    );
			g = gx(quadPoint[0], h_ay(dirY, quadPoint[1], quadPoint[0],ga), ga);
			
			G = G + quadWeight*f*exp(comp(0.0,1.0)*g) * dh_ay(dirY, quadPoint[1], quadPoint[0], ga);
		}
		return G;
	}

	if (fabs(a1) < 100*eps) //If a1 is 0, oscillator is linear  
    {
		comp F1 = 0.0; 
		comp F2 = 0.0;
	
		for q = 1:size(qWL,1)
		{
			%CALCULATE F1
			u = hl(a, 1, qPL(q,1), a2);
			f = fx(u, h_ay(dirY, qPL(q,2), u,ga), fa);
			F1 = F1 + qWL(q)*f*dh_ay(dirY, qPL(q,2), u,ga)*dhl(1, a2);

			%CALCULATE F4
			u = hl(b, 1, qPL(q,1), a2);
			f = fx(u, h_ay(dirY, qPL(q,2), u,ga), fa);
			F2 = F2 + qWL(q)*f*dh_ay(dirY, qPL(q,2), u,ga)*dhl(1, a2);
		}

		G = (exp(i*polyval(gax,a ))*F1 - exp(i*polyval(gax,b ))*F2);
    }     
	else //If a1 is not 0, its quadratic
	{ 
		//Switch path direction based on direction of h on imaginary axis from real axis
		int dir;
		if (2.0*a*a1+a2 < 0.0) dir = -1; else dir = 1;
    
		if (sp < a) || (sp > b) %If sp is outside [a,b], don't use it 
        F1=0; F4=0;
        %decompAnoSP = 1
        for q = 1:size(qWL,1)
            %CALCULATE F1
            u = hab( dir, qPL(q,1), a,a1,a2);
            f = fx(u, h_ay(dirY, qPL(q,2), u,ga), fa);
            F1 = F1 + qWL(q)*f*dh_ay(dirY, qPL(q,2), u,ga)*dhab( dir, qPL(q,1), a,a1,a2);

            %CALCULATE F4
            u = hab( dir, qPL(q,1), b,a1,a2);
            f = fx(u, h_ay(dirY, qPL(q,2), u,ga), fa);
            F4 = F4 + qWL(q)*f*dh_ay(dirY, qPL(q,2), u,ga)*dhab( dir, qPL(q,1), b,a1,a2);
        end

        G = (exp(i*polyval(gax,a ))*F1 - ...
             exp(i*polyval(gax,b ))*F4);
         
		else %If sp is inside [a,b], calculate NSD integral with sp  
        F1=0; F2=0; F3=0; F4=0;
        %decompASP = 1
        for q = 1:size(qWL,1)
            %CALCULATE F1
            u = hab( dir, qPL(q,1), a,a1,a2);
            f = fx(u, h_ay(dirY, qPL(q,2), u,ga), fa);
            F1 = F1 + qWL(q)*f*dh_ay(dirY, qPL(q,2), u,ga)*dhab( dir, qPL(q,1), a,a1,a2);

            %CALCULATE F2 and F3
            u = hsp( dir, qPH(q,1), a1,a2);
            f = fx(u, h_ay(dirY, qPH(q,2), u,ga), fa);
            F2 = F2 + qWH(q)*f*dh_ay(dirY, qPH(q,2), u,ga)*dhsp( dir,a1);

            u = hsp(-dir, qPH(q,1), a1,a2);
            f = fx(u, h_ay(dirY, qPH(q,2), u,ga), fa);
            F3 = F3 + qWH(q)*f*dh_ay(dirY, qPH(q,2), u,ga)*dhsp(-dir,a1);

            %CALCULATE F4
            u = hab(-dir, qPL(q,1), b,a1,a2);
            f = fx(u, h_ay(dirY, qPL(q,2), u,ga), fa);
            F4 = F4 + qWL(q)*f*dh_ay(dirY, qPL(q,2), u,ga)*dhab(-dir, qPL(q,1), b,a1,a2);
        end

        G = (exp(i*polyval(gax,a ))*F1 - ...
             exp(i*polyval(gax,sp))*F2 + ...
             exp(i*polyval(gax,sp))*F3 - ...
             exp(i*polyval(gax,b ))*F4);
    end
end

}

comp vtkfemriUnstructuredGridNSDTime::innerDecompB(vtkCell* cell, int dirY, double* ga, double a, double b)
{
	return comp(1.0, 0.0);
}





