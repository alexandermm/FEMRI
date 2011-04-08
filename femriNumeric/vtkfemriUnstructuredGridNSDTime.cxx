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
  
  this->gaussQuadrature = NULL;
}

vtkfemriUnstructuredGridNSDTime::~vtkfemriUnstructuredGridNSDTime()
{
	if (this->gaussQuadrature)
    {
		this->gaussQuadrature->Delete();
		this->gaussQuadrature = NULL;
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
	
	double xs[] = {0.0,0.0,0.0, 0.0,0.0,0.0};
	double ys[] = {0.0,0.0,0.0, 0.0,0.0,0.0};
	double zs[] = {0.0,0.0,0.0, 0.0,0.0,0.0};
	double p[3];
	
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

	
	if (i == 9)
	{
		cout << endl << ga[0] << " " << ga[1] << " " << ga[2] 
			  << " " << ga[3] << " " << ga[4] << " " << ga[5] << endl;
	}
	
	
	cellValue[0] = 0.0;
	cellValue[1] = 0.0;
	
	int numberOfQuadraturePoints = 4;
	
	//Get time
    double theTime = (double) (clock() - startTime) / CLOCKS_PER_SEC * 1000.0;
	
	//Record values for each element (cell)
	if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
		{ 
		writer << "0 0  " << frequency[0] << " " << frequency[1] << " " << frequency[2] << "  " << i << "  " << numberOfQuadraturePoints << "  " << theTime << "  " << cellValue[0] << " " << cellValue[1] << endl;
		}
	else
		{
		writer << "0 1  " << frequency[0] << " " << frequency[1] << " " << frequency[2] << "  " << i << "  " << numberOfQuadraturePoints << "  " << theTime << "  " << cellValue[0] << " " << cellValue[1] << endl;
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

