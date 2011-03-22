/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator.cxx,v $
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

#include "vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator.h"
#include "vtkfemriGaussArbQuadrature.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkQuadraticTriangle.h"

vtkStandardNewMacro(vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator, "$Revision: 1.13 $");

vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator::vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator()
{
  this->MagnetizationValue = 1.0;
  this->SetNumberOfInputPorts(1);
  
  this->gaussQuadrature = NULL;
}

vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator::~vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator()
{
	if (this->gaussQuadrature)
    {
		this->gaussQuadrature->Delete();
		this->gaussQuadrature = NULL;
    }
}

int vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

void vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator::Initialize()
{
	//Declare and initialize quadrature object
	//ASSUMING ALL CELLS ARE OF SAME TYPE
	if (this->gaussQuadrature)
    {
		this->gaussQuadrature->Delete();
		this->gaussQuadrature = NULL;
    }
	
	this->gaussQuadrature = vtkfemriGaussArbQuadrature::New();
	
	//Get same specified quadrature order for all elements
	this->gaussQuadrature->SetOrder(this->QuadratureOrder);
	
    //Generate 2D or 3D quadrature for square and triangular domains
	//with different Gauss type rules depending on element type 
    this->gaussQuadrature->Initialize(VTK_QUADRATIC_TRIANGLE);
}

//Evaluate fourier transform at given k space values
void vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator::EvaluateFourierFunction(double frequency[3], 
double value[2])
{
  //Get information for all vtk cells (elements)
  vtkUnstructuredGrid* input = vtkUnstructuredGrid::SafeDownCast(this->GetInput());
  
  //Complex signal real and imaginary values at the given kspace value
  value[0] = 0.0;
  value[1] = 0.0;

  //Loop to get signal for each element 
  int numberOfCells = input->GetNumberOfCells();
  
  int i;
  for (i=0; i<numberOfCells; i++)
    {
    vtkCell* cell = input->GetCell(i);
    
	//If dimension of element is not 2, can't apply sufrace femri algorithm 
    if (cell->GetCellDimension() != 2)
      {
      continue;
      }
	
	    
    int numberOfCellPoints = cell->GetNumberOfPoints();
 
    double twoPi = 2.0 * vtkMath::Pi();
    int subId = 0;
    double localQuadPoint[2], globalQuadPoint[3];
    double quadratureWeight;
    double* weights = new double[numberOfCellPoints];
	double* derivs = new double[2*numberOfCellPoints];
    int numberOfQuadraturePoints = 0;
    bool preComputedQuadratureRule = false;
	
	//Loop to go though all quadrature points
	numberOfQuadraturePoints = this->gaussQuadrature->GetNumberOfQuadraturePoints();
	
	int q, j, k;
    for (q=0; q<numberOfQuadraturePoints; q++)
      {
	  
	  //Get quadrature point and weight
      if (!preComputedQuadratureRule)
        {
        this->gaussQuadrature->GetQuadraturePoint(q,localQuadPoint);
        quadratureWeight = this->gaussQuadrature->GetQuadratureWeight(q);
        }
	  
	  //Calculate globalQuadPoint (the quadrature point in global coordinates)
      cell->EvaluateLocation(subId,localQuadPoint,globalQuadPoint,weights);
	  
	  //Calculate normal at local quadrature point
	  double Repsilon[3] = {0.0, 0.0, 0.0};
	  double Reta[3]     = {0.0, 0.0, 0.0};
	  double normal[3];
	  
	  //Get shape function derivatives at local quadrature point 
	  //First part of derivs is delN/delEpsilon, second part delN/delEta
	  //Added error catching when using vtkQuadraticTriangle for the first time
	  if (cell->GetCellType() == VTK_QUADRATIC_TRIANGLE)
		vtkQuadraticTriangle::SafeDownCast(cell)->InterpolationDerivs(localQuadPoint,derivs);
	  else 
		{
			vtkErrorMacro("femri Error: unsupported cell type.");
			return;
		}
	
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
      double twoPik2 = twoPi*(frequency[0]*frequency[0]+frequency[1]*frequency[1]+frequency[2]*frequency[2]);

      if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
        {
        value[0] += this->MagnetizationValue * jacobian * quadratureWeight * 
		(normal[0]*globalQuadPoint[0] + normal[1]*globalQuadPoint[1] + normal[2]*globalQuadPoint[2])/3.0;
        value[1] += 0.0;
        }
      else
        {
        value[0] += 
		this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * sin(twoPi * kdotx);
        value[1] += 
		this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * cos(twoPi * kdotx);
        }
      }
	//Delete any dynamically allocated arrays 
	delete[] weights;
	delete[] derivs;
    }

  this->gaussQuadrature->Delete();
}


//Function needed for EvaluateFourierFunction
double vtkfemriUnstructuredGridArbQuadTriKSpaceGenerator::ComputeJacobian(vtkCell* cell, double pcoords[3])
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

