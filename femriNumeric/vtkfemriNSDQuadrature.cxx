/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriNSDQuadrature.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/03 17:00:30 $
  Version:   $Revision: 1.3 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriNSDQuadrature.h"
#include "vtkObjectFactory.h"
#include "vtkCellType.h"
//alexmbcm: change: added the following line
#include "vtkMath.h"

vtkStandardNewMacro(vtkfemriNSDQuadrature);
vtkCxxRevisionMacro(vtkfemriNSDQuadrature, "$Revision: 1.3 $");

vtkfemriNSDQuadrature::vtkfemriNSDQuadrature()
{
  this->QuadraturePoints = NULL;
  this->QuadratureWeights = NULL;

  this->Order = 1;
  this->PreviousOrder = 0;
  
  this->CellType = VTK_EMPTY_CELL;
}

vtkfemriNSDQuadrature::~vtkfemriNSDQuadrature()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }

  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
}

void vtkfemriNSDQuadrature::Initialize(vtkIdType cellType)
{
  if ((cellType == this->CellType) && (this->Order == this->PreviousOrder))
  {
    return;
  }

  this->CellType = cellType;
  this->PreviousOrder = this->Order;

  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();
 
  switch(cellType)
  {
    case VTK_LINE:
    case VTK_QUADRATIC_EDGE:
    {
      this->QuadraturePoints->SetNumberOfComponents(1);
      this->Initialize1DGauss();
      break;
    }
    case VTK_QUAD:
    case VTK_QUADRATIC_QUAD:
    {
      vtkfemriNSDQuadrature* q1D = vtkfemriNSDQuadrature::New();
      q1D->SetOrder(this->Order);
      q1D->Initialize1DGauss();
      this->TensorProductQuad(q1D);
      q1D->Delete(); 
      break;
    }
    case VTK_TRIANGLE:
    case VTK_QUADRATIC_TRIANGLE:
    {
      if (this->Order == 0 || this->Order ==1)
      {
        this->QuadraturePoints->SetNumberOfComponents(2);
        this->QuadraturePoints->SetNumberOfTuples(1);
        this->QuadratureWeights->SetNumberOfTuples(1);
        double point[2];
        double weight;
        point[0] = 0.33333333333333333333333333333333;
        point[1] = 0.33333333333333333333333333333333;
        weight = 0.5;
        this->QuadraturePoints->SetTuple(0,point);
        this->QuadratureWeights->SetValue(0,weight);
        break;
      }
      vtkfemriNSDQuadrature* gauss1D = vtkfemriNSDQuadrature::New();
      gauss1D->SetOrder(this->Order);
      gauss1D->Initialize1DGauss();
      vtkfemriNSDQuadrature* jacA1D = vtkfemriNSDQuadrature::New();
      jacA1D->SetOrder(this->Order);
      jacA1D->Initialize1DJacobi(1,0);
      this->TensorProductTriangle(gauss1D,jacA1D);
      gauss1D->Delete();
      jacA1D->Delete();
      break;
    }
    case VTK_HEXAHEDRON:
    case VTK_QUADRATIC_HEXAHEDRON:
    {
      vtkfemriNSDQuadrature* q1D = vtkfemriNSDQuadrature::New();
      q1D->SetOrder(this->Order);
      q1D->Initialize1DGauss();
      this->TensorProductHexahedron(q1D);
      q1D->Delete(); 
      break;
    }
    case VTK_WEDGE:
    case VTK_QUADRATIC_WEDGE:
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
    case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
#endif
    {
      vtkfemriNSDQuadrature* q1D = vtkfemriNSDQuadrature::New();
      q1D->SetOrder(this->Order);
      q1D->Initialize1DGauss();
      vtkfemriNSDQuadrature* q2D = vtkfemriNSDQuadrature::New();
      q2D->SetOrder(this->Order);
      q2D->Initialize(VTK_TRIANGLE);
      this->TensorProductWedge(q1D,q2D);
      q1D->Delete();
      q2D->Delete();
      break;
    }
    case VTK_TETRA:
    case VTK_QUADRATIC_TETRA:
    {
      if (this->Order == 0 || this->Order ==1)
      {
        this->QuadraturePoints->SetNumberOfComponents(3);
        this->QuadraturePoints->SetNumberOfTuples(1);
        this->QuadratureWeights->SetNumberOfTuples(1);
        double point[3];
        double weight;
        point[0] = 0.25;
        point[1] = 0.25;
        point[2] = 0.25;
        weight = .1666666666666666666666666666666666666666666667;
        this->QuadraturePoints->SetTuple(0,point);
        this->QuadratureWeights->SetValue(0,weight);
        break;
      }
      vtkfemriNSDQuadrature* gauss1D = vtkfemriNSDQuadrature::New();
      gauss1D->SetOrder(this->Order);
      gauss1D->Initialize1DGauss();
      vtkfemriNSDQuadrature* jacA1D = vtkfemriNSDQuadrature::New();
      jacA1D->SetOrder(this->Order);
      jacA1D->Initialize1DJacobi(1,0);
      vtkfemriNSDQuadrature* jacB1D = vtkfemriNSDQuadrature::New();
      jacB1D->SetOrder(this->Order);
      jacB1D->Initialize1DJacobi(2,0);
      this->TensorProductTetra(gauss1D,jacA1D,jacB1D);
      gauss1D->Delete();
      jacA1D->Delete();
      jacB1D->Delete();
      break;
    }
    default:
    {
      vtkErrorMacro("Unsupported element for Gauss quadrature.");
      return;
    }
  }
}

void vtkfemriNSDQuadrature::ScaleTo01()
{
  if (this->QuadraturePoints->GetNumberOfComponents() != 1)
  {
    vtkErrorMacro("Error: scaling assumes Dimensionality == 1.");
    return;
  }
 
  double point[1];
  double weight;
  int numberOfQuadraturePoints = this->GetNumberOfQuadraturePoints();
  int i;
  for (i=0; i<numberOfQuadraturePoints; i++)
  {
    this->QuadraturePoints->GetTuple(i,point);
    point[0] = 0.5 * (point[0] + 1.0);
    this->QuadraturePoints->SetTuple(i,point);

    weight = this->QuadratureWeights->GetValue(i);
    weight *= 0.5;
    this->QuadratureWeights->SetValue(i,weight);
  }
}

//alexmbcm: changes: added the following inline functions used in the file 
//vtkThinPlateSplineTransform.cxx:

//------------------------------------------------------------------------
// some dull matrix things

inline double** vtkNewMatrix(int rows, int cols) 
{
  double *matrix = new double[rows*cols];
  double **m = new double *[rows];
  for(int i = 0; i < rows; i++) 
    {
    m[i] = &matrix[i*cols];
    }
  return m;
}

//------------------------------------------------------------------------
inline void vtkDeleteMatrix(double **m) 
{
  delete [] *m;
  delete [] m;
}

//------------------------------------------------------------------------
inline void vtkZeroMatrix(double **m, int rows, int cols) 
{
  for(int i = 0; i < rows; i++) 
    {
    for(int j = 0; j < cols; j++) 
      {
      m[i][j] = 0.0;
      }
    }
}

//------------------------------------------------------------------------

//This function finds the quadrature points and weights for 
//Gauss-Legendre quadrature
//Interval [0,1]
//Weight function: 1
void vtkfemriNSDQuadrature::Initialize1DGauss()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();

	//alexmbcm: added Golub's method to get an arbitary number of gauss points
	int numPs;  //Number of points
	int i;      //Variable used for for loops
	//Zeroth moment of weight function 
	//(the integral of the weight function on quad interval) 
	double m_0 = 2.0; 
	double lengthSquared; //Variable used to calculate weight
  
    //Calculate max number of points needed based on max quad order 
	//needed
	if(this->Order % 2 == 0)
		numPs = (this->Order + 2)/2;
	else
		numPs = (this->Order + 1)/2;
	
	this->QuadraturePoints ->SetNumberOfTuples(numPs);
	this->QuadratureWeights->SetNumberOfTuples(numPs);
	
	//Fill points and weights using the Golub–Welsch algorithm
	double* points = this->QuadraturePoints->GetPointer(0);
	double* weights = this->QuadratureWeights->GetPointer(0);
    
	//Make matrix to apply Golub-Welsch algorithm
	//using the JacobiN function in vtkMath.h
	//Using inline functions above current function declaration
	double **A = vtkNewMatrix(numPs,numPs);
	vtkZeroMatrix(A,numPs,numPs);
		
	//Include in A the square root of 2 term recurrence coefficients
	//In the lower and upper diagonals
	for(i = 0; i < (numPs - 1); i++) 
		A[i+1][i] = A[i][i+1] = (i+1.0)/sqrt(4.0*(i+1.0)*(i+1.0) - 1.0);
		
	//Find points using JacobiN function
	double *qPoints       = new double[numPs];
	double **eigenVectors = vtkNewMatrix(numPs,numPs);
	vtkMath::JacobiN(A, numPs, qPoints, eigenVectors);
		
	//Put found points in point array taking into account 
	//they are ordered from positive to negative 
	for(i = 0; i < numPs; i++)
		points[i]  = qPoints[numPs-1 - i];
		
	//Calculate weights (first component of normalized eVect ^2 times fisrt moment of W)
	for(i = 0; i < numPs; i++)
		weights[i] = m_0*eigenVectors[0][numPs-1 - i]*eigenVectors[0][numPs-1 - i];
	
	//Delete matrices when done using them 
	vtkDeleteMatrix(A);
	delete [] qPoints; 
	vtkDeleteMatrix(eigenVectors);

  //Scale points and weights from [-1,1] to [0,1]
  this->ScaleTo01();
}

//This function finds the quadrature points and weights for 
//Gauss-Jacobi quadrature with beta = 0
//Interval [0,1]. Transformation done using a change of variables 
//For more information consult Gautschi's book:
//Numerical Analysis: An Introduction, section 2.4
//Weight function: (1-x)^alpha * (1+x)^beta
void vtkfemriNSDQuadrature::Initialize1DJacobi(int alpha, int beta)
{
  //Check if beta is 0
  if (beta != 0)
  {
	vtkErrorMacro("femri: Cannot generate Gauss-Jacobi quadrature with Beta != 0.");
    return; 
  }
  
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();
 
	//alexmbcm: added Golub's method to get an arbitary number of gauss points
	int numPs;  //Number of points
	int i;      //Variable used for for loops
	//Zeroth moment of weight function 
	//(just the integral of the weight function on quad interval) 
	double m_0 = 1.0/(alpha+1.0); 
	double lengthSquared; //Variable used to calculate weight
  
    //Calculate max number of points needed based on max quad order 
	//needed
	if(this->Order % 2 == 0)
		numPs = (this->Order + 2)/2;
	else
		numPs = (this->Order + 1)/2;
	
	this->QuadraturePoints ->SetNumberOfTuples(numPs);
	this->QuadratureWeights->SetNumberOfTuples(numPs);
	
	//Fill points and weights using the Golub–Welsch algorithm
	double* points = this->QuadraturePoints->GetPointer(0);
	double* weights = this->QuadratureWeights->GetPointer(0);
    
	//Make 0 matrix to apply Golub-Welsch algorithm
	//using the JacobiN function in vtkMath.h
	//Using inline functions above current function declaration
	double **A = vtkNewMatrix(numPs,numPs);
	vtkZeroMatrix(A,numPs,numPs);
		
	//Include in A the square root of 2 term recurrence coefficients
	//In the lower and upper diagonals
	for(i = 0; i < numPs; i++) 
	{
		A[i][i] = -alpha*alpha/(2.0*(i+1.0) + alpha)/(2.0*(i+1.0) - 2.0 + alpha);
			
		if(i < numPs-1)
			A[i+1][i] = A[i][i+1] = 
			sqrt( 4.0*(i+1.0)*(i+1.0)*(i+1.0+alpha)*(i+1.0+alpha) /
			((2.0*(i+1.0)+alpha)*(2.0*(i+1.0)+alpha)*(2.0*(i+1.0)+alpha+1.0)*(2.0*(i+1.0)+alpha-1.0)));
	}
		
	//Find points using JacobiN function
	double *qPoints       = new double[numPs];
	double **eigenVectors = vtkNewMatrix(numPs,numPs);
	vtkMath::JacobiN(A, numPs, qPoints, eigenVectors);
		
	//Put found points (from 0 to 1) in point array taking into account 
	//they are ordered from positive to negative 
	for(i = 0; i < numPs; i++)
		points[i] = (qPoints[numPs-1 - i] + 1.0)/2.0;
		
	//Calculate weights (first component of normalized eVect ^2 times first moment of W)
	for(i = 0; i < numPs; i++) 
		weights[i] = m_0*eigenVectors[0][numPs-1 - i]*eigenVectors[0][numPs-1 - i];
		
	//Delete matrices when done using them 
	vtkDeleteMatrix(A);
	delete [] qPoints; 
	vtkDeleteMatrix(eigenVectors);
}

void vtkfemriNSDQuadrature::TensorProductQuad(vtkfemriNSDQuadrature* q1D)
{
  int numberOf1DQuadraturePoints = q1D->GetNumberOfQuadraturePoints();

  this->QuadraturePoints->SetNumberOfComponents(2);
  this->QuadraturePoints->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);
  this->QuadratureWeights->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);

  double point[3];
  double weight;
  int id;
  int i, j;
  for (j=0; j<numberOf1DQuadraturePoints; j++)
  {
    for (i=0; i<numberOf1DQuadraturePoints; i++)
    {
      id = j * numberOf1DQuadraturePoints + i;
        
      point[0] = q1D->GetQuadraturePoint(i)[0];
      point[1] = q1D->GetQuadraturePoint(j)[0];
      point[2] = 0.0;

      weight = q1D->GetQuadratureWeight(i) * q1D->GetQuadratureWeight(j);
      
      this->QuadraturePoints->SetTuple(id,point);
      this->QuadratureWeights->SetValue(id,weight);
    }
  }
}
  
void vtkfemriNSDQuadrature::TensorProductTriangle(vtkfemriNSDQuadrature* gauss1D, vtkfemriNSDQuadrature* jacA1D)
{
  if (gauss1D->GetNumberOfQuadraturePoints() != jacA1D->GetNumberOfQuadraturePoints())
  {
    vtkErrorMacro("Error: cannot build tensor product rule if rules have different order.");
  }

  int numberOf1DQuadraturePoints = gauss1D->GetNumberOfQuadraturePoints();

  this->QuadraturePoints->SetNumberOfComponents(2);
  this->QuadraturePoints->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);
  this->QuadratureWeights->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);

  double point[2];
  double weight;
  int id;
  int i, j;
  
  // This loop computes the so called conical product rule for the standard triangle [0,0 1,0 0,1]
  // See also the code in LibMesh's src/quadrature/quadrature_conical.C
  for (j=0; j<numberOf1DQuadraturePoints; j++)
  {
    for (i=0; i<numberOf1DQuadraturePoints; i++)
    {
      id = j * numberOf1DQuadraturePoints + i;
        
      point[0] = jacA1D->GetQuadraturePoint(j)[0];
      point[1] = gauss1D->GetQuadraturePoint(i)[0] * (1.0 - jacA1D->GetQuadraturePoint(j)[0]);

      weight = gauss1D->GetQuadratureWeight(i) * jacA1D->GetQuadratureWeight(j);
      
      this->QuadraturePoints->SetTuple(id,point);
      this->QuadratureWeights->SetValue(id,weight);
    }
  }
}

void vtkfemriNSDQuadrature::TensorProductHexahedron(vtkfemriNSDQuadrature* q1D)
{
  int numberOf1DQuadraturePoints = q1D->GetNumberOfQuadraturePoints();

  this->QuadraturePoints->SetNumberOfComponents(3);
  this->QuadraturePoints->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);
  this->QuadratureWeights->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);

  double point[3];
  double weight;
  int id;
  int i, j, k;
  for (k=0; k<numberOf1DQuadraturePoints; k++)
  {
    for (j=0; j<numberOf1DQuadraturePoints; j++)
    {
      for (i=0; i<numberOf1DQuadraturePoints; i++)
      {
        id = k * numberOf1DQuadraturePoints * numberOf1DQuadraturePoints + j * numberOf1DQuadraturePoints + i;
          
        point[0] = q1D->GetQuadraturePoint(i)[0];
        point[1] = q1D->GetQuadraturePoint(j)[0];
        point[2] = q1D->GetQuadraturePoint(k)[0];
  
        weight = q1D->GetQuadratureWeight(i) * q1D->GetQuadratureWeight(j) * q1D->GetQuadratureWeight(k);
        
        this->QuadraturePoints->SetTuple(id,point);
        this->QuadratureWeights->SetValue(id,weight);
      }
    }
  }
}
  
void vtkfemriNSDQuadrature::TensorProductWedge(vtkfemriNSDQuadrature* q1D, vtkfemriNSDQuadrature* q2D)
{
  int numberOf1DQuadraturePoints = q1D->GetNumberOfQuadraturePoints();
  int numberOf2DQuadraturePoints = q2D->GetNumberOfQuadraturePoints();

  this->QuadraturePoints->SetNumberOfComponents(3);
  this->QuadraturePoints->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf2DQuadraturePoints);
  this->QuadratureWeights->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf2DQuadraturePoints);

  double point[3];
  double weight;
  int id;
  int i, j;
  for (j=0; j<numberOf1DQuadraturePoints; j++)
  {
    for (i=0; i<numberOf2DQuadraturePoints; i++)
    {
      id = j * numberOf2DQuadraturePoints + i;
        
      point[0] = q2D->GetQuadraturePoint(i)[0];
      point[1] = q2D->GetQuadraturePoint(i)[1];
      point[2] = q1D->GetQuadraturePoint(j)[0];

      weight = q2D->GetQuadratureWeight(i) * q1D->GetQuadratureWeight(j);
      
      this->QuadraturePoints->SetTuple(id,point);
      this->QuadratureWeights->SetValue(id,weight);
    }
  }
}
  
void vtkfemriNSDQuadrature::TensorProductTetra(vtkfemriNSDQuadrature* gauss1D, vtkfemriNSDQuadrature* jacA1D, vtkfemriNSDQuadrature* jacB1D)
{
  if (gauss1D->GetNumberOfQuadraturePoints() != jacA1D->GetNumberOfQuadraturePoints())
  {
    vtkErrorMacro("Error: cannot build tensor product rule if rules have different order.");
  }

  if (jacA1D->GetNumberOfQuadraturePoints() != jacB1D->GetNumberOfQuadraturePoints())
  {
    vtkErrorMacro("Error: cannot build tensor product rule if rules have different order.");
  }

  int numberOf1DQuadraturePoints = gauss1D->GetNumberOfQuadraturePoints();

  this->QuadraturePoints->SetNumberOfComponents(3);
  this->QuadraturePoints->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);
  this->QuadratureWeights->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);
 
  double point[3];
  double weight;
  int id;
  int i, j, k;
  for (k=0; k<numberOf1DQuadraturePoints; k++)
  {
    for (j=0; j<numberOf1DQuadraturePoints; j++)
    {
      for (i=0; i<numberOf1DQuadraturePoints; i++)
      {
        id = k * numberOf1DQuadraturePoints * numberOf1DQuadraturePoints + j * numberOf1DQuadraturePoints + i;
          
        point[0] = jacB1D->GetQuadraturePoint(k)[0];
        point[1] = jacA1D->GetQuadraturePoint(j)[0] * (1.0 - jacB1D->GetQuadraturePoint(k)[0]);
        point[2] = gauss1D->GetQuadraturePoint(i)[0] * (1.0 - jacA1D->GetQuadraturePoint(j)[0]) * (1.0 - jacB1D->GetQuadraturePoint(k)[0]);
  
        weight = gauss1D->GetQuadratureWeight(i) * jacA1D->GetQuadratureWeight(j) * jacB1D->GetQuadratureWeight(k);
        
        this->QuadraturePoints->SetTuple(id,point);
        this->QuadratureWeights->SetValue(id,weight);
      }
    }
  }
}

