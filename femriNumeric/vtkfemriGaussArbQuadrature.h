/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriGaussArbQuadrature.h,v $
  Language:  C++
  Date:      $Date: 2007/07/07 15:38:00 $
  Version:   $Revision: 1.2 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __vtkfemriGaussArbQuadrature_h
#define __vtkfemriGaussArbQuadrature_h

#include "vtkObject.h"
#include "vtkfemriNumericWin32Header.h"

#include "vtkDoubleArray.h"

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriGaussArbQuadrature : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkfemriGaussArbQuadrature,vtkObject);
  static vtkfemriGaussArbQuadrature* New();

  vtkGetObjectMacro(QuadraturePoints,vtkDoubleArray);
  vtkGetObjectMacro(QuadratureWeights,vtkDoubleArray);

  vtkSetMacro(Order,int);
  vtkGetMacro(Order,int);
  
  int GetNumberOfQuadraturePoints()
  {
    return this->QuadraturePoints->GetNumberOfTuples();
  }
  
  double* GetQuadraturePoint(vtkIdType id)
  {
    return this->QuadraturePoints->GetTuple(id);
  }
 
  void GetQuadraturePoint(vtkIdType id, double* quadraturePoint)
  {
    this->QuadraturePoints->GetTuple(id,quadraturePoint);
  }
  
  double GetQuadratureWeight(vtkIdType id)
  {
    return this->QuadratureWeights->GetValue(id);
  }
 
  void Initialize(vtkIdType cellType);
 
  void Initialize1DGauss();
  void Initialize1DJacobi(int alpha, int beta);
  void ScaleTo01();
 
protected:
  vtkfemriGaussArbQuadrature();
  ~vtkfemriGaussArbQuadrature();

  void TensorProductQuad(vtkfemriGaussArbQuadrature* q1D);
  void TensorProductTriangle(vtkfemriGaussArbQuadrature* gauss1D, vtkfemriGaussArbQuadrature* jacA1D);
  
  void TensorProductHexahedron(vtkfemriGaussArbQuadrature* q1D);
  void TensorProductWedge(vtkfemriGaussArbQuadrature* q1D, vtkfemriGaussArbQuadrature* q2D);
  void TensorProductTetra(vtkfemriGaussArbQuadrature* gauss1D, vtkfemriGaussArbQuadrature* jacA1D, vtkfemriGaussArbQuadrature* jacB1D);
  
  vtkDoubleArray* QuadraturePoints;
  vtkDoubleArray* QuadratureWeights;
 
  int Order;
  int QuadratureType;
  vtkIdType CellType;
  int PreviousOrder;

private:  
  vtkfemriGaussArbQuadrature(const vtkfemriGaussArbQuadrature&);  // Not implemented.
  void operator=(const vtkfemriGaussArbQuadrature&);  // Not implemented.

};

#endif
