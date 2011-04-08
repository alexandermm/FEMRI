/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriNSDQuadrature.h,v $
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

#ifndef __vtkfemriNSDQuadrature_h
#define __vtkfemriNSDQuadrature_h

#include "vtkObject.h"
#include "vtkfemriNumericWin32Header.h"

#include "vtkDoubleArray.h"

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriNSDQuadrature : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkfemriNSDQuadrature,vtkObject);
  static vtkfemriNSDQuadrature* New();

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
  vtkfemriNSDQuadrature();
  ~vtkfemriNSDQuadrature();

  void TensorProductQuad(vtkfemriNSDQuadrature* q1D);
  void TensorProductTriangle(vtkfemriNSDQuadrature* gauss1D, vtkfemriNSDQuadrature* jacA1D);
  
  void TensorProductHexahedron(vtkfemriNSDQuadrature* q1D);
  void TensorProductWedge(vtkfemriNSDQuadrature* q1D, vtkfemriNSDQuadrature* q2D);
  void TensorProductTetra(vtkfemriNSDQuadrature* gauss1D, vtkfemriNSDQuadrature* jacA1D, vtkfemriNSDQuadrature* jacB1D);
  
  vtkDoubleArray* QuadraturePoints;
  vtkDoubleArray* QuadratureWeights;
 
  int Order;
  int QuadratureType;
  vtkIdType CellType;
  int PreviousOrder;

private:  
  vtkfemriNSDQuadrature(const vtkfemriNSDQuadrature&);  // Not implemented.
  void operator=(const vtkfemriNSDQuadrature&);  // Not implemented.

};

#endif
