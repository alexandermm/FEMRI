/*=========================================================================

  Program:   femri
  Module:    vtkfemriPolyDataQuadTriKSpaceGenerator.h
  Language:  C++
  Date:      $Date: 2006/04/06 16:46:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkfemriPolyDataQuadTriKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemriPolyDataQuadTriKSpaceGenerator_h
#define __vtkfemriPolyDataQuadTriKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriNumericWin32Header.h"

class vtkPolyData;
class vtkCell;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriPolyDataQuadTriKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemriPolyDataQuadTriKSpaceGenerator,vtkImageAlgorithm);
  static vtkfemriPolyDataQuadTriKSpaceGenerator* New();

  vtkSetMacro(QuadratureOrder,int);
  vtkGetMacro(QuadratureOrder,int);

  vtkSetMacro(MagnetizationValue,double);
  vtkGetMacro(MagnetizationValue,double);  

  virtual void Initialize();

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
  void EvaluateCellLocation(vtkCell* cell, double pcoords[3], double x[3], double* weights);
  double ComputeJacobian(vtkCell* cell, double pcoords[3]);
  
  
protected:
  vtkfemriPolyDataQuadTriKSpaceGenerator();
  ~vtkfemriPolyDataQuadTriKSpaceGenerator();
  
  virtual int FillInputPortInformation(int, vtkInformation*);

  virtual double ComputeVolume()
  {
    return 1.0;
  }
  
  int QuadratureOrder;
  double MagnetizationValue;

private:
  vtkfemriPolyDataQuadTriKSpaceGenerator(const vtkfemriPolyDataQuadTriKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriPolyDataQuadTriKSpaceGenerator&);  // Not implemented.
};

#endif
