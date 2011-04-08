/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridNSDTime.h,v $
  Language:  C++
  Date:      $Date: 2008/11/04 11:23:41 $
  Version:   $Revision: 1.5 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkfemriUnstructuredGridToImageFilter - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Luca Antiga, PhD \n
// Imaging Research Labs - Robarts Research Institute \n
// Bioengineering Dept - Mario Negri Institute for Pharmacological Research \n 
// email: lantiga@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~lantiga

#ifndef __vtkfemriUnstructuredGridNSDTime_h
#define __vtkfemriUnstructuredGridNSDTime_h

#include "vtkfemriKSpaceGeneratorNSDTime.h"
#include "vtkfemriNumericWin32Header.h"

class vtkUnstructuredGrid;
class vtkCell;
class vtkfemriNSDQuadrature;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriUnstructuredGridNSDTime : public vtkfemriKSpaceGeneratorNSDTime
{
public:
  vtkTypeRevisionMacro(vtkfemriUnstructuredGridNSDTime,vtkfemriKSpaceGeneratorNSDTime);
  static vtkfemriUnstructuredGridNSDTime* New();

  vtkSetMacro(QuadratureOrder,int);
  vtkGetMacro(QuadratureOrder,int);

  vtkSetMacro(MagnetizationValue,double);
  vtkGetMacro(MagnetizationValue,double);  

  vtkGetMacro(NumberOfGaussPointEvaluations,int);
 
  virtual void Initialize();

  virtual void EvaluateFourierFunction(double frequency[3], double value[2], ofstream& writer);

protected:
  vtkfemriUnstructuredGridNSDTime();
  ~vtkfemriUnstructuredGridNSDTime();
  
  virtual int FillInputPortInformation(int, vtkInformation*);

  virtual double ComputeVolume()
  {
    return 1.0;
  }

  double ComputeJacobian(vtkCell* cell, double pcoords[3]);
  
  int QuadratureOrder;
  double MagnetizationValue;
  int NumberOfGaussPointEvaluations;
  
  vtkfemriNSDQuadrature* qLab;
  vtkfemriNSDQuadrature* qHab;
  vtkfemriNSDQuadrature* qGab;
 
 private:
  // Not implemented
  vtkfemriUnstructuredGridNSDTime(const vtkfemriUnstructuredGridNSDTime&);
  // Not implemented
  void operator=(const vtkfemriUnstructuredGridNSDTime&);
};

#endif
