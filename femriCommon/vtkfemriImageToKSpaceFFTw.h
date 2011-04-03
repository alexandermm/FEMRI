/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriImageToKSpaceFFTw.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkfemriImageToKSpaceFFTw - ..
// .SECTION Description
// ..

#ifndef __vtkfemriImageToKSpaceFFTw_h
#define __vtkfemriImageToKSpaceFFTw_h

#include "vtkImageAlgorithm.h"
#include "vtkfemriCommonWin32Header.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemriImageToKSpaceFFTw : public vtkImageAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkfemriImageToKSpaceFFTw,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemriImageToKSpaceFFTw* New();

  vtkSetVectorMacro(Translation,double,3);
  vtkGetVectorMacro(Translation,double,3);

  static void ShiftPhase(double value[2], double frequency[3], double translation[3], double shiftedValue[2]);
  
protected:
  vtkfemriImageToKSpaceFFTw();
  ~vtkfemriImageToKSpaceFFTw();

  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  double Translation[3];

private:
  vtkfemriImageToKSpaceFFTw(const vtkfemriImageToKSpaceFFTw&);  // Not implemented.
  void operator=(const vtkfemriImageToKSpaceFFTw&);  // Not implemented.
};

#endif
