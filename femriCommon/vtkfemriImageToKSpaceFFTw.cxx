/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriImageToKSpaceFFTw.cxx,v $
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

#include "vtkfemriImageToKSpaceFFTw.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

//alexmbcm: to use use this file you need to install fftw
// http://www.fftw.org/
#include <fftw3.h>

#include <iostream>
using namespace std;


vtkStandardNewMacro(vtkfemriImageToKSpaceFFTw);
vtkCxxRevisionMacro(vtkfemriImageToKSpaceFFTw, "$Revision: 1.1.1.1 $");

vtkfemriImageToKSpaceFFTw::vtkfemriImageToKSpaceFFTw()
{
  this->Translation[0] = this->Translation[1] = this->Translation[2] = 0.0;
  this->SetNumberOfInputPorts(1);
  //TODO: force working with number of components = 2
}

vtkfemriImageToKSpaceFFTw::~vtkfemriImageToKSpaceFFTw()
{
}


int vtkfemriImageToKSpaceFFTw::RequestData(
    vtkInformation* vtkNotUsed( request ),
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  int updateExtent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),updateExtent);
  output->SetExtent(updateExtent);
  output->AllocateScalars();

  vtkDoubleArray* newScalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
  
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataArray* inScalars = input->GetPointData()->GetScalars();
  
  int kSpaceDimensions[3];
  double kSpaceSpacing[3];
  int extent[6];
  double* value;
  
  output->GetWholeExtent(extent);
  output->GetSpacing(kSpaceSpacing);
  output->GetDimensions(kSpaceDimensions);
  
  //cout << *(inScalars->GetTuple(26)) << endl;
  //cout << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] << endl;
  
	//Initialize the data array AFTER generating plan (data might be overwriten)
    int ijk[3];
	for (ijk[2]=extent[2]; ijk[2]<=extent[5]; ijk[2]++)
	{
		for (ijk[1]=extent[1]; ijk[1]<=extent[3]; ijk[1]++)
		{
			for (ijk[0]=extent[0]; ijk[0]<=extent[1]; ijk[0]++)
			{
				vtkIdType pointId = input->ComputePointId(ijk);
				value = inScalars->GetTuple(pointId);
				
				newScalars->SetComponent(pointId,0,*value);
				newScalars->SetComponent(pointId,1,0.0);
			}
		}
    }


  return 1;             
}

void vtkfemriImageToKSpaceFFTw::ShiftPhase(double value[2], double frequency[3], double translation[3], double shiftedValue[2])
{
  double phaseShift = -2.0 * vtkMath::Pi() * (frequency[0] * translation[0] + frequency[1] * translation[1] + frequency[2] * translation[2]);
  double shiftComplex[2];
  shiftComplex[0] = cos(phaseShift);
  shiftComplex[1] = sin(phaseShift);
  shiftedValue[0] = value[0] * shiftComplex[0] - value[1] * shiftComplex[1];
  shiftedValue[1] = value[0] * shiftComplex[1] + value[1] * shiftComplex[0];
}

void vtkfemriImageToKSpaceFFTw::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
