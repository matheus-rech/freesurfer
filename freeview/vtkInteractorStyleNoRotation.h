#ifndef __vtkInteractorStyleNoRotation_h
#define __vtkInteractorStyleNoRotation_h

#include "vtkInteractorStyleTrackballCamera.h"


class vtkInteractorStyleNoRotation : public vtkInteractorStyleTrackballCamera
{
public:
  static vtkInteractorStyleNoRotation *New();
  vtkTypeMacro(vtkInteractorStyleNoRotation,vtkInteractorStyleTrackballCamera);

  void Rotate() override {}
  void Spin() override {}

protected:
 vtkInteractorStyleNoRotation();
 virtual ~vtkInteractorStyleNoRotation();

private:
  vtkInteractorStyleNoRotation(const vtkInteractorStyleNoRotation&);
  void operator=(const vtkInteractorStyleNoRotation&);
};

#endif
