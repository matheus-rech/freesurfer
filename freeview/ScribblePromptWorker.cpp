#include "ScribblePromptWorker.h"
#include "TorchScriptModule.h"
#include "LayerMRI.h"
#include "vtkImageData.h"
#include "vtkImageCast.h"
#include "vtkImageShiftScale.h"
#include "vtkImageCast.h"
#include "vtkImageInterpolator.h"
#include "vtkImageResize.h"
#include <QElapsedTimer>
#include <QDebug>
#include <QProcessEnvironment>
#include <QFile>

ScribblePromptWorker::ScribblePromptWorker(QObject *parent)
  : QObject(parent)
{
  m_module = new TorchScriptModule;
  connect(this, SIGNAL(ComputeTriggered()), SLOT(DoCompute()));
  connect(this, SIGNAL(ApplyTriggered()), SLOT(DoApply()));
  connect(this, SIGNAL(InitializationTriggered(QString)), SLOT(DoInitialization(QString)));
  // QString fn = QProcessEnvironment::systemEnvironment().value( "FREESURFER_HOME" ) + "/pytorch_models/traced_ScribblePrompt_UNet_nf192_res128.pt";
  // if (QFile::exists(fn))
  //   Initialize(fn);
  // else
  //   qDebug() << "Could not locate module file";
}

ScribblePromptWorker::~ScribblePromptWorker()
{
  m_module->deleteLater();
}

void ScribblePromptWorker::DoInitialization(const QString &fn)
{
  m_module->Load(fn);
  m_strModelFilename = fn;
}

void ScribblePromptWorker::Compute(LayerMRI *mri_ref, LayerMRI* seg, LayerMRI* seeds, int nPlane, int nSlice, double fill_val, bool include_existing, LayerMRI* mri_exit)
{
  m_ref = mri_ref;
  m_seg = seg;
  m_seeds = seeds;
  m_nInputPlane = nPlane;
  m_nInputSlice = nSlice;
  m_dFillValue = fill_val;
  m_bIncludeExisting = include_existing;
  m_curEdit = mri_exit;
  emit ComputeTriggered();
}

void ScribblePromptWorker::Apply(LayerMRI *seg, LayerMRI *filled, double fill_val)
{
  m_seg = seg;
  m_filled = filled;
  m_dFillValue = fill_val;
  emit ApplyTriggered();
}

vtkImageData* ScribblePromptWorker::GetResizedMriImage(float *ptr, int *dim, int *x_range, int *y_range, int* z_range, int nMag)
{
  vtkSmartPointer<vtkImageData> image_expand = vtkSmartPointer<vtkImageData>::New();
  image_expand->SetSpacing(1, 1, 1);
  int nMagSize = m_nMatSize*nMag, nSizeZ = m_b3D?nMagSize:1;
  image_expand->SetDimensions(nMagSize, nMagSize, nSizeZ);
  image_expand->AllocateScalars(VTK_FLOAT, 1);
  float* img_ptr = (float*)image_expand->GetScalarPointer();
  memset(img_ptr, 0, sizeof(float)*nMagSize*nMagSize*nSizeZ);
  for (int i = x_range[0]; i <= x_range[1]; i++)
  {
    for (int j = y_range[0]; j <= y_range[1]; j++)
    {
      for (int k = z_range[0]; k <= z_range[1]; k++)
      {
        img_ptr[(k-z_range[0])*nMagSize*nMagSize+(j-y_range[0])*nMagSize+(i-x_range[0])] = ptr[k*dim[0]*dim[1]+j*dim[0]+i];
      }
    }
  }
  vtkSmartPointer<vtkImageInterpolator> interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
  interpolator->SetInterpolationModeToCubic();
  vtkSmartPointer<vtkImageResize> resize = vtkSmartPointer<vtkImageResize>::New();
  resize->SetInputData(image_expand);
  resize->InterpolateOn();
  resize->SetInterpolator(interpolator);
  resize->SetResizeMethodToOutputDimensions();
  resize->SetOutputDimensions(m_nMatSize, m_nMatSize, m_b3D?m_nMatSize:1);
  resize->Update();
  vtkImageData* output = resize->GetOutput();
  output->SetReferenceCount(2);
  return output;
}

vtkImageData* ScribblePromptWorker::GetResizedSeedImage(unsigned char* ptr, int *dim, int *x_range, int *y_range, int* z_range, int nMag)
{
  vtkSmartPointer<vtkImageData> new_seed = vtkSmartPointer<vtkImageData>::New();
  new_seed->SetSpacing(1, 1, 1);
  int nSizeZ = m_b3D?m_nMatSize:1;
  new_seed->SetDimensions(m_nMatSize, m_nMatSize, nSizeZ);
  new_seed->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  unsigned char* new_ptr = (unsigned char*)new_seed->GetScalarPointer();
  memset(new_ptr, 0, m_nMatSize*m_nMatSize*nSizeZ);

  for (int n = 1; n <= 3; n++)
  {
    vtkSmartPointer<vtkImageData> image_expand = vtkSmartPointer<vtkImageData>::New();
    image_expand->SetSpacing(1, 1, 1);
    int nMagSize = m_nMatSize*nMag, nMagSizeZ = (m_b3D?nMagSize:1);
    image_expand->SetDimensions(nMagSize, nMagSize, nMagSizeZ);
    image_expand->AllocateScalars(VTK_FLOAT, 1);
    float* img_ptr = (float*)image_expand->GetScalarPointer();
    memset(img_ptr, 0, sizeof(float)*nMagSize*nMagSize*nMagSizeZ);
    for (int i = x_range[0]; i <= x_range[1]; i++)
    {
      for (int j = y_range[0]; j <= y_range[1]; j++)
      {
        for (int k = z_range[0]; k <= z_range[1]; k++)
        {
          if (ptr[k*dim[0]*dim[1]+j*dim[0]+i] == n)
            img_ptr[(k-z_range[0])*nMagSize*nMagSize+(j-y_range[0])*nMagSize+(i-x_range[0])] = 1;
        }
      }
    }
    vtkSmartPointer<vtkImageInterpolator> interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
    interpolator->SetInterpolationModeToCubic();
    vtkSmartPointer<vtkImageResize> resize = vtkSmartPointer<vtkImageResize>::New();
    resize->SetInputData(image_expand);
    resize->InterpolateOn();
    resize->SetInterpolator(interpolator);
    resize->SetResizeMethodToOutputDimensions();
    resize->SetOutputDimensions(m_nMatSize, m_nMatSize, nSizeZ);
    resize->Update();
    vtkSmartPointer<vtkImageData> output = resize->GetOutput();
    float* tmp_ptr = (float*)output->GetScalarPointer();
    for (int i = 0; i < m_nMatSize; i++)
    {
      for (int j = 0; j < m_nMatSize; j++)
      {
        for (int k = 0; k < nSizeZ; k++)
        {
          if (tmp_ptr[k*m_nMatSize*m_nMatSize+j*m_nMatSize+i] >= 1.0/nMag)
            new_ptr[k*m_nMatSize*m_nMatSize+j*m_nMatSize+i] = n;
        }
      }
    }
  }

  new_seed->SetReferenceCount(2);
  return new_seed;
}

void ScribblePromptWorker::ResizeImageData(float *ptr_in, int nx, int ny, int nz, float *ptr_out, int nx_out, int ny_out, int nz_out)
{
  vtkSmartPointer<vtkImageData> input = vtkSmartPointer<vtkImageData>::New();
  input->SetSpacing(1, 1, 1);
  input->SetDimensions(nx, ny, nz);
  input->AllocateScalars(VTK_FLOAT, 1);
  float* ptr = (float*)input->GetScalarPointer();
  memcpy(ptr, ptr_in, sizeof(float)*nx*ny*nz);
  // for (int i = 0; i < nx; i++)
  // {
  //   for (int j = 0; j < ny; j++)
  //   {
  //     for (int k = 0; k < nz; k++)
  //     {
  //       ptr[k*nx*ny+j*nx+i] = ptr_in[k*nx*ny+j*nx+i];
  //     }
  //   }
  // }
  vtkSmartPointer<vtkImageInterpolator> interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
  interpolator->SetInterpolationModeToCubic();
  vtkSmartPointer<vtkImageResize> resize = vtkSmartPointer<vtkImageResize>::New();
  resize->SetInputData(input);
  resize->InterpolateOn();
  resize->SetInterpolator(interpolator);
  resize->SetResizeMethodToOutputDimensions();
  resize->SetOutputDimensions(nx_out, ny_out, nz_out);
  resize->Update();
  ptr = (float*)resize->GetOutput()->GetScalarPointer();
  // for (int i = 0; i < nx_out; i++)
  // {
  //   for (int j = 0; j < ny_out; j++)
  //   {
  //      ptr_out[j*nx_out+i] = ptr[j*nx_out+i];
  //   }
  // }
  memcpy(ptr_out, ptr, sizeof(float)*nx_out*ny_out*nz_out);
}

void ScribblePromptWorker::DoCompute()
{
  QElapsedTimer timer;
  timer.start();
  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
  if (m_b3D)
    cast->SetInputData(m_ref->GetImageData());
  else
    cast->SetInputData(m_ref->GetSliceImageData(m_nInputPlane));
  cast->SetOutputScalarTypeToFloat();
  cast->Update();
  vtkSmartPointer<vtkImageData> img_mri = cast->GetOutput();
  vtkSmartPointer<vtkImageData> mri_seed = vtkSmartPointer<vtkImageData>::New();
  if (m_b3D)
    mri_seed->DeepCopy(m_seeds->GetImageData());
  else
    mri_seed->DeepCopy(m_seeds->GetSliceImageData(m_nInputPlane));
  int* dim = img_mri->GetDimensions();
  if (m_bIncludeExisting)
  {
    cast = vtkSmartPointer<vtkImageCast>::New();
    if (m_b3D)
      cast->SetInputData(m_curEdit->GetImageData());
    else
      cast->SetInputData(m_curEdit->GetSliceImageData(m_nInputPlane));
    cast->SetOutputScalarTypeToFloat();
    cast->Update();
    vtkSmartPointer<vtkImageData> mri_cur = cast->GetOutput();
    float* cur_seg_ptr = (float*)mri_cur->GetScalarPointer();
    unsigned char* seeds_ptr = (unsigned char*)mri_seed->GetScalarPointer();
    for (int i = 0; i < dim[0]; i++)
    {
      for (int j = 0; j < dim[1]; j++)
      {
        for (int k = 0; k < dim[2]; k++)
        {
          if (cur_seg_ptr[k*dim[0]*dim[1]+j*dim[0]+i] == m_dFillValue)
          {
            seeds_ptr[k*dim[0]*dim[1]+j*dim[0]+i] = 1;
          }
        }
      }
    }
  }
  float* mri_ptr = (float*)img_mri->GetScalarPointer();
  unsigned char* seeds_ptr = (unsigned char*)mri_seed->GetScalarPointer();
  bool bOverSize = (dim[0] > m_nMatSize || dim[1] > m_nMatSize);
  if (m_b3D)
    bOverSize = (bOverSize || dim[2] > m_nMatSize);
  int x_range[2] = {1000000,-1000000}, y_range[2] = {1000000,-1000000}, z_range[2] = {1000000,-1000000};
  int start_x = 0, start_y = 0, start_z = 0;
  int nMag = 1;
  vtkImageData *new_mri = NULL, *new_seed = NULL;
  if (bOverSize)
  {
    for (int i = 0; i < dim[0]; i++)
    {
      for (int j = 0; j < dim[1]; j++)
      {
        for (int k = 0; k < dim[2]; k++)
        {
          if (seeds_ptr[k*dim[0]*dim[1]+j*dim[0]+i] > 0)
          {
            if (i < x_range[0])
              x_range[0] = i;
            if (i > x_range[1])
              x_range[1] = i;
            if (j < y_range[0])
              y_range[0] = j;
            if (j > y_range[1])
              y_range[1] = j;
            if (k < z_range[0])
              z_range[0] = k;
            if (k > z_range[1])
              z_range[1] = k;
          }
        }
      }
    }

//    qDebug() << x_range[0] << x_range[1] << y_range[0] << y_range[1] << z_range[0] << z_range[1];

    if (x_range[1]-x_range[0] < 0 || y_range[1]-y_range[0] < 0 || z_range[1]-z_range[0] < 0)
    {
      qDebug() << "No foreground or background seeds selected";
      emit ComputeFinished(timer.elapsed()/1000.0);
      return;
    }
    else if (x_range[1]-x_range[0] > m_nMatSize || y_range[1]-y_range[0] > m_nMatSize ||
             (m_b3D && z_range[1]-z_range[0] > m_nMatSize))
    {
      nMag = qMax((x_range[1]-x_range[0])/m_nMatSize, (y_range[1]-y_range[0])/m_nMatSize);
      if (m_b3D)
        nMag = qMax(nMag, (z_range[1]-z_range[0])/m_nMatSize);
      nMag = nMag + 1;
      new_mri = GetResizedMriImage(mri_ptr, dim, x_range, y_range, z_range, nMag);
      new_seed = GetResizedSeedImage(seeds_ptr, dim, x_range, y_range, z_range, nMag);
      dim = new_mri->GetDimensions();
      mri_ptr = (float*)new_mri->GetScalarPointer();
      seeds_ptr = (unsigned char*)new_seed->GetScalarPointer();
    }
    else
    {
      start_x = qMax(0, (x_range[1]+x_range[0])/2-m_nMatSize/2);
      start_y = qMax(0, (y_range[1]+y_range[0])/2-m_nMatSize/2);
      if (m_b3D)
        start_z = qMax(0, (z_range[1]+z_range[0])/2-m_nMatSize/2);
    }
  }

  QVector<float*> inputs;
  int nSizeZ = (m_b3D?m_nMatSize:1);
  for (int i = 0; i < 5; i++)
  {
    float* ptr = new float[m_nMatSize*m_nMatSize*nSizeZ];
    memset(ptr, 0, sizeof(float)*m_nMatSize*m_nMatSize*nSizeZ);
    inputs << ptr;
  }

  double value_r[2];
  img_mri->GetScalarRange(value_r);
  for (int i = 0; i < m_nMatSize; i++)
  {
    for (int j = 0; j < m_nMatSize; j++)
    {
      for (int k = 0; k < nSizeZ; k++)
      {
        int x = start_x + i, y = start_y + j, z = start_z + k;
        if (x >= dim[0] || y >= dim[1] || z >= dim[2])
          continue;

        inputs[0][k*m_nMatSize*m_nMatSize+j*m_nMatSize+i] = mri_ptr[z*dim[0]*dim[1]+y*dim[0]+x]/value_r[1];
        if (seeds_ptr[z*dim[0]*dim[1]+y*dim[0]+x] == 1)       // foreground
          inputs[2][k*m_nMatSize*m_nMatSize+j*m_nMatSize+i] = 1;
        else if (seeds_ptr[z*dim[0]*dim[1]+y*dim[0]+x] == 2)    // background
          inputs[3][k*m_nMatSize*m_nMatSize+j*m_nMatSize+i] = 1;
        else if (seeds_ptr[z*dim[0]*dim[1]+y*dim[0]+x] == 3)   // box
          inputs[1][k*m_nMatSize*m_nMatSize+j*m_nMatSize+i] = 1;
      }
    }
  }

  float* output = new float[m_nMatSize*m_nMatSize*nSizeZ];
  m_module->Run(inputs, output, m_nMatSize, m_b3D);
  int nMagSize = m_nMatSize*nMag;
  int nMagSizeZ = (m_b3D?nMagSize:1);
  if (nMag > 1)
  {
    float* new_output = new float[nMagSize*nMagSize*nMagSizeZ];
    ResizeImageData(output, m_nMatSize, m_nMatSize, nSizeZ, new_output, nMagSize, nMagSize, nMagSizeZ);
    float* old = output;
    output = new_output;
    delete[] old;
    start_x = x_range[0];
    start_y = y_range[0];
    start_z = z_range[0];
  }
  void* p = m_seg->GetImageData()->GetScalarPointer();
  int nDataType = m_seg->GetImageData()->GetScalarType();
  double fillValue = m_seg->GetFillValue();
  dim = m_seg->GetImageData()->GetDimensions();
  for (int i = 0; i < nMagSize; i++)
  {
    for (int j = 0; j < nMagSize; j++)
    {
      for (int k = 0; k < nMagSizeZ; k++)
      {
        if (output[k*nMagSize*nMagSize + j*nMagSize + i] <= 0)
          continue;

        int x = start_x + i, y = start_y + j, z = start_z + k;
        int n = 0;
        if (m_b3D)
        {
          if (x >= dim[0] || y >= dim[1] || z >= dim[2])
            continue;
          n = z*dim[1]*dim[0] + y*dim[0] + x;
        }
        else
        {
          switch (m_nInputPlane)
          {
          case 0:
            if (y >= dim[2] || x >= dim[1])
              continue;
            n = y*dim[1]*dim[0] + x*dim[0] + m_nInputSlice;
            break;
          case 1:
            if (y >= dim[2] || x >= dim[0])
              continue;
            n = y*dim[1]*dim[0] + m_nInputSlice*dim[0] + x;
            break;
          case 2:
            if (y >= dim[1] || x >= dim[0])
              continue;
            n = m_nInputSlice*dim[1]*dim[0] + y*dim[0] + x;
            break;
          }
        }

        switch (nDataType)
        {
        case VTK_INT:
          ((int*)p)[n] = (int)fillValue;
          break;
        case VTK_UNSIGNED_CHAR:
          ((unsigned char*)p)[n] = (unsigned char)fillValue;
          break;
        case VTK_FLOAT:
          ((float*)p)[n] = (float)fillValue;
          break;
        case VTK_DOUBLE:
          ((double*)p)[n] = (double)fillValue;
          break;
        }
      }
    }
  }

  for (int i = 0; i < inputs.size(); i++)
    delete[] inputs[i];
  delete[] output;

  if (new_seed)
    new_seed->Delete();
  if (new_mri)
    new_mri->Delete();

  m_seg->SetModified();
  emit ComputeFinished(timer.elapsed()/1000.0);
}

void ScribblePromptWorker::DoApply()
{
  m_seg->SaveForUndo();
  void* p = m_seg->GetImageData()->GetScalarPointer();
  int nDataType = m_seg->GetImageData()->GetScalarType();
  double fillValue = m_seg->GetFillValue();
  int* dim = m_filled->GetImageData()->GetDimensions();
  unsigned char* p_filled = (unsigned char*)m_filled->GetImageData()->GetScalarPointer();
  size_t vol_size = dim[0]*dim[1]*dim[2];
  for (size_t i = 0; i < vol_size; i++)
  {
    if (p_filled[i])
    {
      switch (nDataType)
      {
      case VTK_INT:
        ((int*)p)[i] = (int)fillValue;
        break;
      case VTK_UNSIGNED_CHAR:
        ((unsigned char*)p)[i] = (unsigned char)fillValue;
        break;
      case VTK_FLOAT:
        ((float*)p)[i] = (float)fillValue;
        break;
      case VTK_DOUBLE:
        ((double*)p)[i] = (double)fillValue;
        break;
      }
    }
  }
  m_seg->SetModified();
  emit ApplyFinished();
}
