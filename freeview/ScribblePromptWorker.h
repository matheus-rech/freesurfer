#ifndef SCRIBBLEPROMPTWORKER_H
#define SCRIBBLEPROMPTWORKER_H

#include <QObject>
#include "vtkSmartPointer.h"

class TorchScriptModule;
class LayerMRI;
class vtkImageData;

class ScribblePromptWorker : public QObject
{
  Q_OBJECT
public:
  explicit ScribblePromptWorker(QObject *parent = nullptr);
  ~ScribblePromptWorker();

  void Initialize(const QString& fn)
  {
    emit InitializationTriggered(fn);
  }

  QString GetModelFilename()
  {
    return m_strModelFilename;
  }

signals:
  void InitializationTriggered(const QString& fn);
  void ComputeTriggered();
  void ApplyTriggered();
  void ComputeFinished(double elapsed_time);
  void ApplyFinished();

public slots:
  void Compute(LayerMRI *mri_ref, LayerMRI* seg, LayerMRI* seeds, int nPlane, int nSlice, double fill_val, bool include_existing, LayerMRI* mri_edit);
  void Apply(LayerMRI *seg, LayerMRI *filled, double fill_val);
  void LoadModel(const QString& fn, int mat_size, bool b3D = false)
  {
    m_b3D = b3D;
    m_nMatSize = mat_size;
    Initialize(fn);
  }
  void Set3D(bool b)
  {
    m_b3D = b;
  }
  void SetMatrixSize(int n)
  {
    m_nMatSize = n;
  }

private slots:
  void DoInitialization(const QString& fn);
  void DoCompute();
  void DoApply();

private:
  vtkImageData* GetResizedMriImage(float* ptr, int* dim, int* x_range, int* y_range, int* z_range, int nMag);
  vtkImageData* GetResizedSeedImage(unsigned char* ptr, int* dim, int* x_range, int* y_range, int* z_range, int nMag);
  void ResizeImageData(float* ptr_in, int nx, int ny, int nz, float* ptr_out, int nx_out, int ny_out, int nz_out);

  LayerMRI* m_seeds;
  LayerMRI* m_ref;
  LayerMRI* m_seg;
  LayerMRI* m_filled;
  LayerMRI* m_curEdit;
  int     m_nInputPlane;
  int     m_nInputSlice;
  double    m_dFillValue;
  bool    m_bIncludeExisting;
  QString m_strModelFilename;
  TorchScriptModule* m_module;
  bool    m_b3D;
  int     m_nMatSize;
};

#endif // SCRIBBLEPROMPTWORKER_H
