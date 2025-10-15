#include "TorchScriptModule.h"
#undef slots
#undef Byte
#include "torch/script.h"
#define slots Q_SLOTS
#include <QDebug>

TorchScriptModule::TorchScriptModule(QObject *parent)
    : QObject(parent)
{
  m_module = new torch::jit::script::Module;
}

TorchScriptModule::~TorchScriptModule()
{
  delete ((torch::jit::script::Module*)m_module);
}

void TorchScriptModule::Load(const QString &fn)
{
  *((torch::jit::script::Module*)m_module) = torch::jit::load(qPrintable(fn));
}

void FillInTensorFromBuffer(torch::Tensor &t, int n, float *buf_in, int mat_size, bool b3D)
{
  if (b3D)
  {
    auto ptr = t.accessor<float,5>();
    for (int i = 0; i < mat_size; i++)
    {
      for (int j = 0; j < mat_size; j++)
      {
        for (int k = 0; k < mat_size; k++)
        {
          ptr[0][n][k][j][i] = buf_in[k*mat_size*mat_size+j*mat_size+i];
        }
      }
    }
  }
  else
  {
    auto ptr = t.accessor<float,4>();
    for (int i = 0; i < mat_size; i++)
    {
      for (int j = 0; j < mat_size; j++)
      {
        ptr[0][n][j][i] = buf_in[j*mat_size+i];
      }
    }
  }
}

void FillInBufferFromTensor(float *buffer, torch::Tensor &t_in, int n, int mat_size, bool b3D)
{
  if (b3D)
  {
    auto tptr = t_in.accessor<float,5>();
    for (int i = 0; i < mat_size; i++)
    {
      for (int j = 0; j < mat_size; j++)
      {
        for (int k = 0; k < mat_size; k++)
        {
          buffer[k*mat_size*mat_size+j*mat_size+i] = tptr[0][n][k][j][i];
        }
      }
    }
  }
  else
  {
    auto tptr = t_in.accessor<float,4>();
    for (int i = 0; i < mat_size; i++)
    {
      for (int j = 0; j < mat_size; j++)
      {
        buffer[j*mat_size+i] = tptr[0][n][j][i];
      }
    }
  }
}

void TorchScriptModule::Run(QVector<float*> in_ptr, float *output, int mat_size, bool b3D)
{
  try {
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor t;
    if (b3D)
      t = torch::zeros({1, in_ptr.size(), mat_size, mat_size, mat_size});
    else
      t = torch::zeros({1, in_ptr.size(), mat_size, mat_size});
    for (int i = 0; i < in_ptr.size(); i++)
    {
      if (in_ptr[i])
        FillInTensorFromBuffer(t, i, in_ptr[i], mat_size, b3D);
    }
    inputs.push_back(t);
    // Execute the model and turn its output into a tensor.
    torch::Tensor t_out = ((torch::jit::script::Module*)m_module)->forward(inputs).toTensor();
    FillInBufferFromTensor(output, t_out, 0, mat_size, b3D);
  }
  catch (const c10::Error& e) {
    qDebug() << "error running the module\n";
  }
  emit Finished();
}

