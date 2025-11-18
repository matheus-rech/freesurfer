#include "DialogLabelStatsInRange.h"
#include "ui_DialogLabelStatsInRange.h"
#include "LayerMRI.h"
#include <QClipboard>

DialogLabelStatsInRange::DialogLabelStatsInRange(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLabelStatsInRange)
{
  ui->setupUi(this);
}

DialogLabelStatsInRange::~DialogLabelStatsInRange()
{
  delete ui;
}

void DialogLabelStatsInRange::SetInfo(LayerMRI *mri, const QString &name, int val, int plane)
{
  int dim[3];
  double vs[3];
  mri->GetVolumeInfo(dim, vs);
  if (mri != m_mri || plane != m_nPlane)
  {
    ui->spinBoxStart->setRange(0, dim[plane]-1);
    ui->spinBoxStart->setValue(0);
    ui->spinBoxEnd->setRange(0, dim[plane]-1);
    ui->spinBoxEnd->setValue(dim[plane]-1);
  }
  m_mri = mri;
  m_strLabelName = name;
  m_nLabelValue = val;
  m_nPlane = plane;
  QStringList plane_names;
  plane_names << "Sagittal" << "Coronal" << "Axial";
  ui->labelMain->setText(QString("Stats of %1 between %2 slices:").arg(name).arg(plane_names[plane]));
  OnButtonUpdate();
}

void DialogLabelStatsInRange::OnButtonCopy()
{
  if (sender() == ui->pushButtonCopyVoxelCount)
    QApplication::clipboard()->setText(ui->labelVoxelCount->text());
  else
    QApplication::clipboard()->setText(ui->labelVolume->text().split(" ").first());
}

void DialogLabelStatsInRange::OnButtonUpdate()
{
  int dim[3];
  double vs[3];
  m_mri->GetVolumeInfo(dim, vs);
  int nCount = m_mri->GetLabelCount(m_nLabelValue, m_nPlane, ui->spinBoxStart->value(), ui->spinBoxEnd->value());
  ui->labelVoxelCount->setText(QString::number(nCount));
  ui->labelVolume->setText(QString("%1 mm3").arg(nCount*vs[0]*vs[1]*vs[2]));
}
