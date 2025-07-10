#include "DialogWayPointsAngle.h"
#include "ui_DialogWayPointsAngle.h"
#include "LayerPointSet.h"
#include "MainWindow.h"
#include <QApplication>
#include <QClipboard>
#include "MyUtils.h"
#include "vtkMath.h"
#include <QDebug>

DialogWayPointsAngle::DialogWayPointsAngle(QWidget *parent) :
  QDialog(parent), m_pointset(NULL),
  ui(new Ui::DialogWayPointsAngle)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Window | Qt::WindowTitleHint | Qt::WindowStaysOnTopHint |
                        Qt::CustomizeWindowHint | Qt::WindowCloseButtonHint );

  ui->spinBoxA->installEventFilter(this);
  ui->spinBoxB->installEventFilter(this);
  ui->spinBoxC->installEventFilter(this);
}

DialogWayPointsAngle::~DialogWayPointsAngle()
{
  delete ui;
}

bool DialogWayPointsAngle::eventFilter(QObject *watched, QEvent *event)
{
  if (event->type() == QEvent::FocusIn)
  {
    TriggerMoveToPoint(watched);
  }

  return QDialog::eventFilter(watched, event);
}

void DialogWayPointsAngle::SetPointSet(Layer* layer)
{
  LayerPointSet* ps = qobject_cast<LayerPointSet*>(layer);
  m_pointset = ps;
  if (!ps)
  {
    hide();
    return;
  }

  ui->spinBoxA->setRange(1, ps->GetNumberOfPoints());
  ui->spinBoxB->setRange(1, ps->GetNumberOfPoints());
  ui->spinBoxC->setRange(1, ps->GetNumberOfPoints());

  UpdateAngle();
}

void DialogWayPointsAngle::OnSpinBoxValueChanged(int n)
{
  TriggerMoveToPoint(sender());
  UpdateAngle();
}

void DialogWayPointsAngle::TriggerMoveToPoint(QObject *obj)
{
  QSpinBox* spinbox = qobject_cast<QSpinBox*>(obj);
  if (spinbox && m_pointset)
  {
    double pt[3];
    if (m_pointset->GetPoint(spinbox->value()-1, pt))
      MainWindow::GetMainWindow()->SetSlicePosition(pt);
  }
}

void DialogWayPointsAngle::UpdateAngle()
{
  if (!m_pointset)
    return;

  int nA = ui->spinBoxA->value()-1;
  int nB = ui->spinBoxB->value()-1;
  int nC = ui->spinBoxC->value()-1;
  if (nB == nA || nB == nC)
  {
    ui->labelAngle->setText("N/A");
    return;
  }
  double ptA[3], ptB[3], ptC[3], v1[3], v2[3];
  m_pointset->GetPoint(nA, ptA);
  m_pointset->GetPoint(nB, ptB);
  m_pointset->GetPoint(nC, ptC);

  MyUtils::GetVector(ptB, ptA, v1, true);
  MyUtils::GetVector(ptB, ptC, v2, true);

  double val = vtkMath::AngleBetweenVectors(v1, v2)*180/vtkMath::Pi();
  ui->labelAngle->setText(QString("%1").arg(val, 0, 'f', 2));
}

void DialogWayPointsAngle::OnButtonCopy()
{
  QClipboard* clipboard = QApplication::clipboard();
  clipboard->setText(ui->labelAngle->text());
}
