#include "DialogGifMaker.h"
#include "ui_DialogGifMaker.h"
#include <QFileDialog>
#include "MainWindow.h"
#include <QImage>
#include <QFile>
#include <QDateTime>
#include "RenderView.h"
#include "GifWriterWrapper.h"
#include <QMessageBox>


DialogGifMaker::DialogGifMaker(QWidget *parent) :
  QDialog(parent), m_nNumberOfFrames(0),
  ui(new Ui::DialogGifMaker)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );

  m_gif = new GifWriterWrapper;
  Reset();
}

DialogGifMaker::~DialogGifMaker()
{
  delete m_gif;
  delete ui;
}

void DialogGifMaker::hideEvent(QHideEvent *e)
{
  Reset();
}

void DialogGifMaker::OnButtonAdd()
{
  int w = ui->spinBoxWidth->value();
  int h;
  QString fn = QDir::tempPath() + "/freeview-temp-" + QString::number(QDateTime::currentMSecsSinceEpoch()) + ".png";
  MainWindow::GetMainWindow()->GetMainView()->SaveScreenShot(fn, false);
  QImage img(fn);
  if (ui->checkBoxRescale->isChecked())
  {
    img = img.scaledToWidth(w, Qt::SmoothTransformation);
  }
  w = img.width();
  h = img.height();
  int ndelay = ui->spinBoxDelay->value();
  if (m_nNumberOfFrames == 0)
  {
    m_strTempFilename = QDir::tempPath() + "/freeview-temp-" + QString::number(QDateTime::currentMSecsSinceEpoch()) + ".gif";
    m_gif->Initialize(m_strTempFilename, img.size(), ndelay);
    ui->pushButtonSave->setEnabled(true);
  }
  m_gif->AddToGif(img, ndelay);
  m_nNumberOfFrames++;
  ui->checkBoxRescale->setEnabled(false);
  ui->spinBoxWidth->setEnabled(false);
  ui->labelNumFrames->setText(QString::number(m_nNumberOfFrames));
}

void DialogGifMaker::Reset()
{
  if (m_nNumberOfFrames > 0)
    m_gif->EndGif();
  m_nNumberOfFrames = 0;
  ui->checkBoxRescale->setEnabled(true);
  ui->spinBoxWidth->setEnabled(ui->checkBoxRescale->isChecked());
  ui->labelNumFrames->setText("0");
  ui->pushButtonSave->setEnabled(false);
}

void DialogGifMaker::OnButtonSave()
{
  QString fn = QFileDialog::getSaveFileName(this, "Save to GIF",
                                            "",
                                            "GIF files (*.gif)");
  if (!fn.isEmpty())
  {
    if (m_nNumberOfFrames > 0)
      Reset();
    if (QFile::exists(fn))
      QFile::remove(fn);
    QFile::copy(m_strTempFilename, fn);
    hide();
  }
}

void DialogGifMaker::OnButtonClose()
{
  if (m_nNumberOfFrames > 0)
  {
    if (QMessageBox::question(this, "Warning", "Current sequence has not been saved. Do you want to close it without saving?",
                              QMessageBox::No, QMessageBox::Yes) != QMessageBox::Yes)
      return;
  }
  hide();
}
