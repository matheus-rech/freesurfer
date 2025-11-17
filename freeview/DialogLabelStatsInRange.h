#ifndef DIALOGLABELSTATSINRANGE_H
#define DIALOGLABELSTATSINRANGE_H

#include <QDialog>

namespace Ui {
class DialogLabelStatsInRange;
}

class LayerMRI;

class DialogLabelStatsInRange : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLabelStatsInRange(QWidget *parent = nullptr);
  ~DialogLabelStatsInRange();

  void SetInfo(LayerMRI* mri, const QString& name, int val, int plane);

public slots:
  void OnButtonCopy();
  void OnButtonUpdate();

private:
  Ui::DialogLabelStatsInRange *ui;

  LayerMRI* m_mri;
  QString  m_strLabelName;
  int m_nLabelValue;
  int m_nPlane;
};

#endif // DIALOGLABELSTATSINRANGE_H
