#ifndef DIALOGWAYPOINTSANGLE_H
#define DIALOGWAYPOINTSANGLE_H

#include <QDialog>

class Layer;
class LayerPointSet;

namespace Ui {
class DialogWayPointsAngle;
}

class DialogWayPointsAngle : public QDialog
{
  Q_OBJECT

public:
  explicit DialogWayPointsAngle(QWidget *parent = nullptr);
  ~DialogWayPointsAngle();

  bool eventFilter(QObject *watched, QEvent *event);

public slots:
  void OnSpinBoxValueChanged(int n);
  void SetPointSet(Layer* layer);
  void UpdateAngle();
  void OnButtonCopy();

private:
  void TriggerMoveToPoint(QObject* obj);

  Ui::DialogWayPointsAngle *ui;
  LayerPointSet*   m_pointset;
};

#endif // DIALOGWAYPOINTSANGLE_H
