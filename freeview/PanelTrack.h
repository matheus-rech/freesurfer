/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#ifndef PANELTRACK_H
#define PANELTRACK_H

#include "PanelLayer.h"
#include <QList>

class LUTDataHolder;

namespace Ui
{
class PanelTrack;
}

class PanelTrack : public PanelLayer
{
  Q_OBJECT

public:
  explicit PanelTrack(QWidget *parent = 0);
  ~PanelTrack();

protected:
  void DoUpdateWidgets();
  void DoIdle();
  virtual void ConnectLayer( Layer* layer );

protected slots:
  void OnSliderOpacity(int);
  void OnLineEditOpacity(const QString&);
  void OnSliderScalarThreshold(int);
  void OnLineEditScalarThreshold(const QString&);
  void OnComboLookupTable(int nSel);

private:
  Ui::PanelTrack *ui;

  QList<QWidget*> m_widgetlistDirectionalColor;
  QList<QWidget*> m_widgetlistSolidColor;
  QList<QWidget*> m_widgetlistScalarColor;
  QList<QWidget*> m_widgetlistScalarLut;
  QList<QWidget*> m_widgetlistScalarThreshold;

  LUTDataHolder* m_luts;
};

#endif // PANELTRACK_H
