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
#include "colortab.h"

class LUTDataHolder;
class QTreeWidgetItem;

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
  void OnColorTableItemChanged(QTreeWidgetItem *item);
  void OnCheckBoxSelectAllLabels(int nState);

private:
  void PopulateColorTable( COLOR_TABLE* ctab, bool bForce = false );

  Ui::PanelTrack *ui;

  QList<QWidget*> m_widgetlistDirectionalColor;
  QList<QWidget*> m_widgetlistSolidColor;
  QList<QWidget*> m_widgetlistScalarColor;
  QList<QWidget*> m_widgetlistScalarLut;
  QList<QWidget*> m_widgetlistScalarThreshold;

  LUTDataHolder* m_luts;
  COLOR_TABLE* m_curCTAB;
};

#endif // PANELTRACK_H
