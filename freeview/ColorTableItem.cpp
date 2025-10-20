#include "ColorTableItem.h"

int ColorTableItem::SortType = ColorTableItem::ST_VALUE;
bool ColorTableItem::SortAscending = true;

bool ColorTableItem::operator<(const QTreeWidgetItem &other) const
{
  QString txt = text(0);
  QString other_txt = other.text(0);
  bool bRet = false;
  if (SortType == ColorTableItem::ST_VALUE)
  {
    bRet = (data(0, Qt::UserRole+1).toInt() >
            other.data(0, Qt::UserRole+1).toInt());
  }
  else
  {
    //    if (txt.trimmed().contains(" "))
    //      txt = txt.split(" ", MD_SkipEmptyParts).at(1);
    //    if (other_txt.trimmed().contains(" "))
    //      other_txt = other_txt.split(" ", MD_SkipEmptyParts).at(1);
    if (txt.toLower() != other_txt.toLower())
    {
      txt = txt.toLower();
      other_txt = other_txt.toLower();
    }
    bRet = (txt > other_txt);
  }
  if (!SortAscending)
    bRet = !bRet;
  return bRet;
}
