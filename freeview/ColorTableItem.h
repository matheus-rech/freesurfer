#ifndef COLORTABLEITEM_H
#define COLORTABLEITEM_H

#include <QTreeWidgetItem>

class ColorTableItem : public QTreeWidgetItem
{
public:
  explicit ColorTableItem(int type = Type) : QTreeWidgetItem(type) {}
  explicit ColorTableItem(QTreeWidget* tree) : QTreeWidgetItem(tree) {}

  enum SORT_TYPE  { ST_VALUE = 0, ST_NAME };

  virtual bool operator < ( const QTreeWidgetItem& other ) const;

  static int  SortType;
  static bool SortAscending;
};

#endif // COLORTABLEITEM_H
