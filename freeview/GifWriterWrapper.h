#ifndef GIFWRITERWRAPPER_H
#define GIFWRITERWRAPPER_H

#include <QString>
#include <QSize>

class QImage;

class GifWriterWrapper
{
public:
  GifWriterWrapper();
  ~GifWriterWrapper();

  void Initialize(const QString& fn, const QSize& sz, int nDelay);  // nDelay in ms

  void AddToGif(const QString& fn, int nDelay = 0);
  void AddToGif(const QImage& img, int nDelay = 0);
  void EndGif();

  void* m_gifWriter;
  int m_nDelay;
};

#endif // GIFWRITERWRAPPER_H
