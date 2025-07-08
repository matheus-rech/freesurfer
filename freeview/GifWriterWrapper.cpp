#include "GifWriterWrapper.h"
#include "gif.h"
#include <QImage>

GifWriterWrapper::GifWriterWrapper()
{
  m_gifWriter = NULL;
}

GifWriterWrapper::~GifWriterWrapper()
{
  EndGif();
}

void GifWriterWrapper::Initialize(const QString &fn, const QSize &sz, int nDelay)
{
  if (!m_gifWriter)
    m_gifWriter = new GifWriter;
  m_nDelay = nDelay/10;
  GifBegin((GifWriter*)m_gifWriter, qPrintable(fn), sz.width(), sz.height(), m_nDelay);
}

void GifWriterWrapper::AddToGif(const QImage &img, int nDelay_in)
{
  int nDelay = m_nDelay;
  if (nDelay_in > 0)
    nDelay = nDelay_in/10;
  GifWriter* gw = (GifWriter*)m_gifWriter;
  GifWriteFrame(gw, img.bits(), img.width(), img.height(), nDelay);
}

void GifWriterWrapper::AddToGif(const QString &fn, int nDelay_in)
{
  AddToGif(QImage(fn), nDelay_in);
}

void GifWriterWrapper::EndGif()
{
  if (m_gifWriter)
  {
    GifEnd((GifWriter*)m_gifWriter);
    delete (GifWriter*)m_gifWriter;
    m_gifWriter = NULL;
  }
}
