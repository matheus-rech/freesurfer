# third-party packages
project(packages)

# turn off all warnings and make sure libraries are position-independent
add_compile_options(-w -fPIC)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectories(
  jpeg
  expat
  glut
  xml2
  minc
  dicom
  dcm2niix
  cephes
  netcdf
  tetgen
  nrrdio
  nifti
  gifti
  svm
)
if(NOT TIFF_SYSLIBS)
  add_subdirectories(tiff)
endif()
