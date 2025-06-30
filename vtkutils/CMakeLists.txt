project(vtkutils)

if(VTK_FOUND)

  if(NOT TIFF_SYSLIBS)
    include_directories(
      ${FS_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_SOURCE_DIR}/packages/tiff      
      ${VTK_INCLUDE_DIRS}
      ${X11_INCLUDE_DIR}
    )
  else()
    include_directories(
      ${FS_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${VTK_INCLUDE_DIRS}
      ${X11_INCLUDE_DIR}
    )    
  endif()

  set(SOURCES
    vtkFreesurferLookupTable.cxx
    vtkArrowPipeline.cxx
    vtkRGBATransferFunction.cxx
    vtkRGBAColorTransferFunction.cxx
    vtkInflatePolyData.cxx
    # vtkFDTensorGlyph.cxx
    # vtkODFGlyph.cxx
  )

  add_compile_options(-Wno-deprecated -Wno-inconsistent-missing-override)

  add_library(vtkutils STATIC ${SOURCES})
  target_link_libraries(vtkutils tiff)

endif()
