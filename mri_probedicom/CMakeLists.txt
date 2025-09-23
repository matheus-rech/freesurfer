project(mri_probedicom)

if(OPENGL_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/packages/dicom
    ${GLUT_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/packages/dcm2niix
    ${X11_INCLUDE_DIR}
  )

  add_executable(mri_probedicom mri_probedicom.cpp)
  # dcm2niix issue867 reduce stack pressure (https://github.com/rordenlab/dcm2niix/commit/d8aacc58446c5dee08ebc8a0b83a264c401e8e88)
  #if(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  #  # -Wl,-dead_strip already added from top level
  #  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-stack_size -Wl,0x1000000")
  #endif()
  if(HOST_OS MATCHES "Rocky9" OR HOST_OS MATCHES "Ubuntu24")
     target_link_libraries(mri_probedicom utils ${OPENGL_LIBRARIES} -L/usr/lib64 -lXmu ${GLUT_LIBRARIES} ${X11_LIBRARIES})
  else()
     target_link_libraries(mri_probedicom utils ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${X11_LIBRARIES})
  endif()
  install(TARGETS mri_probedicom DESTINATION bin)

endif()
