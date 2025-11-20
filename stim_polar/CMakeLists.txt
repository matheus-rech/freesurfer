project(stim_polar)

if(BUILD_GUIS AND OPENGL_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${GLUT_INCLUDE_DIR} ${X11_INCLUDE_DIR})

  add_executable(stim_polar stim_polar.cpp)
  target_link_libraries(stim_polar utils ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${X11_LIBRARIES})
  install(TARGETS stim_polar DESTINATION bin)

endif()
