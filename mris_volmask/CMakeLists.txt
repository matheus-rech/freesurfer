project(mris_volmask)

include_directories(${FS_INCLUDE_DIRS})

add_compile_options(-Wno-deprecated)

add_executable(mris_volmask mris_volmask.cpp cmd_line_interface.cpp)
target_compile_definitions(mris_volmask PRIVATE NO_VTK)
add_help(mris_volmask mris_volmask.help.xml)
target_link_libraries(mris_volmask utils)
install(TARGETS mris_volmask DESTINATION bin)

add_test_script(NAME mris_volmask_test SCRIPT test.sh)

add_executable(mris_volmask_novtk mris_volmask_old.cpp cmd_line_interface.cpp)
target_compile_definitions(mris_volmask_novtk PRIVATE NO_VTK)
add_help(mris_volmask_novtk mris_volmask.help.xml)
target_link_libraries(mris_volmask_novtk utils)
install(TARGETS mris_volmask_novtk DESTINATION bin)

if(VTK_FOUND AND ITK_FOUND AND NOT MINIMAL)

  include_directories(${VTK_INCLUDE_DIRS})
  add_executable(mris_volmask_vtk mris_volmask_old.cpp cmd_line_interface.cpp)

  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     set(UB24_BASE "/usr/lib/x86_64-linux-gnu")

     # reset and do not use default glob that pulls in every library
     set(VTK_LIBRARIES "")
     set(MRIS_VOLMASK_VTK_LIBS libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1 libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkFiltersHybrid-9.1.so.1)

     # reset and do not use default glob that pulls in every library
     set(ITK_LIBRARIES "")
     set(MRIS_VOLMASK_ITK_LIBS libitkv3p_netlib-5.3.so.1 libitkvnl-5.3.so.1 libitkvnl_algo-5.3.so.1)

     foreach(LIB ${MRIS_VOLMASK_ITK_LIBS} ${MRIS_VOLMASK_VTK_LIBS})
        target_link_libraries(mris_volmask_vtk ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(mris_volmask_vtk utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(mris_volmask_vtk utils ${VTK_LIBRARIES})
  endif()

  add_help(mris_volmask_vtk mris_volmask.help.xml)
  install(TARGETS mris_volmask_vtk DESTINATION bin)

  if(LINK_SHARED_AS_NEEDED)
     target_link_options(mris_volmask_vtk PRIVATE "LINKER:--as-needed")
  endif()

endif()
