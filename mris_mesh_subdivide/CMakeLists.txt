project(mris_mesh_subdivide)

if(VTK_FOUND AND ITK_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

  add_compile_options(-Wno-deprecated)

  add_executable(mris_mesh_subdivide mris_mesh_subdivide.cxx)

  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     set(UB24_BASE "/usr/lib/x86_64-linux-gnu")

     # reset and do not use default glob that pulls in every library
     set(VTK_LIBRARIES "")
     set(MRIS_MESH_SUBDIVIDE_VTK_LIBS libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1 libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersModeling-9.1.so.1)

     # reset and do not use default glob that pulls in every library
     set(ITK_LIBRARIES "")
     set(MRIS_MESH_SUBDIVIDE_ITK_LIBS libitkv3p_netlib-5.3.so.1 libitkvnl-5.3.so.1 libitkvnl_algo-5.3.so.1)

     foreach(LIB ${MRIS_MESH_SUBDIVIDE_ITK_LIBS} ${MRIS_MESH_SUBDIVIDE_VTK_LIBS})
        target_link_libraries(mris_mesh_subdivide ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(mris_mesh_subdivide utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(mris_mesh_subdivide utils ${VTK_LIBRARIES})
  endif()

  if(LINK_SHARED_AS_NEEDED)
     target_link_options(mris_mesh_subdivide PRIVATE "LINKER:--as-needed")
  endif()

  install(TARGETS mris_mesh_subdivide DESTINATION bin)

endif()
