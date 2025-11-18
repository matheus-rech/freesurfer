project(anatomicuts)

if(ITK_FOUND AND VTK_FOUND)
  include(${ITK_USE_FILE})  
  if(HOST_OS MATCHES "Rocky9" OR HOST_OS MATCHES "Ubuntu24")
     include_directories(${VTK_INCLUDE_DIRS})
  else()
     include(${VTK_USE_FILE})
  endif()

  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     set(UB24_BASE "/usr/lib/x86_64-linux-gnu")

     # set(UB24_COMMON_LINK_OPTIONS "-Wl,--unresolved-symbols=ignore-in-shared-libs -Wl,-rpath,/usr/lib/x86_64-linux-gnu/lib:/usr/lib/jvm/default-java/lib:/usr/lib/jvm/default-java/lib/server:/usr/lib/x86_64-linux-gnu/hdf5/serial")
     set(UB24_COMMON_LINK_OPTIONS "-Wl,--unresolved-symbols=ignore-in-shared-libs")

     # reset and do not use default glob that pulls in every library
     set(ITK_LIBRARIES "")

     set(ANATOMICUTS_COMMON_ITK_LIBS libITKCommon-5.3.so.1 libITKIOBMP-5.3.so.1 libITKIOBioRad-5.3.so.1 libITKIOBruker-5.3.so.1 libITKIOGDCM-5.3.so.1 libITKIOGE-5.3.so.1 libITKIOGIPL-5.3.so.1 libITKIOHDF5-5.3.so.1 libITKIOImageBase-5.3.so.1 libITKIOJPEG-5.3.so.1 libITKIOJPEG2000-5.3.so.1 libITKIOLSM-5.3.so.1 libITKIOMINC-5.3.so.1 libITKIOMRC-5.3.so.1 libITKIOMeta-5.3.so.1 libITKIONIFTI-5.3.so.1 libITKIONRRD-5.3.so.1 libITKIOPNG-5.3.so.1 libITKIOStimulate-5.3.so.1 libITKIOTIFF-5.3.so.1 libITKIOVTK-5.3.so.1 libITKNrrdIO-5.3.so.1 libITKniftiio-5.3.so.1 libITKznz-5.3.so.1 libitkdouble-conversion-5.3.so.1 libitkopenjpeg-5.3.so.1 libitksys-5.3.so.1 libitktiff-5.3.so.1 libitkvnl-5.3.so.1 libitkvnl_algo-5.3.so.1)

     set(DMRI_ANATOMICUTS_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS})
     set(DMRI_COLOREDFA_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS})
     set(DMRI_EXTRACTSURFACEMEASUREMENTS_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS} libitkv3p_netlib-5.3.so.1 libITKIOTransformHDF5-5.3.so.1 libITKIOTransformInsightLegacy-5.3.so.1 libITKIOTransformMatlab-5.3.so.1)
     set(DMRI_GROUPBYENDPOINTS_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS} libitkv3p_netlib-5.3.so.1)
     set(DMRI_MATCH_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS} libitkv3p_netlib-5.3.so.1)
     set(DMRI_NEIGHBORINGREGIONS_ITK_LIBS libitkdouble-conversion-5.3.so.1)
     set(DMRI_PROJECTENDPOINTS_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS} libitkv3p_netlib-5.3.so.1)
     set(DMRI_SAVEHISTOGRAMS_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS})
     set(DMRI_STATS_AC_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS})
     set(STREAMLINEFILTER_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS})
     set(TESTORIENTATIONPLANESFROMPARCELLATION_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS} libitkv3p_netlib-5.3.so.1)
     set(TRK_TOOLS_ITK_LIBS ${ANATOMICUTS_COMMON_ITK_LIBS})

     # reset and do not use default glob that pulls in every library
     set(VTK_LIBRARIES "")

     set(ANATOMICUTS_COMMON_VTK_LIBS_1 libITKIOVTK-5.3.so.1 libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1)
     set(ANATOMICUTS_COMMON_VTK_LIBS_2 libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1)

     set(DMRI_ANATOMICUTS_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkCommonSystem-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkIOCore-9.1.so.1 libvtkIOLegacy-9.1.so.1)
     set(DMRI_COLOREDFA_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkIOLegacy-9.1.so.1)
     set(DMRI_EXTRACTSURFACEMEASUREMENTS_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkIOLegacy-9.1.so.1)
     set(DMRI_GROUPBYENDPOINTS_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkIOLegacy-9.1.so.1)
     set(DMRI_MATCH_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersCore-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkFiltersSources-9.1.so.1 libvtkIOLegacy-9.1.so.1)
     # no vtk lib entries for DMRI_NEIGHBORINGREGIONS
     set(DMRI_PROJECTENDPOINTS_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkIOLegacy-9.1.so.1)
     set(DMRI_SAVEHISTOGRAMS_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersCore-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkFiltersSources-9.1.so.1 libvtkIOLegacy-9.1.so.1)
     set(DMRI_STATS_AC_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_1} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkIOLegacy-9.1.so.1)

     set(STREAMLINEFILTER_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_2} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkIOCore-9.1.so.1 libvtkIOLegacy-9.1.so.1)
     set(TESTORIENTATIONPLANESFROMPARCELLATION_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_2} libvtkCommonExecutionModel-9.1.so.1 libvtkFiltersCore-9.1.so.1 libvtkFiltersSources-9.1.so.1)
     set(TRK_TOOLS_VTK_LIBS ${ANATOMICUTS_COMMON_VTK_LIBS_2} libvtkIOCore-9.1.so.1 libvtkIOLegacy-9.1.so.1)
  else()
     set(ITK_LIBRARIES "ITKCommon;ITKIOImageBase;ITKIOBMP;ITKIOBioRad;ITKIOBruker;ITKIOGDCM;ITKIOGE;ITKIONRRD;ITKIOTIFF;ITKIOVTK;ITKIOStimulate;ITKIONIFTI;ITKIOMINC;ITKIOGIPL;ITKIOHDF5;ITKIOJPEG;ITKIOLSM;ITKIOMRC;ITKIOMeta;ITKIOPNG")
     if(HOST_OS MATCHES "CentOS9stream" OR HOST_OS MATCHES "Rocky9" OR (HOST_OS MATCHES "Ubuntu24" AND NOT UB24_SYSLIBS) OR (HOST_OS MATCHES "Ubuntu22" AND (HOST_OS_GCC_VERSION MATCHES "11" OR HOST_OS_GCC_VERSION MATCHES "13")))
       set(ITK_LIBRARIES "${ITK_LIBRARIES};ITKIOJPEG2000")
     elseif(HOST_OS MATCHES "Rocky8" AND HOST_OS_GCC_VERSION MATCHES "10")
       set(ITK_LIBRARIES "${ITK_LIBRARIES};ITKIOJPEG2000")
     endif()
  endif()

  set(TRACKIO ../freeview/track_io/TrackIO.cpp)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Code
    ${CMAKE_SOURCE_DIR}/freeview/track_io
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )

  if(NOT APPLE)
    add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedefs  -Wno-maybe-uninitialized )
  else()
    # Mac OS 10.13.6 native clang compiler does not support -Wno-maybe-uninitialized
    add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedefs -Wno-inconsistent-missing-override )
  endif()

#AnatomiCuts correspondences
  add_executable(dmri_AnatomiCuts AnatomiCuts.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_ANATOMICUTS_ITK_LIBS} ${DMRI_ANATOMICUTS_VTK_LIBS})
        target_link_libraries(dmri_AnatomiCuts ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_AnatomiCuts utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_AnatomiCuts utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_AnatomiCuts DESTINATION bin)

#AnatomiCuts
  add_executable(dmri_match AnatomiCuts_correspondences.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_MATCH_ITK_LIBS} ${DMRI_MATCH_VTK_LIBS})
        target_link_libraries(dmri_match ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_match utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_match utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_match DESTINATION bin)

#TRK tools for header
  add_executable(trk_tools trk_tools.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${TRK_TOOLS_ITK_LIBS} ${TRK_TOOLS_VTK_LIBS})
        target_link_libraries(trk_tools ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(trk_tools utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(trk_tools utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS trk_tools DESTINATION bin)

#streamline filtering
  add_executable(streamlineFilter StreamlineFilter.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${STREAMLINEFILTER_ITK_LIBS} ${STREAMLINEFILTER_VTK_LIBS})
        target_link_libraries(streamlineFilter ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(streamlineFilter utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(streamlineFilter utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS streamlineFilter DESTINATION bin)
 
#dmri_stats_ac
  add_executable(dmri_stats_ac dmri_stats_ac.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_STATS_AC_ITK_LIBS} ${DMRI_STATS_AC_VTK_LIBS})
        target_link_libraries(dmri_stats_ac ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_stats_ac utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_stats_ac utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_stats_ac DESTINATION bin)

#orientation planes
  add_executable(testOrientationPlanesFromParcellation testOrientationFromParcellationFilter.cxx ) 
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${TESTORIENTATIONPLANESFROMPARCELLATION_ITK_LIBS} ${TESTORIENTATIONPLANESFROMPARCELLATION_VTK_LIBS})
        target_link_libraries(testOrientationPlanesFromParcellation ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(testOrientationPlanesFromParcellation utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(testOrientationPlanesFromParcellation ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS testOrientationPlanesFromParcellation DESTINATION bin)

#Save histograms
  add_executable(dmri_saveHistograms SaveHistograms.cxx ${TRACKIO}) 
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_SAVEHISTOGRAMS_ITK_LIBS} ${DMRI_SAVEHISTOGRAMS_VTK_LIBS})
        target_link_libraries(dmri_saveHistograms ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_saveHistograms utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_saveHistograms ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_saveHistograms DESTINATION bin)

#dmri_groupByEndpoints
  add_executable(dmri_groupByEndpoints dmri_groupByEndpoints.cxx ${TRACKIO}) 
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_GROUPBYENDPOINTS_ITK_LIBS} ${DMRI_GROUPBYENDPOINTS_VTK_LIBS})
        target_link_libraries(dmri_groupByEndpoints ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_groupByEndpoints utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_groupByEndpoints utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_groupByEndpoints DESTINATION bin)

#dmri_extractSurfaceMeasurements
  # also needs ITKIOTransformHDF5;ITKIOTransformMatlab;ITKIOTransformInsightLegacy
  add_executable(dmri_extractSurfaceMeasurements dmri_extractSurfaceMeasurements.cxx ${TRACKIO}) 
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_EXTRACTSURFACEMEASUREMENTS_ITK_LIBS} ${DMRI_EXTRACTSURFACEMEASUREMENTS_VTK_LIBS})
        target_link_libraries(dmri_extractSurfaceMeasurements ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_extractSurfaceMeasurements utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_extractSurfaceMeasurements utils nifti "${ITK_LIBRARIES};ITKIOTransformHDF5;ITKIOTransformMatlab;ITKIOTransformInsightLegacy" ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_extractSurfaceMeasurements DESTINATION bin)

#dmri_coloredFA
  add_executable(dmri_coloredFA dmri_coloredFA.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_COLOREDFA_ITK_LIBS} ${DMRI_COLOREDFA_VTK_LIBS})
        target_link_libraries(dmri_coloredFA ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_coloredFA utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_coloredFA utils nifti ${OPENGL_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_coloredFA DESTINATION bin)

#dmri_projectEndPoints
  add_executable(dmri_projectEndPoints dmri_projectEndPoints.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     foreach(LIB ${DMRI_PROJECTENDPOINTS_ITK_LIBS} ${DMRI_PROJECTENDPOINTS_VTK_LIBS})
        target_link_libraries(dmri_projectEndPoints ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_projectEndPoints utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_projectEndPoints utils nifti ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_projectEndPoints DESTINATION bin)

#dmri_neighboringRegions
  # only needs ITKCommon
  add_executable(dmri_neighboringRegions dmri_neighboringRegions.cxx ${TRACKIO})
  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     # No entries for  DMRI_NEIGHBORINGREGIONS_VTK_LIBS
     foreach(LIB ${DMRI_NEIGHBORINGREGIONS_ITK_LIBS})
        target_link_libraries(dmri_neighboringRegions ${UB24_BASE}/${LIB})
     endforeach()
     target_link_libraries(dmri_neighboringRegions utils nifti)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UB24_COMMON_LINK_OPTIONS}")
  else()
     target_link_libraries(dmri_neighboringRegions utils nifti "ITKCommon" ${VTK_LIBRARIES})
  endif()
  install(TARGETS dmri_neighboringRegions DESTINATION bin)

   if(LINK_SHARED_AS_NEEDED)
     target_link_options(dmri_AnatomiCuts PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_match PRIVATE "LINKER:--as-needed")
     target_link_options(trk_tools PRIVATE "LINKER:--as-needed")
     target_link_options(streamlineFilter PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_stats_ac PRIVATE "LINKER:--as-needed")
     target_link_options(testOrientationPlanesFromParcellation PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_saveHistograms PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_groupByEndpoints PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_extractSurfaceMeasurements PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_coloredFA PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_projectEndPoints PRIVATE "LINKER:--as-needed")
     target_link_options(dmri_neighboringRegions PRIVATE "LINKER:--as-needed")
   endif()

endif()

install(PROGRAMS dmri_ac.sh DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(diffusionUtils anatomiCutsUtils dmri_violinPlots)
else()
   install_pyscript(diffusionUtils anatomiCutsUtils dmri_violinPlots)
endif()

