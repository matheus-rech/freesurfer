project(mri_volcluster)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_volcluster mri_volcluster.cpp)
target_link_libraries(mri_volcluster utils)

install(TARGETS mri_volcluster DESTINATION bin)

add_executable(mri_cluster mri_cluster.cpp)
target_link_libraries(mri_cluster utils)
install(TARGETS mri_cluster DESTINATION bin)
