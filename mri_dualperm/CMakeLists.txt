project(mri_dual_perm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_dualperm mri_dualperm.cpp)
target_link_libraries(mri_dualperm utils)
install(TARGETS mri_dualperm DESTINATION bin)
