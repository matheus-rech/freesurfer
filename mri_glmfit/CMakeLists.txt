project(mri_glmfit)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_glmfit mri_glmfit.cpp)
target_link_libraries(mri_glmfit utils)

add_test_script(NAME mri_glmfit_test SCRIPT test.sh)

install(TARGETS mri_glmfit DESTINATION bin RENAME mri_glmfit)

add_executable(mri_mvglmfit mri_mvglmfit.cpp)
target_link_libraries(mri_mvglmfit utils)
install(TARGETS mri_mvglmfit DESTINATION bin )
