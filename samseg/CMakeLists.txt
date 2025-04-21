project(samseg)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${CMAKE_INSTALL_PREFIX}/\.\./fspython/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

install_configured(samseg samseg-long DESTINATION bin)
if(NOT INTEGRATE_SAMSEG)
   if (FSPYTHON_INSTALL_TREE)
      install_pyscript_fspython_tree(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
   else()
      install_pyscript(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
  endif()
endif()
add_test_script(NAME samseg_test SCRIPT test.sh)

