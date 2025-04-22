project(samseg)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${CMAKE_INSTALL_PREFIX}/\.\./fspython/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

function(install_samseg_cli)
  # ??? todo: retrieve samseg directory using get_package_dir(samseg samsegdir) ???
  # set SAMSEGCLI_DIR, this is hard-coded for now
  if(NOT "${USER_SAMSEG_PATH}" STREQUAL "")
    set(SAMSEGCLI_DIR ${USER_SAMSEG_PATH})
  else()
    set(SAMSEGCLI_DIR ${FSPYTHON_INSTALL_PREFIX}/python/packages)
  endif()

  set(SAMSEG_CLI_SCRIPTS
    run_samseg
    run_samseg_long
    sbtiv
    gems_compute_atlas_probs
    segment_subregions
    prepareAtlasDirectory
    merge_add_mesh_alphas
    computeTissueConcentrations)

  # config samseg cli scripts
  foreach(SCRIPT ${SAMSEG_CLI_SCRIPTS})
    # config fspython wrapper script to wrap around ${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py
    if (FSPYTHON_INSTALL_TREE)
      # modified from config_fspythonwrapper_fspython_tree()
      install(CODE "
        message(STATUS \"Configuring fspython wrapper: ${CMAKE_INSTALL_PREFIX}/../fspython/bin/${SCRIPT}\")
        set(SCRIPTNAME ${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py)
        configure_file(${CMAKE_SOURCE_DIR}/samseg/samseg_cli_wrapper ${CMAKE_INSTALL_PREFIX}/../fspython/bin/${SCRIPT} @ONLY)"
      )      
    else()
      # modified from config_fspythonwrapper()
      install(CODE "
        message(STATUS \"Configuring fspython wrapper: ${CMAKE_INSTALL_PREFIX}/bin/${SCRIPT}\")
        set(SCRIPTNAME ${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py)
        configure_file(${CMAKE_SOURCE_DIR}/samseg/samseg_cli_wrapper ${CMAKE_INSTALL_PREFIX}/bin/${SCRIPT} @ONLY)"
    )
    endif()  
  
    # remove pip-installed cli wrappers
    remove_files(${FSPYTHON_INSTALL_PREFIX}/python/packages/bin/${SCRIPT})
  endforeach()
endfunction()

install_configured(samseg samseg-long DESTINATION bin)
if(NOT INTEGRATE_SAMSEG)
   if (FSPYTHON_INSTALL_TREE)
      install_pyscript_fspython_tree(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
   else()
      install_pyscript(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
   endif()
else()
   install_samseg_cli()   
endif()
add_test_script(NAME samseg_test SCRIPT test.sh)

