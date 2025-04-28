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
    # config fspython wrapper script
    if (FSPYTHON_INSTALL_TREE)
      config_fspythonwrapper_fspython_tree(${SCRIPT})
    else()
      config_fspythonwrapper(${SCRIPT})
    endif()
     
    # config fspython wrapper script to wrap around ${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py
    if(NOT ${USER_SAMSEG_PATH} STREQUAL "")
      # create links in ${FSPYTHON_INSTALL_PREFIX}/python/scripts/ to real python scripts in user samseg package
      symlink(${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py ${FSPYTHON_INSTALL_PREFIX}/python/scripts/${SCRIPT})      
    else()
      install(CODE "
        # copy samseg cli scripts from installed samseg package to ${FSPYTHON_INSTALL_PREFIX}/python/scripts/
        execute_process(COMMAND bash -c \"mkdir -p ${FSPYTHON_INSTALL_PREFIX}/python/scripts && cp ${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py ${FSPYTHON_INSTALL_PREFIX}/python/scripts/${SCRIPT}\" RESULT_VARIABLE retcode)
        if(NOT \${retcode} STREQUAL 0)
          message(FATAL_ERROR \"Could not copy pip-installed samseg cli\")
        endif()"
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

