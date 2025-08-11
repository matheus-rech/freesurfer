project(qatools)

if(NOT INSTALL_PYTHON_DEPENDENCIES)
  set(PIP_FLAGS "--no-dependencies")
endif()

if (FSPYTHON_INSTALL_TREE)
   set(qatools_dest "${FSPYTHON_INSTALL_PREFIX}/python/packages")
   install_pyscript_fspython_tree(qatools.py)
else()
   set(qatools_dest "${CMAKE_INSTALL_PREFIX}/python/packages")
   install_pyscript(qatools.py)
endif()

# run 'pip install' and direct package installs to freesufer/python/packages
# requirements are currently defined by python/requirements-extra.txt

if(HOST_OS MATCHES "Ubuntu24")
   if (NOT PYTHON_EXECUTABLE)
      set(PYTHON_EXECUTABLE "${FS_PACKAGES_DIR}/fspython/3.8/bin/python")
      message(STATUS "For HOST_OS=${HOST_OS} set PYTHON_EXEUTABLE to ${PYTHON_EXECUTABLE}")
   endif()
endif()

install(CODE "
  message(STATUS \"Installing QA tools\")
  execute_process(COMMAND bash -c \"${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-dependencies --upgrade --target=${qatools_dest}\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    # message(FATAL_ERROR \"Could not install qatools python module\")
    message(SEND_ERROR \"Could not install qatools python module\")
  endif()"
)

