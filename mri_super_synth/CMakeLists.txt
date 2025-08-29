project(mri_super_synth)

install_configured(mri_super_synth DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
   install(DIRECTORY SuperSynth DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python/packages)
else()
   install(DIRECTORY SuperSynth DESTINATION python/packages)
endif()

