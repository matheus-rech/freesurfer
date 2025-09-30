project(mri_super_synth)

install_configured(mri_super_synth DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
   install(DIRECTORY SuperSynth DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python/packages)
   #install(FILES SuperSynth/scripts/inference.py DESTINATION ${FSPYTHON_INSTALL_PREFIX}/ERC_bayesian_segmentation/scripts)
else()
   install(DIRECTORY SuperSynth DESTINATION python/packages)
   #install(FILES SuperSynth/scripts/inference.py DESTINATION python/packages/ERC_bayesian_segmentation/scripts)
endif()

