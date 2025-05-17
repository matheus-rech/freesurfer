project(mri_histo_atlas_segment)

install_configured(mri_histo_atlas_segment DESTINATION bin)
install_configured(mri_histo_atlas_segment_fast DESTINATION bin)
install_configured(mri_histo_atlas_segment_fireants DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
   install(DIRECTORY ERC_bayesian_segmentation DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python/packages)
else()
   install(DIRECTORY ERC_bayesian_segmentation DESTINATION python/packages)
endif()

