project(mri_3d_photo_recon)

set(photo_recon_dest "python/packages/3d_photo_recon")

# install the python package and dependencies
if (FSPYTHON_INSTALL_TREE)
    install_directories_fspython_tree(ext mni photo_reconstruction scripts DESTINATION ${photo_recon_dest})
else()
    install_directories(ext mni photo_reconstruction scripts DESTINATION ${photo_recon_dest})
endif()

# install the entry point
install(PROGRAMS mri_3d_photo_recon DESTINATION bin)

