project(mri_pglands_seg)
# install the bash wrapper for the main entry point,
if (FSPYTHON_INSTALL_TREE)
    set(pglands_dest "\.\./fspython/models")
    install_pyscript_fspython_tree(mri_pglands_seg)
else()
    set(pglands_dest "models")
    install_pyscript(mri_pglands_seg)
endif()

# install the cLUT under luts/
install(DIRECTORY pglands_seg DESTINATION ${pglands_dest})

# unpack and install the qa label volumes
install_tarball(pglands_seg_qa_data.tar.gz DESTINATION ${pglands_dest}/pglands_seg)
