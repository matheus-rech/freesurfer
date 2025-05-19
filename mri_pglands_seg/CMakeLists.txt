project(mri_pglands_seg)
#add_subdirectory(pglands_seg)
# install the bash wrapper for the main entry point,
if (FSPYTHON_INSTALL_TREE)
    set(pglands_dest "\.\./fspython/models/pglands_seg")
    install_pyscript_fspython_tree(mri_pglands_seg)
else()
    set(pglands_dest "models/pglands_seg")
    install_pyscript(mri_pglands_seg)
endif()

# install the cLUT, template and model
install_symlinks(TYPE FILES pglands_seg/mni152_label_template.mgz pglands_seg/pglands.ctab pglands_seg/pglands_seg.pth DESTINATION ${pglands_dest})

# unpack and install the qa label volumes
install_tarball(pglands_seg_qa_data.tar.gz DESTINATION ${pglands_dest})
