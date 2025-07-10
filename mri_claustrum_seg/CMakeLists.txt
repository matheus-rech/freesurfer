project(mri_claustrum_seg)

# handle FSPYTHON tree split install
if (FSPYTHON_INSTALL_TREE)
    set(claustrum_python "\.\./fspython")
    set(claustrum_models "\.\./fspython/models")
else()
    set(claustrum_python "python")
    set(claustrum_models "models")
endif()

# bash wrapper always installed the same way, fspython agnostic
install(PROGRAMS mri_claustrum_seg DESTINATION bin)

# handle the main python entry point
install(PROGRAMS model/claustrum-seg.py DESTINATION "${claustrum_python}/scripts")

# install the model file
install(FILES model/claustrum_seg_20250616.h5 DESTINATION "${claustrum_models}")

# isntall the atlas dir
install_tarball(claustrum_seg.tar.gz DESTINATION average)

