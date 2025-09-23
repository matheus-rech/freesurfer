# install all licenses to the doc directory
file(GLOB LICENSES "license*")
install(FILES ${LICENSES} DESTINATION docs)
