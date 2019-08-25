aux_source_directory(. zlib_sources)

add_library (zlib SHARED ${zlib_sources})
target_link_libraries (zlib ${zlib_dependencies})

install (TARGETS zlib LIBRARY DESTINATION lib)
