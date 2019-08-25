set(zstd_sources
./decompress/zstd_decompress.c
./decompress/huf_decompress.c
./compress/zstd_lazy.c
./compress/zstdmt_compress.c
./compress/zstd_double_fast.c
./compress/zstd_fast.c
./compress/hist.c
./compress/fse_compress.c
./compress/zstd_opt.c
./compress/zstd_compress.c
./compress/huf_compress.c
./compress/zstd_ldm.c
./common/xxhash.c
./common/fse_decompress.c
./common/pool.c
./common/zstd_common.c
./common/error_private.c
./common/debug.c
./common/threading.c
./common/entropy_common.c
./deprecated/zbuff_compress.c
./deprecated/zbuff_decompress.c
./deprecated/zbuff_common.c
./dictBuilder/zdict.c
./dictBuilder/divsufsort.c
./dictBuilder/cover.c
./legacy/zstd_v03.c
./legacy/zstd_v07.c
./legacy/zstd_v06.c
./legacy/zstd_v01.c
./legacy/zstd_v02.c
./legacy/zstd_v04.c
./legacy/zstd_v05.c
	)

add_library (zstd SHARED ${zstd_sources})
target_link_libraries (zstd ${zstd_dependencies})

install (TARGETS zstd LIBRARY DESTINATION lib)
