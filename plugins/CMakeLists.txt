# --- Sources

set(ufofilter_SRCS
    ufo-fftconvolution-task.c
    ufo-hessian-kernel-task.c
    ufo-hessian-analysis-task.c
    ufo-ring-pattern-task.c
    ufo-ring-writer-task.c
    ufo-local-maxima-task.c
    ufo-blob-test-task.c
    ufo-label-cluster-task.c
    ufo-combine-test-task.c
    ufo-sum-task.c
    ufo-log-task.c
    ufo-brightness-cut-task.c
    ufo-complex-mult-task.c
    ufo-bilateral-task.c
    ufo-piv-contrast-task.c
    ufo-hough-likelihood-task.c
    ufo-hough-likelihood-image-task.c
    ufo-candidate-filter-task.c
    ufo-multi-search-task.c
    ufo-device-info-task.c
    ufo-azimuthal-test-task.c
    ufo-contrast-task.c
    ufo-ordfilt-task.c
    ufo-filter-particle-task.c
    ufo-of-ring-pattern-task.c
    ufo-get-dup-circ-task.c
    ufo-concatenate-result-task.c
    )

file(GLOB ufofilter_KERNELS "kernels/*.cl")


# --- Dependencies

set(UFO_CORE_MINIMUM "0.6")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(PkgConfigVars)

find_package(PkgConfig REQUIRED)
find_package(OpenCL REQUIRED)
find_package(GSL REQUIRED)

pkg_check_modules(UFO ufo>=${UFO_CORE_MINIMUM} REQUIRED)
pkg_check_modules(OCLFFT oclfft REQUIRED)

pkg_check_variable(ufo plugindir)
pkg_check_variable(ufo kerneldir)

link_directories(
    ${UFO_LIBRARY_DIRS}
    ${OCLFFT_LIBRARY_DIRS})

include_directories(
    ${UFO_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${OCLFFT_INCLUDE_DIRS}
    ${OPENCL_INCLUDE_DIRS})


# --- Variables

if (CMAKE_COMPILER_IS_GNUCC OR ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))
    add_definitions("-std=c99 -pedantic -Wall -Wextra -fPIC")
    add_definitions("-Wcast-align -Wcast-qual -Winline -Wmissing-declarations "
                    "-Wmissing-prototypes -Wnested-externs -Wno-long-long "
                    "-Wno-missing-field-initializers -Wpointer-arith "
                    "-Wredundant-decls -Wshadow -Wstrict-prototypes -Wwrite-strings")
endif()

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")


# --- Targets

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

foreach(_src ${ufofilter_SRCS})
    # find plugin suffix
    string(REGEX REPLACE "ufo-([^ \\.]+)-task.*" "\\1" task "${_src}")
    string(REPLACE "-" "" _targetname ${task})
    set(target "ufofilter${_targetname}")

    # build single shared library per filter
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        add_library(${target} MODULE ufo-priv.c ${_src})
    else()
        add_library(${target} SHARED ufo-priv.c ${_src})
    endif()

    target_link_libraries(${target}
        ${UFO_LIBRARIES}
        oclfft)

    install(TARGETS ${target}
        LIBRARY DESTINATION ${UFO_PLUGINDIR})
endforeach()

foreach(_kernel ${ufofilter_KERNELS})
    install(FILES ${_kernel} DESTINATION ${UFO_KERNELDIR})
endforeach()
