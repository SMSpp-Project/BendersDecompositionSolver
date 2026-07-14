# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ BendersDecompositionSolver                      #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)
smspp_derive_version(SMSPP_MODULE_VERSION)

project(BendersDecompositionSolver
        VERSION ${SMSPP_MODULE_VERSION}
        DESCRIPTION "SMS++ BendersDecompositionSolver module"
        HOMEPAGE_URL https://gitlab.com/smspp/bendersdecompositionsolver
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})

# C-safe form of modName for the preprocessor macros in <modName>Config.h
string(MAKE_C_IDENTIFIER "${modName}" modMacro)
set(modNamespace "SMS++")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# ----- Requirements -------------------------------------------------------- #
# If it's not being called by the umbrella, we need to
# look for the system-installed libraries.
# The master and the subproblem Solver are instantiated through the Solver
# factory at run time, so they are not link-time dependencies here: only the
# core SMS++ (which provides BendersBFunction) is required.
if (NOT hasParent)
    find_package(SMS++ REQUIRED)
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_17)

if (MSVC AND BUILD_SHARED_LIBS)
    target_link_options(${modName} INTERFACE "/INCLUDE:SMSpp_force_load_${modName}")
endif ()

# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
        src/BendersDecompositionSolver.cpp)

target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

target_link_libraries(${modName} PUBLIC
        ${modNamespace}::SMS++)

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- Tests --------------------------------------------------------------- #
# The validation test compares BendersDecompositionSolver against a monolithic
# reference solve, so it needs a master Solver (BundleSolver) and a subproblem
# Solver (a *MILPSolver); it is built only when both modules are available.
if (BUILD_BundleSolver AND BUILD_MILPSolver)
    # enable testing in this directory's scope: the umbrella calls
    # enable_testing() only after this subdirectory has been processed, so the
    # test would otherwise not be registered with CTest
    enable_testing()
    add_executable(${modName}_test test/test.cpp)
    target_link_libraries(${modName}_test PRIVATE
            ${modNamespace}::${modName}
            ${modNamespace}::BundleSolver
            ${modNamespace}::MILPSolver)
    add_test(NAME ${modName}_test
            COMMAND ${modName}_test
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test)

    # extirpation comparison: BDS+BundleSolver vs the ad-hoc Benders (BenForm)
    # of CapacitatedFacilityLocationBlock on a real CFL instance, checking the
    # solutions coincide (same model) and reporting the times
    if (BUILD_CapacitatedFacilityLocationBlock)
        add_executable(${modName}_cfl_test test/test_cfl.cpp)
        target_link_libraries(${modName}_cfl_test PRIVATE
                ${modNamespace}::${modName}
                ${modNamespace}::BundleSolver
                ${modNamespace}::MILPSolver
                ${modNamespace}::CapacitatedFacilityLocationBlock)
        add_test(NAME ${modName}_cfl_test
                COMMAND ${modName}_cfl_test
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test)
    endif ()

    # unit test for PolyhedralFunction::remove_parallel_rows() and
    # PolyhedralFunctionBlock::remove_redundant_rows() (the row pruning that
    # replaced the standalone CutProcessing class)
    add_executable(${modName}_prune_test test/test_prune.cpp)
    target_link_libraries(${modName}_prune_test PRIVATE
            ${modNamespace}::SMS++
            ${modNamespace}::MILPSolver)
    add_test(NAME ${modName}_prune_test
            COMMAND ${modName}_prune_test
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test)
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .dll
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # .so .dylib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}) # .lib

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/BendersDecompositionSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# ----- Add the build tree to RPATH------------------------------------------ #
# Add the binary directory to RPATH, so executables that use
# dynamic loading will look into it.
if (hasParent)
    set(SMSPP_BUILD_RPATH
            "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
            PARENT_SCOPE)
endif ()

# ----- Packaging support --------------------------------------------------- #
# The following commands allow one to build a proper package using CPack.
# See: https://cmake.org/cmake/help/latest/manual/cpack.1.html

# Using this with the umbrella project is too messy.
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "SMS++ Team")
    set(CPACK_PACKAGE_CONTACT "Donato Meoli")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
    set(CPACK_SOURCE_IGNORE_FILES
            \\\\.git/
            \\\\.gitignore$
            \\\\.gitattributes$)
    include(CPack)
endif ()

# --------------------------------------------------------------------------- #
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #
