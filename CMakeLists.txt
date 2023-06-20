cmake_minimum_required(VERSION 3.11)
include(CheckCXXCompilerFlag)

project(analizeAlignments 
	VERSION 0.1 
	LANGUAGES CXX C
	)

set(CMAKE_CXX_EXTENSIONS OFF)
if (NOT DEFINED CMAKE_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD 14)
endif()

# Adding a Profile build type
get_property(isMultiConfig GLOBAL
	PROPERTY GENERATOR_IS_MULTICONFIG
)
if(isMultiConfig)
	if(NOT "Test" IN_LIST CMAKE_CONFIGURATION_TYPES)
		list(APPEND CMAKE_CONFIGURATION_TYPES Test)
	endif()
else()
	set(allowedBuildTypes Debug Release Test)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
		STRINGS "${allowedBuildTypes}"
	)
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE) # default build type is Release
	elseif(NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
		message(FATAL_ERROR
			"Unknown build type ${CMAKE_BUILD_TYPE}"
		)
	endif()
endif()

#
# Set warnings
#
set(CLANG_WARNINGS
	-Wall
	-Wextra
	-Wconversion
	-Wpedantic
	-Wshadow
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Wcast-align
	-Woverloaded-virtual
	-Wsign-conversion
	-Wnull-dereference
	-Wdouble-promotion
	-Wformat=2
	-Wimplicit-fallthrough
)
set(GCC_WARNINGS
	${CLANG_WARNINGS} 
	-Wmisleading-indentation
	-Wduplicated-cond
	-Wduplicated-branches
	-Wlogical-op
	-Wuseless-cast
)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
	set(PROJECT_WARNINGS_CXX ${CLANG_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(PROJECT_WARNINGS_CXX ${GCC_WARNINGS})
else()
	message(AUTHOR_WARNING "No compiler warnings set for CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
#
# Tests with Catch2
#
option(BUILD_TESTS
	"Build tests"
	OFF
)
if(BUILD_TESTS)
	Include(FetchContent)

	FetchContent_Declare(
	  Catch2
	  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	  GIT_TAG        v3.3.2 # or a later release
	)
	FetchContent_MakeAvailable(Catch2)
	list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
#
# Find active available sanitizers
#
set(SANITIZER_LIST "")
if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(SANITIZERS "")
	check_cxx_compiler_flag(-fsanitize=address   HAS_ADDRESS_SANITIZER)
	check_cxx_compiler_flag(-fsanitize=leak      HAS_LEAK_SANITIZER)
	check_cxx_compiler_flag(-fsanitize=undefined HAS_UNDEFINED_SANITIZER)
	if(${HAS_ADDRESS_SANITIZER})
		list(APPEND SANITIZERS "address")
	endif()
	if(${HAS_LEAK_SANITIZER})
		list(APPEND SANITIZERS "leak")
	endif()
	if(${HAS_UNDEFINED_SANITIZER})
		list(APPEND SANITIZERS "undefined")
	endif()
	list(JOIN
		SANITIZERS
		","
		SANITIZER_LIST)
endif()

set(CMAKE_CXX_FLAGS_PROFILE
	"-g -p -O3"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_PROFILE
	"-g -p -O3"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_PROFILE
	"-g -p -O3"
)

set(CMAKE_CXX_FLAGS_DEBUG
	"-g -O1 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_DEBUG
	"-g -O1 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_DEBUG
	"-g -O1 -fsanitize=${SANITIZER_LIST}"
)

set(CMAKE_CXX_FLAGS_TEST
	"-g -O3 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_TEST
	"-g -O3 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_TEST
	"-g -O3 -fsanitize=${SANITIZER_LIST}"
)

include(GNUInstallDirs)
# executables
add_executable(homoruns
	apps/homoruns.cpp
)
target_include_directories(homoruns
	PRIVATE include
)
target_link_libraries(homoruns
	PRIVATE analizeAlignments
)
target_compile_options(homoruns
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
set_target_properties(homoruns PROPERTIES
	CXX_STANDARD_REQURED ON
)
install(TARGETS homoruns
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_executable(extractWindow
	apps/extractWindow.cpp
)
target_include_directories(extractWindow
	PRIVATE include
)
target_link_libraries(extractWindow
	PRIVATE analizeAlignments
)
target_compile_options(extractWindow
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
set_target_properties(extractWindow PROPERTIES
	CXX_STANDARD_REQURED ON
)
install(TARGETS extractWindow
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# library
add_library(smithWaterman
	externals/stripedSW/src/ssw.c
	externals/stripedSW/src/ssw_cpp.cpp
)
target_include_directories(smithWaterman
	PRIVATE externals/stripedSW/src
)
set_target_properties(smithWaterman PROPERTIES
	POSITION_INDEPENDENT_CODE ON
)
target_compile_options(smithWaterman
	PRIVATE -msse2
)
set_target_properties(smithWaterman PROPERTIES
	CXX_STANDARD_REQURED ON
)

add_library(analizeAlignments
	src/fastaParser.cpp
	src/extraFunctions.cpp
)
target_include_directories(analizeAlignments
	PRIVATE include
	PRIVATE externals/stripedSW/src
)
target_link_libraries(analizeAlignments
	PRIVATE smithWaterman
)
set_target_properties(analizeAlignments PROPERTIES
	POSITION_INDEPENDENT_CODE ON
)
target_compile_options(analizeAlignments
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
set_target_properties(analizeAlignments PROPERTIES
	CXX_STANDARD_REQURED ON
)
# no install for now

# tests
if(BUILD_TESTS)
	add_executable(tests
		tests/tests.cpp
	)
	target_link_libraries(tests
		PRIVATE analizeAlignments
		PRIVATE Catch2::Catch2WithMain
	)
	target_include_directories(tests
		PRIVATE include
	)
	include(CTest)
	include(Catch)
	catch_discover_tests(tests)
endif()
