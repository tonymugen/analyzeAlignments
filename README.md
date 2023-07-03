# Overview

A C++14 library and software to extract unique sequences from FASTA alignments. These unique sequences can be found by scanning the entire alignment, providing a position and window length, or providing a query sequence.

# Dependencies

Building the library and binaries requires a C++14 compiler. The build process requires `cmake` version 3.11 or later. Sequence query uses Smith-Watemrman alignment, implemented [here](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library). This code is included as a submodule.

# Download and install

The repository comes with a submodule, so to clone use

```sh
git clone --recurse-submodules https://github.com/tonymugen/analyzeAlignments
```

Next, create a build directory

```sh
cd vash
mkdir build
```

Finally, run `cmake` to build and install the software

```sh
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
```

Installation may require root privileges.

# Tests

Optionally, one can also build the unit tests. These require [Catch2](https://github.com/catchorg/Catch2), although its installation is taken care of by `cmake`. To build the tests, create a `build-Tests` directory, say, and run

```sh
cd build
cmake -DCMAKE_BUILD_TYPE=Test -DBUILD_TESTS=ON ..
cmake --build .
```

To run the tests from the build directory, simply run

```sh
./tests
```

# Binaries

Two binaries are built as part of the project. Command line flags and their descriptions can be printed by running the programs without parameters.

## homoruns

The `homoruns` binary takes an alignment and sliding window parameters (window and step size) and outputs unique sequence counts for each window. Sequences themselves are not saved, but counts are reported for each unique sequence.

## extractWindow

The `extractWindow` binary takes an alignment and either a start window position and length or a query sequence. It returns all unique sequences in the window (or best matches to the query) with their counts. The sequences can be optionally sorted by their counts in descending order.
