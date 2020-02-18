# Ashperix-Palabos coupling

## Preface

This is an illustrative example of how to couple a CFD code to the Aspherix® 
simulation engine. The physics depicted in this example are not correct, however 
the type and order of the instructions are.

## Prerequisites

To compile and run this example, two more git repositories need to be cloned.

### Aspherix-CoSim-Socket-Lib

The CoSim socket lib handles the communication with Aspherix®. To install it, clone 
[this repository](https://github.com/CFDEMproject/Aspherix-CoSim-Socket-Lib). Then, 
compile and install the content with 

    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=path/to/install/destination path/to/socketlib/src
    make
    make install

### Palabos

Palabos, the **Pa**rallel **La**ttice **Bo**ltzmann **S**olver 
([link](https://palabos.unige.ch/)) is a lattice Boltzmann 
library developed initially by FlowKit Ltd., and now maintained by the University of 
Geneva. Its development is ongoing, and the current version can be obtained from UNIGE's 
[gitlab repository](https://gitlab.com/unigespc/palabos).
    
We have created a fork of UNIGE's repository, which preserves the state of Palabos with 
which the content of this repository was developed. The repository can be found 
[here](https://github.com/CFDEMproject/Palabos-fork). To compile and run the Aspherix-Palabos
coupling, clone our fork instead of the official Palabos repository.

No further steps are necessary, since the case compiles the necessary Palabos components on-demand.

## Compilation

Make sure that you compiled the Aspherix CoSim library correctly, and that the Palabos-fork
repository is properly cloned. To compile the case, a few changes to the `Makefile` in this
repository are necessary. The `Makefile` contains several variables to configure the build 
process. The following variables need to be adjusted:

* `palabosRoot`: set this to the full absolute or relative path of your Palabos download.
* `libraryPath`: set this to the full, absolute path where the compiled `libaspherix_cosim_socket.so` is located (usually `cosim/install/dir/lib/`)
* `includePath`: set this to the full, absolute path where the file `aspherix_cosim_socket.h` is located (usually `cosim/install/dir/include/aspherix_cosim_socket/`)

Once these settings are complete, compile the case with

    make

## Execution


To run, the `LD_LIBRARY_PATH` needs to be adjusted to run the case:

    export LD_LIBRARY_PATH=path/to/cosim/lib:$LD_LIBRARY_PATH
    
where the path is the same as `libraryPath` in the `Makefile`. Once this is done, just run

    ./Allrun.sh
    
in the case directory.

