# claw1dArena compilation and installation

``claw1dArena`` runs at least on Linux and MacOS. It has not been tested under other OS's.
If you try it out on another OS, we'd be interested in a sucess report and/or specific instructions to make it work, so that we may include it in this document.

## Dependencies
In order to build ``claw1dArena`` you'll need

* ``cmake``
* a C++ compiler
* the BOOST ``program_option`` library (for command line option processing)
* a working blas/lapack installation and C-headers (lapacke)

## Configuring ``claw1dArena``
``claw1dArena`` uses ``cmake`` as building tool.
Please create a folder for the build, e.g. in the ``claw1dArena``'s root folder:

``
    claw1dArena$ mkdir build
``

Enter the `build` folder

``
    claw1dArena$ cd build
``

 In order to configure ``claw1dArena``, run:

``
    claw1dArena/build$ ccmake ../
``

and then press "c" to configure, adjust the settings as needed, and finally press "g" to generate the Makefiles.

In this way, ``claw1dArena`` can be configured for debugging or production.

In particular you can customize your build by choosing Debug or Release (optimized) compilation or the type of floating point numbers that you wish to use.

## Building and Installing
Once ``cmake`` has configured ``claw1dArena``, you can give the command

``
    claw1dArena/build$ make
``

to build it.

## Building Documentation
The code has a rather complete documentation that can be built with

``
    claw1dArena/build$ make doc
``

provided that a working Doxygen installation is present on the build machine.

## Help
For any problem, please contact <a href="http://www.dipmatematica.unito.it/do/docenti.pl/Alias?matteo.semplice">Matteo Semplice</a>.
