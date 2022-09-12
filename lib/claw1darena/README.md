# claw1dArena

``claw1dArena`` is a C++ library for the numerical integration of 1d conservation and balance laws.

It's main design goal is to allow experimenting with different combinations of reconstruction procedures, numerical fluxes, time integrators and so on.

``claw1dArena`` provides different libraries for each task (grid management, recontructions, numerical fluxes, time integration, sundry tools), each of which can be easily extended by adding your own favourite methods.

The chosen balance between efficiency and convenience is the following: a different executable is built for each conservation law (so that the compiler has some room for optimization) but all choices regarding the numerical method can be made at runtime.

The ``src/claw1d.cpp`` program is an example on how the library may be used to build a simulation code. To discover all the available options, please browse through the source file ``claw1d.cpp`` and ``utils/parseOptions``.

Of course you are very welcome to contribute your own code to any of the libraries so that your methods may be available for everyone to test.
