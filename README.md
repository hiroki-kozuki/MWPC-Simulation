# Simlation of Multi-Wire Proportional Chamber (MWPC) as part of 2024 CERN Summer Internship
This repository holds .C files with CMake Makefiles to build an application which simulates and visualises the signal outputs of MWPC, a gaseous particle detector developed at CERN.
Four configurations are presented with different number of layers and cathode architectures (foil vs. wires).

The program runs of Garfield++, which requires C++, CMake, GNU tools, and ROOT. OpenMP is optional for accelerating the simulation (remove OpenMP requirements from CMakeLists.txt if you wish to not use it).
You can either install Garfiled++ or build it from source following this Garfield++ [User Guide](https://garfieldpp.web.cern.ch/documentation/).
