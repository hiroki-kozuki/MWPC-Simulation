# Simlation of Multi-Wire Proportional Chamber (MWPC) as part of 2024 CERN Summer Internship
This repository holds .C files with CMake Makefiles to build an application which simulates and visualises the signal outputs of MWPC, a gaseous particle detector developed at CERN.
Four configurations are presented with different numbers of layers and cathode architectures (foils vs. wires).

The program runs on Garfield++, which requires C++, CMake, GNU tools, and ROOT. OpenMP is optional for accelerating the simulation (remove OpenMP requirement from CMakeLists.txt if you wish to not use it).
You can either install Garfiled++ or build it from source by following this Garfield++ [User Guide](https://garfieldpp.web.cern.ch/documentation/), "Chapter 2 - Getting Started" or from [here](https://garfieldpp.web.cern.ch/garfieldpp/getting-started/).

More information about Garfield++, dependencies, and build system can be found below:
- Garfield++: [Home](https://garfieldpp.web.cern.ch/garfieldpp/), [GitLab](https://gitlab.cern.ch/garfield/garfieldpp)
- CMake: [Basic CMake](https://root.cern/install/basic_cmake/), [Documentation](https://cmake.org/documentation/)
- ROOT: [Dependencies](https://root.cern/install/dependencies/), [Installing ROOT](https://root.cern/install/), [Building ROOT from source](https://root.cern/install/build_from_source/), [Documentation](https://root.cern/doc/master/)
- GNU tools: [Documentation](https://www.gnu.org/manual/manual.en.html)
- OpenMP: [Documentation](https://www.openmp.org/resources/refguides/)


<img width="800" alt="MWPC_configs" src="https://github.com/user-attachments/assets/6e4eaeac-ddbe-437e-a096-00877f9af15a">
<img width="400" alt="config1Geometry" src="https://github.com/user-attachments/assets/2f20ddc7-ea21-4106-a8b3-b5d793cd867b">
<img width="400" alt="Screenshot 2024-11-11 at 15 43 11" src="https://github.com/user-attachments/assets/48af6040-dc4e-44f1-8547-4a0f4d506210">
<img width="800" alt="Screenshot 2024-11-11 at 15 37 57" src="https://github.com/user-attachments/assets/74b902b8-6820-48b8-888e-1e19c5d5bba3">

On the last plot, orange represents the induced current in one of the two anode sensing wires due to electron drifts, red represents the induced current due to ion drifts, and blue represents the total induced current (electrons + ions). The program (Heed) treats the ionisation of the gaseous medium in clusters resulting in numerous sharp peaks corresponding to different electron and ion clusters. 
