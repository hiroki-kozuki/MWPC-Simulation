# Simlation of Multi-Wire Proportional Chamber (MWPC) as part of 2024 CERN Summer Internship
This repository holds .C files with CMake Makefiles to build an application which simulates and visualises the signal outputs of MWPC, a gaseous particle detector developed at CERN.
Four configurations are presented with different number of layers and cathode architectures (foil vs. wires).
<img width="887" alt="MWPC_configs" src="https://github.com/user-attachments/assets/6e4eaeac-ddbe-437e-a096-00877f9af15a">
![config1Geometry](https://github.com/user-attachments/assets/2f20ddc7-ea21-4106-a8b3-b5d793cd867b)
<img width="1618" alt="Screenshot 2024-11-11 at 15 37 57" src="https://github.com/user-attachments/assets/74b902b8-6820-48b8-888e-1e19c5d5bba3">

The program runs of Garfield++, which requires C++, CMake, GNU tools, and ROOT. OpenMP is optional for accelerating the simulation (remove OpenMP requirements from CMakeLists.txt if you wish to not use it).
You can either install Garfiled++ or build it from source by following this Garfield++ [User Guide](https://garfieldpp.web.cern.ch/documentation/), "Chapter 2 - Getting Started".
