# Simlation of Multi-Wire Proportional Chamber (MWPC) as part of 2024 CERN Summer Internship
This repository holds .C files with CMake Makefiles to build an application which simulates and visualises the signal outputs of MWPC, a gaseous particle detector developed at CERN.
Four configurations are presented with different number of layers and cathode architectures (foil vs. wires).
<img width="887" alt="MWPC_configs" src="https://github.com/user-attachments/assets/6e4eaeac-ddbe-437e-a096-00877f9af15a">
![oldMWPC](https://github.com/user-attachments/assets/52df8412-05eb-4197-bf34-0b05e8d90d80)
[Config1SignalX6_190724.pdf](https://github.com/user-attachments/files/17704464/Config1SignalX6_190724.pdf)

The program runs of Garfield++, which requires C++, CMake, GNU tools, and ROOT. OpenMP is optional for accelerating the simulation (remove OpenMP requirements from CMakeLists.txt if you wish to not use it).
You can either install Garfiled++ or build it from source by following this Garfield++ [User Guide](https://garfieldpp.web.cern.ch/documentation/), "Chapter 2 - Getting Started".
