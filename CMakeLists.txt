cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(neBEM)

#--- Find Garfield package
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

#--- Enable parallel processing with OpenMP
find_package(OpenMP REQUIRED)

#--- Build executables
add_executable(generate generate.C)
target_link_libraries(generate Garfield::Garfield)
target_link_libraries(generate OpenMP::OpenMP_CXX)

add_executable(oldMWPC oldMWPC.C)
target_link_libraries(oldMWPC Garfield::Garfield)
target_link_libraries(oldMWPC OpenMP::OpenMP_CXX)

add_executable(oldMWPC2 oldMWPC2.C)
target_link_libraries(oldMWPC2 Garfield::Garfield)
target_link_libraries(oldMWPC2 OpenMP::OpenMP_CXX)

add_executable(oldMWPCManyParticles oldMWPCManyParticles.C)
target_link_libraries(oldMWPCManyParticles Garfield::Garfield)
target_link_libraries(oldMWPCManyParticles OpenMP::OpenMP_CXX)

add_executable(newMWPC newMWPC.C)
target_link_libraries(newMWPC Garfield::Garfield)
target_link_libraries(newMWPC OpenMP::OpenMP_CXX)

add_executable(config1OldFoil5LayerMWPC config1OldFoil5LayerMWPC.C)
target_link_libraries(config1OldFoil5LayerMWPC Garfield::Garfield)
target_link_libraries(config1OldFoil5LayerMWPC OpenMP::OpenMP_CXX)

add_executable(config2DWCWire6LayerMWPC config2DWCWire6LayerMWPC.C)
target_link_libraries(config2DWCWire6LayerMWPC Garfield::Garfield)
target_link_libraries(config2DWCWire6LayerMWPC OpenMP::OpenMP_CXX)

add_executable(config3NewWire5LayerMWPC config3NewWire5LayerMWPC.C)
target_link_libraries(config3NewWire5LayerMWPC Garfield::Garfield)
target_link_libraries(config3NewWire5LayerMWPC OpenMP::OpenMP_CXX)

add_executable(config4Foil6LayerMWPC config4Foil6LayerMWPC.C)
target_link_libraries(config4Foil6LayerMWPC Garfield::Garfield)
target_link_libraries(config4Foil6LayerMWPC OpenMP::OpenMP_CXX)

