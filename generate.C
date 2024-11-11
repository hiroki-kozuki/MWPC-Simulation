#include <iostream>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  const double pressure = 1.11458 * AtmosphericPressure; // [Torr] 
  const double temperature = 298.15;                     // [K]
 
  // Setup the gas.
  MediumMagboltz gas("Ar", 50., "CO2", 50.);
  gas.SetTemperature(temperature);
  gas.SetPressure(pressure);

  // Set the field range to be covered by the gas table. 
  const size_t nE = 20;                                  // # of field points between emin and emax.
  const double emin = 1.;                                // Minimum E field [V/cm]
  const double emax = 1.0e6;                             // Maximum E field [V/cm] 
  // Flag to request logarithmic spacing.
  constexpr bool useLog = true;
  gas.SetFieldGrid(emin, emax, nE, useLog); 
  // Run Magboltz to generate the gas table.
  const int ncoll = 10;                                  // # of collisions in multiples of 10**7.
  gas.GenerateGasTable(ncoll);
  // Save the table. 
  gas.WriteGasFile("ar_50_co2_50_V5.gas");

}
