#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

#include "Garfield/SolidBox.hh" 
#include "Garfield/SolidWire.hh" 
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Make a gas medium.
  MediumMagboltz gas;
  gas.LoadGasFile("ar_50_co2_50_V5.gas");
  auto installdir = std::getenv("GARFIELD_INSTALL");
  if (!installdir) {
    std::cerr << "GARFIELD_INSTALL variable not set.\n";
    return 1;
  }
  const std::string path = installdir;
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt"); // Also try "/share/Garfield/Data/IonMobility_CO2+_CO2.txt"
  
  //Plot drift velocity vs. E field of ions:
  ViewMedium mediumView;
  mediumView.SetMedium(&gas);
  mediumView.PlotElectronVelocity('e');

  // Define conductor material. Change to tungsten.
  MediumConductor W;

  // Geometry.
  GeometrySimple geo;

  // Coordinates of centre of gravity of film 1 [cm]:
  const double cog_x1 = 0., cog_y1 = 0., cog_z1 = 0.; // 5mm appart from anode plane.
  // Coordinates of centre of gravity of film 2 [cm]:
  const double cog_x2 = 0., cog_y2 = 0., cog_z2 = 1.;  // 5mm appart from anode plane.
  // Coordinates of centre of gravity of film 3 [cm]:
  const double cog_x3 = 0., cog_y3 = 0., cog_z3 = 2.;  // 5mm appart from anode plane.
  // Half lengths of film 1 [cm]:
  const double hl_x1 = 5., hl_y1 = 5., hl_z1 = 0.00125; // 25 microns thick
  // Half lengths of film 2 [cm]:
  const double hl_x2 = 5., hl_y2 = 5., hl_z2 = 0.00125; // 25 microns thick
  // Half lengths of film 3 [cm]:
  const double hl_x3 = 5., hl_y3 = 5., hl_z3 = 0.00125; // 25 microns thick
  SolidBox cathode1(cog_x1, cog_y1, cog_z1, hl_x1, hl_y1, hl_z1);
  cathode1.SetBoundaryPotential(-4000.);
  // cathode1.SetLabel(const std::string& label);
  SolidBox cathode2(cog_x2, cog_y2, cog_z2, hl_x2, hl_y2, hl_z2);
  cathode2.SetBoundaryPotential(-4000.);
  // cathode2.SetLabel(const std::string& label);
  SolidBox cathode3(cog_x3, cog_y3, cog_z3, hl_x3, hl_y3, hl_z3);
  cathode3.SetBoundaryPotential(-4000.);
  // cathode3.SetLabel(const std::string& label);
  geo.AddSolid(&cathode1, &W);
  geo.AddSolid(&cathode2, &W);
  geo.AddSolid(&cathode3, &W);

  // Dimensions of single anode wire [cm]:
  const double radius = 0.001;                   // 20 microns thick.
  const double halflength = 5.;

  // Generate 100 anode wires for x position.
  std::vector<SolidWire*> anodeWiresX(100, nullptr);
  for (int i = 0; i < 100; ++i) {                // 100 wires.
    int j = i - 50;                              // Offset by 50 wires to correctly centre the anode plane.
    double xpos = 0.05 + (double)j*0.1;          // 1 mm spacing along x axis. Offset by 0.5 mm to centre the anode plane.
    anodeWiresX[i] = new SolidWire(xpos, 0.0, 0.5, radius, halflength, 0, 1, 0); // Last 3 inputs define the axis of wires.
    anodeWiresX[i]->SetBoundaryPotential(0.);          // Grounded anode plane.
    anodeWiresX[i]->SetLabel("anodeX");
    geo.AddSolid(anodeWiresX[i], &W);
  }
  // Similarly, generate 100 anode wires for y position.
  std::vector<SolidWire*> anodeWiresY(100, nullptr);
  for (int i = 0; i < 100; ++i) {                
    int j = i - 50;                              
    double ypos = 0.05 + (double)j*0.1;          
    anodeWiresY[i] = new SolidWire(0.0, ypos, 1.5, radius, halflength, 1, 0, 0); 
    anodeWiresY[i]->SetBoundaryPotential(0.);
    anodeWiresY[i]->SetLabel("anodeY");
    geo.AddSolid(anodeWiresY[i], &W);
  }
  // Add boundary volume filled with gas.
  geo.SetMedium(&gas);

  // Plot device geometry in 2D
  //ViewGeometry geomView2d(&geo);
  //geomView2d.SetArea(-10, -10, -1, 10, 10, 2);
  //geomView2d.SetPlane(0, 1, 0, 0, 0, 0.0);
  //geomView2d.Plot2d();

  // Plot device geometry in 3D
  ViewGeometry geomView3d(&geo);
  geomView3d.Plot();

	// ALTERNATIVELY, create a lambda function of the analytical potential and analytical weighting potential (much faster than neBEM!). 
	// Condition the equation on the z position of the electron. This lambda function should be written using ComponentUser and call it as cmp. 
        // Them replace all &nebem with &cmp in the following section of the code.
        // Might need to create two separate sensors (not anodeX and anodeY planes) if using cmp instead of neBEM. 
  
  // Compute fields with neBEM. By default, neBEM calculates both the actual fields 
  // and weighting fields (used for signal calculation).
  ComponentNeBem3d nebem;
  nebem.SetGeometry(&geo);
  nebem.SetTargetElementSize(0.1);              // 0.1 // Linear size of elements measured along their edges.
  nebem.SetMinMaxNumberOfElements(3, 15);        // (3, 15) // Smallest and largest # of elements.
  nebem.SetNumberOfThreads(8);
  // box1.SetDiscretisationLevel();
  // box2.SetDiscretisationLevel();
  // wires.SetDiscretisationLevel();
  // nebem.EnableDebugging();
  nebem.Initialise();
 
  //Medium* medium = nullptr; 
  //double ex = 0., ey = 0., ez = 0., v = 0.;
  //int status = 0;
  //nebem.ElectricField(0, 0, 0, ex, ey, ez, v, medium, status);
  //std::printf("E = (%15.8f, %15.8f %15.8f), V = %15.8f, status = %d\n", ex, ey, ez, v, status);

  
  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&nebem);
  sensor.AddElectrode(&nebem, "anodeX");
  sensor.AddElectrode(&nebem, "anodeY");
  // Set the time bins [ns]
  const unsigned int nTimeBins = 1000;  // 1000 bins.
  const double tmin = 0.;              // From 0 ns.
  const double tmax = 1000.;            // Upto 1000 ns.
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Simulate ionisation with Heed.
  std::cout<<"-- Particle track simulation --"<<std::endl;
  TrackHeed track;
  track.SetParticle("muon");
  track.SetEnergy(170e9);
  track.SetSensor(&sensor);
  // Compute the track.
  track.EnableDebugging();

  // Compute electron drift lines with RKF integration.
  DriftLineRKF drift;
  drift.SetSensor(&sensor);
  drift.UseWeightingPotential(true);
  drift.EnableAvalanche();
  drift.EnableSignalCalculation();
  drift.EnableIonTail();

  // Set up the viewer of the simulated drift lines.
  // Create a new canvas.
  TCanvas* cD = new TCanvas("cD"," ",600,600);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  drift.EnablePlotting(&driftView);
  track.EnablePlotting(&driftView);

  // Set the initial position of particle: starting at z0 = 0.1 mm along perpendicular axis (z axis) of MWPC.
  const double x0 = 0.05, y0 = 0.05, z0 = 0.01, t0 = 0.;
  // Set the track direction and define the track: 
  // Passing near the centre of the MWPC and perpendicularly to its planes.
  const double dx0 = 0., dy0 = 0., dz0 = 1.;
  track.NewTrack(x0,y0,z0,t0,dx0,dy0,dz0);
  
  // Loop over the clusters.
  int loopcounter = 0;
  double enx = 0., eny = 0., enz = 0., ent = 0.;
  int st = 0;  // st is a status code that gets modified when error occurs. 
  for (const auto& cluster : track.GetClusters()) {
  	for (const auto& electron : cluster.electrons) {
		drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);	// Causes the sharp peaks in the signal.
		drift.DriftIon(electron.x, electron.y, electron.z, electron.t);		// Causes the smooth tail in the signal.
		drift.GetEndPoint(enx, eny, enz, ent, st);  
                std::cout<<"x: "<<enx<<" y: "<<eny<<" z: "<<enz<<" t: "<<ent<<" stat: "<<st<<std::endl;
	}
  }
  

  
  std::cout<<"Induced charge: "<<sensor.GetTotalInducedCharge("anode")<<std::endl;

  // Plot the potential and weighting potential contours.
  // Actual potential in x-y plane.
  ViewField potViewZ(&nebem);
  potViewZ.SetArea(-5., -5., 0., 5., 5., 2.);
  potViewZ.SetPlane(0., 0., 1., 0., 0., 0.5); // Plot at (0,0,0.5) with +ve z axis out of the page.
  potViewZ.PlotContour("v");
  // Actual potential in x-z plane.
  ViewField potViewX(&nebem);
  potViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  potViewX.SetPlane(0., -1., 0., 0., 0., 1.); // Plot at (0,0,1.) with +ve y axis into the page.
  potViewX.PlotContour("v");
  // Weighting potential around anodeX.
  //ViewField wpotViewX(&nebem);
  //wpotViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  //wpotViewX.SetPlane(0., -1., 0., 0., 0., 1.); // Plot at (0,0,1.) with +ve y axis into the page.
  //wpotViewX.PlotContourWeightingField("anodeX","v");
  // Actual potential in y-z plane.
  ViewField potViewY(&nebem);
  potViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  potViewY.SetPlane(-1., 0., 0., 0., 0., 1.); // Plot at (0,0,1.) with +ve x axis into the page.
  potViewY.PlotContour("v");
  // Weighting potential around anodeY.
  //ViewField wpotViewY(&nebem);
  //wpotViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  //wpotViewY.SetPlane(-1., 0., 0., 0., 0., 0.); // Plot at (0,0,1.) with +ve x axis into the page.
  //wpotViewY.PlotContourWeightingField("anodeY","v");

  // Plot field and weighting field profiles along z.
  ViewField fieldView(&nebem);
  fieldView.PlotProfile(0., 0., 0., 0., 0., 2.);
  // Actual field in x-y plane.
  ViewField fieldViewZ(&nebem);
  fieldViewZ.SetArea(-5., -5., 0., 5., 5., 2.);
  fieldViewZ.SetPlane(0., 0., 1., 0., 0., 0.5);  // Plot at (0,0,0.5) with +ve z axis out of the page.
  fieldViewZ.PlotContour("e");
  // Actual field in x-z plane.
  ViewField fieldViewX(&nebem);
  fieldViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  fieldViewX.SetPlane(0., -1., 0., 0., 0., 1.);  // Plot at (0,0,1) with +ve y axis into the page.
  fieldViewX.PlotContour("e");
  // Weighting field around anodeX.
  //ViewField wfieldViewX(&nebem);
  //wfieldViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  //wfieldViewX.SetPlane(0., -1., 0., 0., 0., 1.);  // Plot at (0,0,1) with +ve y axis into the page.
  //wfieldViewX.PlotContour("e");
  //wfieldViewX.PlotContourWeightingField("anodeX","e");
  // Actual field in y-z plane.
  ViewField fieldViewY(&nebem);
  fieldViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  fieldViewY.SetPlane(-1., 0., 0., 0., 0., 1.5);  // Plot at (0,0,1.5) with +ve x axis out of the page.
  fieldViewY.PlotContour("e");
  // Weighting field around anodeY.
  //ViewField wfieldViewY(&nebem);
  //wfieldViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  //wfieldViewY.SetPlane(-1., 0., 0., 0., 0., 1.5);  // Plot at (0,0,1.5) with +ve x axis out of the page.
  //wfieldViewY.PlotContourWeightingField("anodeY","e");

  // Plot the induced current in the two anode wire planes.
  ViewSignal signalViewX;            // Anode wires for x coordinate
  signalViewX.SetSensor(&sensor);
  signalViewX.PlotSignal("anodeX");
  ViewSignal signalViewY;            // Anode wires for y coordinate
  signalViewY.SetSensor(&sensor);
  signalViewY.PlotSignal("anodeY");

  //Set the object with which to plot and canvas.
  ViewGeometry geoView;
  geoView.SetGeometry(&geo);
  geoView.SetCanvas(cD);
  geoView.Plot2d();

  constexpr bool twod = true;
  constexpr bool drawaxis = false;
  driftView.Plot(twod,drawaxis);

  std::cout<<"End of Program"<<std::endl;

  app.Run(true);
  return 0;
}


