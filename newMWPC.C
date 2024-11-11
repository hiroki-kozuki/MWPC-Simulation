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
  gas.LoadGasFile("ar_50_co2_50_V4.gas");
  auto installdir = std::getenv("GARFIELD_INSTALL");
  if (!installdir) {
    std::cerr << "GARFIELD_INSTALL variable not set.\n";
    return 1;
  }
  const std::string path = installdir;
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt"); // Try "/share/Garfield/Data/IonMobility_CO2+_CO2.txt"
  
  //Plot drift velocity vs. E field of ions:
  ViewMedium mediumView;
  mediumView.SetMedium(&gas);
  mediumView.PlotElectronVelocity('e');

  // Define conductor material. Change to tungsten.
  MediumConductor W;

  // Geometry.
  GeometrySimple geo;

  // Dimensions of single cathode wire [cm]:
  const double cathodeRadius = 0.005;                   // 100 microns thick.
  const double cathodeHalflength = 5.;


  // Generate 100 wires for 1st cathode layer (repeating in y).
  std::vector<SolidWire*> cathodeWires1(100, nullptr);
  for (int i = 0; i < 100; ++i) {                // 100 wires.
    int j = i - 50;                              // Offset by 50 wires to correctly centre the cathode plane.
    double ypos = 0.05 + (double)j*0.1;          // 1 mm spacing along y axis. Offset by 0.5 mm to centre the cathode plane.
    cathodeWires1[i] = new SolidWire(0.0, ypos, -1.0, cathodeRadius, cathodeHalflength, 1, 0, 0); // Last 3 inputs define the axis of wires.
    cathodeWires1[i]->SetBoundaryPotential(-4000.);
    cathodeWires1[i]->SetLabel("cathode1");
    geo.AddSolid(cathodeWires1[i], &W);
  }
  // Generate 100 wires for 2nd cathode layer (repeating in y).
  std::vector<SolidWire*> cathodeWires2(100, nullptr);
  for (int i = 0; i < 100; ++i) {
    int j = i - 50;
    double ypos = 0.05 + (double)j*0.1;
    cathodeWires2[i] = new SolidWire(0.0, ypos, 0.0, cathodeRadius, cathodeHalflength, 1, 0, 0);
    cathodeWires2[i]->SetBoundaryPotential(-4000.);
    cathodeWires2[i]->SetLabel("cathode2");
    geo.AddSolid(cathodeWires2[i], &W);
  }
  // Generate 100 wires for 3rd cathode layer (repeating in x).
  std::vector<SolidWire*> cathodeWires3(100, nullptr);
  for (int i = 0; i < 100; ++i) {
    int j = i - 50;
    double xpos = 0.05 + (double)j*0.1;
    cathodeWires3[i] = new SolidWire(xpos, 0.0, 1.0, cathodeRadius, cathodeHalflength, 0, 1, 0);
    cathodeWires3[i]->SetBoundaryPotential(-4000.);
    cathodeWires3[i]->SetLabel("cathode3");
    geo.AddSolid(cathodeWires3[i], &W);
  }


  // Dimensions of single anode wire [cm]:
  const double anodeRadius = 0.001;                   // 20 microns thick.
  const double anodeHalflength = 5.;

  // Generate 100 anode wires for x position.
  std::vector<SolidWire*> anodeWiresX(100, nullptr);
  for (int i = 0; i < 100; ++i) {                // 100 wires.
    int j = i - 50;                              // Offset by 50 wires to correctly centre the anode plane.
    double xpos = 0.05 + (double)j*0.1;          // 1 mm spacing along x axis. Offset by 0.5 mm to centre the anode plane.
    anodeWiresX[i] = new SolidWire(xpos, 0.0, -0.5, anodeRadius, anodeHalflength, 0, 1, 0); // Last 3 inputs define the axis of wires.
    anodeWiresX[i]->SetBoundaryPotential(0.);          // Grounded anode plane.
    anodeWiresX[i]->SetLabel("anodeX");
    geo.AddSolid(anodeWiresX[i], &W);
  }
  // Generate 100 anode wires for y position.
  std::vector<SolidWire*> anodeWiresY(100, nullptr);
  for (int i = 0; i < 100; ++i) {                
    int j = i - 50;                              
    double ypos = 0.05 + (double)j*0.1;          
    anodeWiresY[i] = new SolidWire(0.0, ypos, 0.5, anodeRadius, anodeHalflength, 1, 0, 0); 
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

  // Compute fields with neBEM. By default, neBEM calculates both the actual fields 
  // and weighting fields (used for signal calculation).
  ComponentNeBem3d nebem;
  nebem.SetGeometry(&geo);
  nebem.SetTargetElementSize(0.1);               // Linear size of elements measured along their edges.
  nebem.SetMinMaxNumberOfElements(3, 15);        // Smallest and largest # of elements.
  nebem.SetNumberOfThreads(8);
  // box1.SetDiscretisationLevel();
  // box2.SetDiscretisationLevel();
  // wires.SetDiscretisationLevel();
  //nebem.EnableDebugging();
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
  const unsigned int nTimeBins = 500;  // 500 bins.
  const double tmin = 0.;              // From 0 ns.
  const double tmax = 1000.;            // Upto 350 ns.
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Simulate ionisation with Heed.
  std::cout<<"-- Particle track simulation --"<<std::endl;
  TrackHeed track;
  track.SetParticle("muon");
  track.SetEnergy(170e9);
  track.SetSensor(&sensor);
  // Compute the track.
  // track.EnableDebugging();

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

  // Set the initial position of particle: starting at z0 = -9.9 mm along perpendicular axis (z axis) of MWPC.
  const double x0 = 0.05, y0 = 0.05, z0 = -0.99, t0 = 0.;
  // Set the track direction and define the track: 
  // Passing through centre of MWPC perpendiculat to its plane.
  const double dx0 = 0., dy0 = 0., dz0 = 1.;
  track.NewTrack(x0,y0,z0,t0,dx0,dy0,dz0);
  
  // Loop over the clusters.
  double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
  int nc = 0;
  while(track.GetCluster(xc,yc,zc,tc,nc,ec,extra)){
	  // Loop over the electrons in the cluster.
	  for(int k = 0; k<nc; ++k){
		  // Electron data vars. are defined
		  double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
		  double dx = 0., dy = 0., dz = 0.;
		  double enx = 0., eny = 0., enz = 0., ent = 0.;
		  int st = 0.;
		  track.GetElectron(k, xe, ye, ze, te, ee, dx, dy, dz);
		  drift.DriftElectron(xe, ye, ze, te);
		  drift.GetEndPoint(enx, eny, enz, ent, st);
		  std::cout<<"x: "<<enx<<" y: "<<eny<<" z: "<<enz<<" t: "<<ent<<" stat: "<<st<<std::endl;
	  }
  }
  
  std::cout<<"Induced charge: "<<sensor.GetTotalInducedCharge("anode")<<std::endl;

  // Plot the potential and weighting potential contours.
  // Actual potential in x-y plane.
  ViewField potViewZ(&nebem);
  potViewZ.SetArea(-5., -5., -1., 5., 5., 1.);
  potViewZ.SetPlane(0., 0., 1., 0., 0., -0.5); // Plot at (0,0,-0.5) with +ve z axis out of the page.
  potViewZ.PlotContour("v");
  // Actual potential in x-z plane.
  ViewField potViewX(&nebem);
  potViewX.SetArea(-5., -5., -1., 5., 5., 1.);
  potViewX.SetPlane(0., -1., 0., 0., 0., 0.); // Plot at (0,0,0) with +ve y axis into the page.
  potViewX.PlotContour("v");
  // Weighting potential around anodeX.
  //ViewField wpotViewX(&nebem);
  //wpotViewX.SetArea(-5., -5., -1., 5., 5., 1.);
  //wpotViewX.SetPlane(0., -1., 0., 0., 0., 0.); // Plot at (0,0,0) with +ve y axis into the page.
  //wpotViewX.PlotContourWeightingField("anodeX","v");
  // Actual potential in y-z plane.
  ViewField potViewY(&nebem);
  potViewY.SetArea(-5., -5., -1., 5., 5., 1.);
  potViewY.SetPlane(-1., 0., 0., 0., 0., 0.); // Plot at (0,0,0) with +ve x axis into the page.
  potViewY.PlotContour("v");
  // Weighting potential around anodeY.
  //ViewField wpotViewY(&nebem);
  //wpotViewY.SetArea(-5., -5., -1., 5., 5., 1.);
  //wpotViewY.SetPlane(-1., 0., 0., 0., 0., 0.); // Plot at (0,0,0) with +ve x axis into the page.
  //wpotViewY.PlotContourWeightingField("anodeY","v");

  // Plot field and weighting field profiles along z.
  ViewField fieldView(&nebem);
  fieldView.PlotProfile(0., 0., -1.0, 0., 0., 1.0);
  // Actual field in x-y plane.
  ViewField fieldViewZ(&nebem);
  fieldViewZ.SetArea(-5., -5., -1.0, 5., 5., 1.0);
  fieldViewZ.SetPlane(0., 0., 1., 0., 0., -0.5);  // Plot at (0,0,-0.5) with +ve z axis out of the page.
  fieldViewZ.PlotContour("e");
  // Actual field in x-z plane.
  ViewField fieldViewX(&nebem);
  fieldViewX.SetArea(-5., -5., -1.0, 5., 5., 1.0);
  fieldViewX.SetPlane(0., -1., 0., 0., 0., 0.);  // Plot at (0,0,0) with +ve y axis into the page.
  fieldViewX.PlotContour("e");
  // Weighting field around anodeX.
  //ViewField wfieldViewX(&nebem);
  //wfieldViewX.SetArea(-5., -5., -1.0, 5., 5., 1.0);
  //wfieldViewX.SetPlane(0., -1., 0., 0., 0., 0.);  // Plot at (0,0,0) with +ve y axis into the page.
  //wfieldViewX.PlotContour("e");
  //wfieldViewX.PlotContourWeightingField("anodeX","e");
  // Actual field in y-z plane.
  ViewField fieldViewY(&nebem);
  fieldViewY.SetArea(-5., -5., -1.0, 5., 5., 1.0);
  fieldViewY.SetPlane(-1., 0., 0., 0., 0., 0.5);  // Plot at (0,0,0.5) with +ve x axis out of the page.
  fieldViewY.PlotContour("e");
  // Weighting field around anodeY.
  //ViewField wfieldViewY(&nebem);
  //wfieldViewY.SetArea(-5., -5., -1.0, 5., 5., 1.0);
  //wfieldViewY.SetPlane(-1., 0., 0., 0., 0., 0.5);  // Plot at (0,0,0.5) with +ve x axis out of the page.
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


