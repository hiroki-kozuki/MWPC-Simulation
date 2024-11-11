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
  gas.LoadGasFile("ar_50_co2_50.gas");
  auto installdir = std::getenv("GARFIELD_INSTALL");
  if (!installdir) {
    std::cerr << "GARFIELD_INSTALL variable not set.\n";
    return 1;
  }
  const std::string path = installdir;
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_CO2+_CO2.txt");
  
  //Plot drift velocity vs. E field of ions:
  ViewMedium mediumView;
  mediumView.SetMedium(&gas);
  mediumView.PlotElectronVelocity('e');

  // Define conductor material. Change to tungsten.
  MediumConductor W;

  // Geometry.
  GeometrySimple geo;

  // Coordinates of centre of gravity of film 1 [cm]:
  const double cog_x1 = 0., cog_y1 = 0., cog_z1 = -0.5; // 5mm appart from anode plane.
  // Coordinates of centre of gravity of film 2 [cm]:
  const double cog_x2 = 0., cog_y2 = 0., cog_z2 = 0.5;  // 5mm appart from anode plane.
  // Half lengths of film 1 [cm]:
  const double hl_x1 = 5., hl_y1 = 5., hl_z1 = 0.00125; // 25 microns thick
  // Half lengths of film 2 [cm]:
  const double hl_x2 = 5., hl_y2 = 5., hl_z2 = 0.00125; // 25 microns thick
  SolidBox cathode1(cog_x1, cog_y1, cog_z1, hl_x1, hl_y1, hl_z1);
  cathode1.SetBoundaryPotential(-3500.);
  // cathode1.SetLabel(const std::string& label);
  SolidBox cathode2(cog_x2, cog_y2, cog_z2, hl_x2, hl_y2, hl_z2);
  // cathode2.SetLabel(const std::string& label);
  cathode2.SetBoundaryPotential(-3500.);
  geo.AddSolid(&cathode1, &W);
  geo.AddSolid(&cathode2, &W);

  // Dimensions of single anode wire [cm]:
  const double radius = 0.001;                   // 20 microns thick.
  const double halflength = 5.;

  std::vector<SolidWire*> wires(100, nullptr);
  for (int i = 0; i < 100; ++i) {                // 100 wires.
    int j = i - 50;                              // Offset by 50 wires to correctly centre the anode plane.
    double xpos = 0.05 + (double)j*0.1;          // 1 mm spacing along x. Offset by 0.5 mm to centre the anode plane.
    wires[i] = new SolidWire(xpos, 0.0, 0.0, radius, halflength, 0, 1, 0);
    wires[i]->SetBoundaryPotential(0.);          // Grounded anode plane.
    wires[i]->SetLabel("anode");
    geo.AddSolid(wires[i], &W);
  }
  // Add boundary volume filled with gas.
  geo.SetMedium(&gas);

  // Plot device geometry in 2D
  ViewGeometry geomView2d(&geo);
  geomView2d.SetArea(-10, -10, -1, 10, 10, 2);
  geomView2d.SetPlane(0, 1, 0, 0, 0, 0.0);
  geomView2d.Plot2d();

  // Plot device geometry in 3D
  ViewGeometry geomView3d(&geo);
  geomView3d.Plot();

  // Compute field with neBEM.
  ComponentNeBem3d nebem;
  nebem.SetGeometry(&geo);
  nebem.SetTargetElementSize(0.1);               // Linear size of elements measured along their edges.
  nebem.SetMinMaxNumberOfElements(3, 15);        // Smallest and largest # of elements.
  nebem.SetNumberOfThreads(8);
  // box1.SetDiscretisationLevel();
  // box2.SetDiscretisationLevel();
  // wires.SetDiscretisationLevel();
  nebem.EnableDebugging();
  nebem.Initialise();
 
  Medium* medium = nullptr; 
  double ex = 0., ey = 0., ez = 0., v = 0.;
  int status = 0;
  nebem.ElectricField(0, 0, 0, ex, ey, ez, v, medium, status);
  std::printf("E = (%15.8f, %15.8f %15.8f), V = %15.8f, status = %d\n", ex, ey, ez, v, status);

  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&nebem);
  sensor.AddElectrode(&nebem, "anode");
  // Set the time bins [ns]
  const unsigned int nTimeBins = 500;  // 500 bins.
  const double tmin = 0.;              // From 0 ns.
  const double tmax = 350.;            // Upto 350 ns.
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
  drift.EnableAvalanche();
  drift.EnableSignalCalculation();

  // Set up the viewer of the simulated drift lines.
  // Create a new canvas.
  TCanvas* cD = new TCanvas("cD"," ",600,600);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  drift.EnablePlotting(&driftView);
  track.EnablePlotting(&driftView);

  // Set the initial position of particle: starting at z0 = -4.9 mm along perpendicular axis of MWPC.
  const double x0 = 0., y0 = 0., z0 = -0.49, t0 = 0.;
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

  // Plot potential contour
  ViewField potView(&nebem);
  potView.SetArea(-5., -5., -0.5, 5., 5., 0.5);
  potView.SetPlane(0., -1., 0., 0., 0., 0.0); // Plot at (0,0,0) with +ve y axis into the page.
  potView.PlotContour();

  // Plot field profile along z and contour in x-z plane.
  ViewField fieldView(&nebem);
  fieldView.PlotProfile(0., 0., -0.5, 0., 0., 0.5);
  fieldView.SetArea(-5., -5., -0.5, 5., 5., 0.5);
  fieldView.SetPlane(0., -1., 0., 0., 0., 0.0);  // Plot at (0,0,0) with +ve y axis into the page.
  fieldView.PlotContour("e");

  // Plot the induced current in the anode wires.
  ViewSignal signalView;
  signalView.SetSensor(&sensor);
  signalView.PlotSignal("anode");

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


