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
#include "Garfield/ComponentUser.hh"
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
  
  // Compute potentials and weighting potentials (for signal calculation) with analtic expressions.   
  // Make a component for the drift field;
  ComponentUser cmp;
  cmp.SetGeometry(&geo);
  // Parametrisation of the potential.
  auto potential = [](const double x, const double y, const double z, double& V) {
	const double eps0 = 8.85e-12 * (10**15) / (10**2);  // permittivity of free space in [fC V^-1 cm^-1].
	const double l = 0.5;                               // distance between anode plane and cathode plane [cm].
	const double s = 0.1;                               // interval bwtween anode wires [cm].
	const double a = 0.001;                             // anode wire radius [cm].
	const double C = 2 * M_PI * eps0 / (M_PI*l/s - log(2*M_PI*a/s));  // capacitace per unit length.
	const double V0 = 4000.;                              // potential difference between cathode and anode.
	const double z0 = 0., const double z1 = 0.5, const double z2 = 1., const double z3 = 1.5, const double z4 = 2.; // z positions of planes.

	if (z0 <= z < z1) { 
		V = C*V0/(4*M_PI*eps0) * (2*M_PI*l/s - log(4*(sin(M_PI*x/s)*sin(M_PI*x/s) + sinh(M_PI*(z - l)/s)*sinh(M_PI*(z - l)/s))));
	} else if (z1 <= z < z2) {
		V = C*V0/(4*M_PI*eps0) * (2*M_PI*l/s - log(4*(sin(M_PI*x/s)*sin(M_PI*x/s) + sinh(M_PI*(z2 - z - l)/s)*sinh(M_PI*(z2 - z - l)/s))));
	} else if (z2 <= z < z3) {
		V = C*V0/(4*M_PI*eps0) * (2*M_PI*l/s - log(4*(sin(M_PI*y/s)*sin(M_PI*y/s) + sinh(M_PI*(z - 2*l)/s)*sinh(M_PI*(z - 2*l)/s))));
	} else if (z3 <= z <= z4) {
		V = C*V0/(4*M_PI*eps0) * (2*M_PI*l/s - log(4*(sin(M_PI*x/s)*sin(M_PI*x/s) + sinh(M_PI*(z4 - z - l)/s)*sinh(M_PI*(z4 - z - l)/s))));
	} else {
		V = 0.;
	}	
  };
  cmp.SetPotential(potential);

  // Parametrisation of the weighting potential.
  auto weightingPotential = [](const double x, const double y, const double z, double& weightingV) {
	const double V0 = 4000.;  // potential difference between anode and cathode.
	const double s = 0.1;     // internval between anode wires [cm].
	const double a = 0.001;   // anode wire radius [cm].
	const double l = 0.5;     // distance between anode plane and cathode plane [cm].
	const double rc = s/(2*M_PI) * exp(M_PI*l/s);   // radius of tube surrounding each anode wire [cm].
	const double z0 = 0., const double z1 = 0.5, const double z2 = 1., const double z3 = 1.5, const double z4 = 2.; // z positions of planes.
	double x0 = 0.;


	if (x > 0.) {
		x0 = x - floor(x);
	} else {
		x0 = x - ceil(x);
	}


	if (z0 <= z < z2) {
		r = ;
		if (r < s/2) {
			weightingV = ;
		} else {
			weightingV = ;
		}
        } else if (z2 <= z <= z4) {
		r = ;
		if (r < s/2) {
                        weightingV = ;
                } else {
                        weightingV = ;
                }
        } else {
                weightingV = 0.;
        }
  };
  cmp.SetWeightingPotential(weightingPotential);


  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "anodeX");
  sensor.AddElectrode(&cmp, "anodeY");
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
  int st = 0;
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
  ViewField potViewZ(&cmp);
  potViewZ.SetArea(-5., -5., 0., 5., 5., 2.);
  potViewZ.SetPlane(0., 0., 1., 0., 0., 0.5); // Plot at (0,0,0.5) with +ve z axis out of the page.
  potViewZ.PlotContour("v");
  // Actual potential in x-z plane.
  ViewField potViewX(&cmp);
  potViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  potViewX.SetPlane(0., -1., 0., 0., 0., 1.); // Plot at (0,0,1.) with +ve y axis into the page.
  potViewX.PlotContour("v");
  // Weighting potential around anodeX.
  //ViewField wpotViewX(&cmp);
  //wpotViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  //wpotViewX.SetPlane(0., -1., 0., 0., 0., 1.); // Plot at (0,0,1.) with +ve y axis into the page.
  //wpotViewX.PlotContourWeightingField("anodeX","v");
  // Actual potential in y-z plane.
  ViewField potViewY(&cmp);
  potViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  potViewY.SetPlane(-1., 0., 0., 0., 0., 1.); // Plot at (0,0,1.) with +ve x axis into the page.
  potViewY.PlotContour("v");
  // Weighting potential around anodeY.
  //ViewField wpotViewY(&cmp);
  //wpotViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  //wpotViewY.SetPlane(-1., 0., 0., 0., 0., 0.); // Plot at (0,0,1.) with +ve x axis into the page.
  //wpotViewY.PlotContourWeightingField("anodeY","v");

  // Plot field and weighting field profiles along z.
  ViewField fieldView(&cmp);
  fieldView.PlotProfile(0., 0., 0., 0., 0., 2.);
  // Actual field in x-y plane.
  ViewField fieldViewZ(&cmp);
  fieldViewZ.SetArea(-5., -5., 0., 5., 5., 2.);
  fieldViewZ.SetPlane(0., 0., 1., 0., 0., 0.5);  // Plot at (0,0,0.5) with +ve z axis out of the page.
  fieldViewZ.PlotContour("e");
  // Actual field in x-z plane.
  ViewField fieldViewX(&cmp);
  fieldViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  fieldViewX.SetPlane(0., -1., 0., 0., 0., 1.);  // Plot at (0,0,1) with +ve y axis into the page.
  fieldViewX.PlotContour("e");
  // Weighting field around anodeX.
  //ViewField wfieldViewX(&cmp);
  //wfieldViewX.SetArea(-5., -5., 0., 5., 5., 2.);
  //wfieldViewX.SetPlane(0., -1., 0., 0., 0., 1.);  // Plot at (0,0,1) with +ve y axis into the page.
  //wfieldViewX.PlotContour("e");
  //wfieldViewX.PlotContourWeightingField("anodeX","e");
  // Actual field in y-z plane.
  ViewField fieldViewY(&cmp);
  fieldViewY.SetArea(-5., -5., 0., 5., 5., 2.);
  fieldViewY.SetPlane(-1., 0., 0., 0., 0., 1.5);  // Plot at (0,0,1.5) with +ve x axis out of the page.
  fieldViewY.PlotContour("e");
  // Weighting field around anodeY.
  //ViewField wfieldViewY(&cmp);
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


