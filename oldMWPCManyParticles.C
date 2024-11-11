#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>
#include <random>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewCell.hh"

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

  // Basic file operations
  #include <iostream>
  #include <fstream>
  using namespace std;  

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
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_CO2+_CO2.txt");
  
  //Plot drift velocity vs. E field of ions:
  ViewMedium mediumView;
  mediumView.SetMedium(&gas);
  mediumView.PlotElectronVelocity('e');
  // mediumView.PlotElectronTownsend();

  // Define conductor material. Change to tungsten.
  MediumConductor W;

  // Geometry.
  GeometrySimple geo;

  // Coordinates of centre of gravity of film 1 [cm]:
  const double cog_x1 = 0., cog_y1 = 0., cog_z1 = -1.; // 5mm appart from anode plane.
  // Coordinates of centre of gravity of film 2 [cm]:
  const double cog_x2 = 0., cog_y2 = 0., cog_z2 = 0.;  // 5mm appart from anode plane.
  // Coordinates of centre of gravity of film 3 [cm]:
  const double cog_x3 = 0., cog_y3 = 0., cog_z3 = 1.;  // 5mm appart from anode plane.
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
    anodeWiresX[i] = new SolidWire(xpos, 0.0, -0.5, radius, halflength, 0, 1, 0); // Last 3 inputs define the axis of wires.
    anodeWiresX[i]->SetBoundaryPotential(0.);          // Grounded anode plane.
    anodeWiresX[i]->SetLabel("anodeX" + std::to_string(i));
    geo.AddSolid(anodeWiresX[i], &W);
  }
  // Similarly, generate 100 anode wires for y position.
  std::vector<SolidWire*> anodeWiresY(100, nullptr);
  for (int i = 0; i < 100; ++i) {                
    int j = i - 50;                              
    double ypos = 0.05 + (double)j*0.1;          
    anodeWiresY[i] = new SolidWire(0.0, ypos, 0.5, radius, halflength, 1, 0, 0); 
    anodeWiresY[i]->SetBoundaryPotential(0.);
    anodeWiresY[i]->SetLabel("anodeY" + std::to_string(i));
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
  for (int i = 0; i < 100; ++i) {                // 100 wires.
    sensor.AddElectrode(&nebem, "anodeX" + std::to_string(i));
    sensor.AddElectrode(&nebem, "anodeY" + std::to_string(i));
  }
  // Set the time bins [ns]
  const unsigned int nTimeBins = 500;  // 500 bins.
  const double tmin = 0.;              // From 0 ns.
  const double tmax = 500.;            // Upto 500 ns.
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);
  sensor.ClearSignal();

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
  drift.UseWeightingPotential(true); // false
  drift.EnableAvalanche();
  drift.EnableSignalCalculation();
  drift.EnableIonTail();

  // Set up the viewer of the simulated drift lines.
  // Create a new canvas.
  TCanvas* cD = nullptr;
  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) { 
  	 cD = new TCanvas("cD","",600,600);
  	 driftView.SetCanvas(cD);
  	 drift.EnablePlotting(&driftView);
  	 track.EnablePlotting(&driftView);
  }

  //TCanvas* cS = nullptr;
  //ViewSignal signalView;
  //constexpr bool plotSignal = true;
  //if (plotSignal) {
	  //cS = new TCanvas("cS", "", 600, 600);
	  //signalView.SetSensor(&sensor);
	  //signalView.SetLabelY("Signal [fc/ns]")
  //}
  
  // Create a new canvas to draw multiple signals from anodeX.
  TCanvas* SignalX = nullptr;
  ViewSignal signalViewX;
  constexpr bool plotSignal = true;
  if (plotSignal) {
          SignalX = new TCanvas("SignalX", "", 600, 600);
          signalViewX.SetCanvas(SignalX);
          signalViewX.SetSensor(&sensor);
  }

  // Set the initial position of particle: starting at z0 = -9.9 mm along perpendicular axis (z axis) of MWPC.
  double x0 = 0.05, y0 = 0.05, z0 = -0.99, t0 = 0.;
  // Set the track direction and define the track: 
  // Passing through centre of MWPC perpendiculat to its plane.
  const double dx0 = 0., dy0 = 0., dz0 = 1.;
  const unsigned int nTracks = 1; // 10000 // test with 1 track only.
  const double Average_timeSpace = 25; //[ns] // 25 ns = average time space between proton bunches in the LHC ///////////////////////////////////////////////////////////////////// Check with Robert why this is 0.

  double elapsed = 0;
  double mu = 0;
  double std = 5;

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(mu, std);

  for (unsigned int j = 0; j < nTracks; ++j) {

	// Generate initial position of particle.
  	x0 = 0.05, y0 = 0.05; // x0 = distribution(generator), y0 = distribution(generator); // For nTracks > 1.
  	// If initial particle is outside 9.5 cm x 9.5 cm active area regenerate coordinates.
  	while (abs(x0) > 9.5 && abs(y0) > 9.5) {x0 = distribution(generator), y0 = distribution(generator);}

	double time_to_event = -Average_timeSpace * log(static_cast<double>(std::rand()) / RAND_MAX); ////////////////////////////////////////////////////////////////////////////////////// CHECK THIS LINE with Robert. What is it exactly?
	elapsed = elapsed + time_to_event;

	double PrimaryElectronCount = 0;
	double SecondaryElectronCount = 0;

	if (abs(x0) <= 9.5 && abs(y0) <= 9.5) {
		
		// Create and open an empty file to store signal data of "jth" particle for all anodeX and anodeY wires.
		//std::string filename = "hits_particle_" + std::to_string(j) + ".txt";
        	std::ofstream outFile1;
        	outFile1.open("hits_particle_" + std::to_string(j) + ".txt", std::ios::out);

		track.NewTrack(x0, y0, z0, elapsed, dx0, dy0, dz0);
		int loopcounter = 0;
		for (const auto& cluster : track.GetClusters()) {
			  
			for (const auto& electron : cluster.electrons) {
				drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);
				PrimaryElectronCount = PrimaryElectronCount + 1;
				SecondaryElectronCount = SecondaryElectronCount + drift.GetGain()*drift.GetLoss();
				if (loopcounter == 0) {
					std::cout << "Gain for electron: " << drift.GetGain() << "\n";
					std::cout << "Loss for electron: " << drift.GetLoss() << "\n";
				}
	
				loopcounter = loopcounter + 1;
			}	  
		}

		// Loop over 100 wires in anodeX and anodeY:
       		// Each data file should contain results for kth wire of both anodeX and anodeY.
       		for (int k = 0; k < 100; ++k) {
               		sensor.ClearSignal();

			double MeanSignalX = 0, MeanSignalY = 0;
			for (unsigned int i = 0; i < nTimeBins; ++i) {
				double feX = (sensor.GetElectronSignal("anodeX" + std::to_string(k), i));
				double feY = (sensor.GetElectronSignal("anodeX" + std::to_string(k), i));
				MeanSignalX = MeanSignalX + feX, MeanSignalY = MeanSignalY + feX;
			} 
			MeanSignalX = MeanSignalX / nTimeBins, MeanSignalY = MeanSignalY / nTimeBins;
			  
			double ChargeIntegralX = 0, ChargeIntegralY = 0;
			for (unsigned int i = 0; i < nTimeBins; ++i) {
				double fX = (sensor.GetSignal("anodeX" + std::to_string(k), i));
				double fY = (sensor.GetSignal("anodeY" + std::to_string(k), i));
				ChargeIntegralX = ChargeIntegralX + (tstep * fX), ChargeIntegralY = ChargeIntegralY + (tstep * fY);
			}
			ChargeIntegralX = ChargeIntegralX / 1e-15, ChargeIntegralY = ChargeIntegralY / 1e-15; // Transform to C from fC.
			  
			std::cout << "Run #: " << j << " of " << nTracks << "\n";
			std::cout << "Wires #: " << k << " of " << 99 << "\n";
			std::cout << "Mean signal in X: " << MeanSignalX << "\n";
			std::cout << "Mean signal in Y: " << MeanSignalY << "\n";
			std::cout << "Primary Electrons: " << PrimaryElectronCount << "\n";
			std::cout << "Secondary Electrons: " << SecondaryElectronCount << "\n";
			// std::cout << "Gain: " << SecondaryElectronCount / PrimaryElectronCount << "\n";
			std::cout << "Total Charge in X: " << ChargeIntegralX << "C \n";
			std::cout << "Total Charge in Y: " << ChargeIntegralY << "C \n";
			std::cout << "Total # of electrons moved at anodeX: " << -ChargeIntegralX / 1.60217663e-19 << "\n";
			std::cout << "Total # of electrons moved at anodeY: " << -ChargeIntegralY / 1.60217663e-19 << "\n";
			// std::cout << "Gain at anodeX (primary charge to charge on capacitor): " 
			//           << (-ChargeIntegralX / 1.60217663e-19) / PrimaryElectronCount << "\n";
			// std::cout << "Gain at anodeY (primary charge to charge on capacitor): " 
			//           << (-ChargeIntegralY / 1.60217663e-19) / PrimaryElectronCount << "\n";
			// std::cout << "Gain at anodeX (secondary charge to charge on capacitor): " 
			// 	     << (-ChargeIntegralX / 1.60217663e-19) / SecondaryElectronCount << "\n";
			// std::cout << "Gain at anodeY (secondary charge to charge on capacitor): " 
			// 	     << (-ChargeIntegralY / 1.60217663e-19) / SecondaryElectronCount << "\n";
			std::cout << k << "," << x0 << "," << y0 << "," << time_to_event << "\n";
			outFile1 << k << "," << x0 << "," << y0 << "," << time_to_event << "," << PrimaryElectronCount << "," 
				 << ChargeIntegralX << "," << ChargeIntegralY << "," << SecondaryElectronCount << "\n";

			if (plotSignal) {
				// Draw the first signal and keep the canvas.
  				SignalX->cd(); // Select the canvas.
				signalViewX.SetCanvas(SignalX);
  				signalViewX.PlotSignal("anodeX" + std::to_string(k));
  				// Draw the remaining signals without clearing the canvas.
  				SignalX->Update(); // Update the canvas to display the plots.

  			}			
	  	}
	        // int nt = 0;
         	// if (!sensor.ComputeThresholdCrossings(-2., , nt)) continue;

		//Close data file of "jth" particle for all anodeX and anodeY wires.
		outFile1.close();

		//std::ofstream outfile2;
	  }

  }




  
  
  // Loop over the clusters.
  //double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
  //int nc = 1; //0
  //while(track.GetCluster(xc,yc,zc,tc,nc,ec,extra)){
	  // Loop over the electrons in the cluster.
          //for(int k = 0; k<nc; ++k){
		  // Electron data vars. are defined
		  //double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
                  //double dx = 0., dy = 0., dz = 0.;
		  //double enx = 0., eny = 0., enz = 0., ent = 0.;
		  //int st = 0.;
		  //track.GetElectron(k, xe, ye, ze, te, ee, dx, dy, dz);
		  //drift.DriftElectron(xe, ye, ze, te);
		  //drift.GetEndPoint(enx, eny, enz, ent, st);
		  //std::cout<<"x: "<<enx<<" y: "<<eny<<" z: "<<enz<<" t: "<<ent<<" stat: "<<st<<std::endl;
	  //}
  //}
  
  //const double rTrak = ;
  //double x0 = rTrack;
  //double y0 = ;
  //const unsigned int nTracks = 100;
  //const double AverageTimeSpace = 0; // [ns]

  //double elapsed = 0;
  //double mu = 0;
  //double std = 1;

  //std::default_random_engine generator;
  //std::normal_distribution<double> distribution(mu, std);

  //std::ofstream outfile1;
  //outfile1.open("hits.txt", std::ios::out);

  //for (unsigned int j = 0; j < nTracks; ++j) {
	  //sensor.ClearSignal();

	  //x0 =  
	  //while 

  //std::cout<<"Induced charge: "<<sensor.GetTotalInducedCharge("anode")<<std::endl;

  
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

  //Set the object with which to plot and canvas.
  ViewGeometry geoView;
  geoView.SetGeometry(&geo);
  geoView.SetCanvas(cD);
  geoView.Plot2d();

  constexpr bool twod = true;
  constexpr bool drawaxis = false;
  driftView.Plot(twod,drawaxis);

  // Create a new canvas to draw multiple signals from anodeX.
  //TCanvas* SignalX = nullptr;
  //ViewSignal signalViewX;
  //constexpr bool plotSignal = true;
  //if (plotSignal) {
	  //SignalX = new TCanvas("SignalX", "", 600, 600);
	  //signalViewX.SetCanvas(SignalX);
	  //signalViewX.SetSensor(&sensor);
  //}
  // Draw the first signal and keep the canvas.
  //SignalX->cd(); // Select the canvas.
  //signalViewX.PlotSignal("anodeX0");
  // Draw the remaining signals without clearing the canvas. 
  //for (int i = 1; i < 100; ++i) {
    //signalViewX.PlotSignal("anodeX" + std::to_string(i), "same");
  //}
  //SignalX->Update(); // Update the canvas to display the plots.
  

  // Create a new canvas to draw multiple signals from anodeX.
  //TCanvas* SignalX = nullptr;
  //ViewSignal signalViewX;
  //constexpr bool plotSignal = true;
  //if (plotSignal) {
          //SignalX = new TCanvas("SignalX", "", 600, 600);
          //signalViewX.SetCanvas(SignalX);
          //signalViewX.SetSensor(&sensor);
  //}
  // Draw the first signal and keep the canvas.
  //SignalX->cd(); // Select the canvas.
  //signalViewX.PlotSignal("anodeX0");
  // Draw the remaining signals without clearing the canvas. 
  //for (int i = 1; i < 100; ++i) {
    //signalViewX.PlotSignal("anodeX" + std::to_string(i), "same");
  //}
  //SignalX->Update(); // Update the canvas to display the plots.

  // Create a new canvas to draw multiple signals from anodeY.
  //TCanvas* SignalY = new TCanvas("SignalY", " ", 600, 600);
  //ViewSignal signalViewY;
  //signalViewY.SetSensor(&sensor);
  // Draw the first signal and keep the canvas.
  //SignalY->cd(); // Select the canvas.
  //signalViewX.PlotSignal("anodeY0");
  // Draw the remaining signals without clearning the canvas. 
  //for (int i = 1; i < 100; ++i) {
    //signalViewY.PlotSignal("anodeY" + std::to_string(i), "same");
  //}
  //SignalY->Update(); // Update the canvas to display the plots.


  std::cout<<"End of Program"<<std::endl;

  app.Run(true);
  return 0;
}


