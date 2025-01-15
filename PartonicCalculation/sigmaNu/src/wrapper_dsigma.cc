#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <array>
#include <cstring>
#include <getopt.h>
// General files
#include <cuba.h>
#include "integration.hh"
#include "var.hh"
#include "cuts.hh"
#include "tools.hh"
// Process files
#include "dsigma.hh"
//Timing
#include <sys/time.h>

using namespace std;
using namespace Variables;

// Define some globals, to ease the set-up
// pdf_info pdf_cache;
int seed_cache;
int grid_cache = 2;

// Write more user-friendly command line reader
void print_usage() {
	cout << endl
		<< "Usage, dSigma:" << endl
		<< " -s --seed			<s>	seed number for integration\n"
		<< " -a --iobs			<a>	(0 = fixed Ecms, 1 = Ecms scan)\n"
		<< " -n --pdg_nu		<n>	pdg of (anti)neutrino projecile\n" 
		<< " -c --chan			<c>	channel selection\n"
		<< " -v --virt		  <v> virtual option for channels 1-9\n"
		<< " -h --help			<h> print available channels "
		<< endl;
	exit(0);
	return;
}

// Introduce a more user friendly interface
void read_arguments(int argc, char* argv[], int &seed, int &analysis, int &ichan, int &flav_nu, bool &virt) {
	// provide it length 4 integer array for PDG codes (these must all be entered)
	const char* const short_options = "s:a:n:c:v:h:";
	const struct option long_options[] = { { "help", 0, NULL, 'h' },
		   { "seed", 1, NULL,  's' },
		   { "analysis", 1, NULL,  'a' },
		   { "pdg_nu",1, NULL, 'n' },
		   { "chan", 1, NULL,  'c' },
		   { "virt", 1, NULL,  'v' },		   
		   { NULL, 0, NULL, 0 } };
	int next_option;
	do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
		switch (next_option) {
			case 's':
				seed = stoi(optarg, NULL);
				break;						
			case 'a':
				analysis = stoi(optarg, NULL);
				break;
			case 'c':
				ichan = stoi(optarg, NULL);
				break;				
			case 'n':
				flav_nu = stoi(optarg, NULL);
				break;
			case 'v':
				virt = stoi(optarg, NULL);
				break;
			case 'h':
				print_channels();
				abort();																
			case '?':
				print_usage();
			case -1: break;
			default: abort();
		}
	}
	while (next_option != -1);
	return;
}


// Rather than an interface to an external numerical integrator, we can also just supply the usual integration variable
// variable e.g. "cos_theta"
double dsigma_Interface_2to2(double s12, std::string variable, double var_value, int channel_number ){

	// Appropriate 2to2 channel? should be < 100
	if( channel > 100 ){
		cerr << "dsigma_Interface_2to2: requested channel " << channel << endl;
		cerr << "not a supported 2to2 process\n";
		abort();
	}

	// Check the chosen variable is suitable/supported
	// The list of variables etc. and be appropriately extended
	if( variable == "costh13_com" ){
		// Check that the values of the supplied variable make physical sense
		if( abs(var_value) > 1 ){
			cerr << "dsigma_Interface_2to2: " << variable << " out of physical range: " << var_value << endl;
			abort();
		}
	}
	else{
		cerr << "dsigma_Interface_2to2: extend interface for variable " << variable << endl;
		abort();
	}

	// 2to2 for neutrino + X > 2 scattering
	int nfinal = 2;

	// Create array of masses of final state particles
	double masses[nfinal] = {0.};
	int pdgs[nfinal] = {0};

	// Hardcode a void function for the final state particles, pass by ref
	init_channel_information( channel, pdgs, masses );

	// Sanity check on threshold of final-state particles
	double mfinal(0.);
	for( int i(0); i < nfinal; i++)
		mfinal += masses[i];
	if( pow(mfinal,2) > s12 ){
		return 0.;
	}

	// Construct the kinematics of the scattering process from the provided inputs
	// Incoming particles are always massless.
	double m12 = sqrt(s12); // Ecms^2 / 2
	// p1 and p2 already defined
	p4vec p1_CoM( m12/2, 0., 0., +m12/2 );
	p4vec p2_CoM( m12/2, 0., 0., -m12/2 );
	// p3 and p4 (four momenta depend on their masses)
	std::array<double,3> CoM_Kins = CoM_1to2_kinematics_sq( s12, pow(masses[0],2), pow(masses[1],2) );
	double E3 = CoM_Kins[0];
	double E4 = CoM_Kins[1];
	double mod_p3 = CoM_Kins[2];

	// Variable dependent implementation here
	double costh_13(0.), sinth_13(0.);
	if( variable == "costh13_com" ){
		costh_13 = var_value;
	}
	else{
		cerr << "unsupported variable " << variable << endl;
		abort();
	}
	sinth_13 = sqrt( (1.+costh_13)*(1.-costh_13) );

  // This now defines the particle four momenta in partonic CoM frame
  p4vec p3_CoM ( E3, 0, mod_p3*sinth_13, mod_p3*costh_13 );
  p4vec p4_CoM ( E4, 0,-mod_p3*sinth_13,-mod_p3*costh_13 );
  // Now place these momenta in the KinematicData structure to evaluate the |M|^2 values for any of processes
  KinematicData Kinematics(4);
  // Add partonic CoM momenta
  Kinematics.set_pi( p1_CoM, 1 );
  Kinematics.set_pi( p2_CoM, 2 );
  Kinematics.set_pi( p3_CoM, 3 );
  Kinematics.set_pi( p4_CoM, 4 );
  double scale = KinematicScales( Kinematics, scale_opt );
  Kinematics.set_muf( scale * muf_var );
  Kinematics.set_mur( scale * mur_var );
  // Have used that lambda^{1/2}(1,s_3,s_4) = 2 |p3| / sqrt(s12) in partonic CoM frame, m12=sqrt(s12)
  double ps_factor = (2.*M_PI) / ( 2. * pow(4.*M_PI,2) ) * ( 2. * mod_p3 / m12 );
  Kinematics.set_weight( ps_factor );
  Kinematics.set_flux( 2.0 * s12 );
  // Apply analysis cuts (check if phase-space point passes cuts or not)
  apply_cuts( Kinematics );
	// Check if the phase-space point passes any differential selections
	if( !Kinematics.get_cuts() ){
		return 0;
	}

	// If the phase-space point passes any restrictions, evaluate d sigma / d costh in pb
	return dsigma_channels( Kinematics, channel ) * hbarc2;
}



int Vegas_Interface(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
	// The integrand value (thing to be integrated), initialise to zero

	double dsigma_summed(0.);

	// Define some variables to control the phase-space and random variable setup
	int nrandom(0); // No additional integrations required for parton level computations

	// 2to3 if target is a photon, otherwise 2to2
	int nfinal = ( channel > 100 )? 3: 2;

	// Create array of masses of final state particles
	double masses[nfinal] = {0.};
	int pdgs[nfinal] = {0};

	// Hardcode a void function for the final state particles, pass by ref
	init_channel_information( channel, pdgs, masses );

	// Sanity check on threshold of final-state particles
	double mfinal(0.);
	for( int i(0); i < nfinal; i++)
		mfinal += masses[i];
	if( pow(mfinal,2) > Ecms2 ){
		ff[0] = 0.0;
		return 0;
	}

	// Create phase-space
	KinematicData Kin = Generate_Phase_Space( xx, nfinal, masses, nrandom, "ee" );
	// Manually set as(muR) currently not used
	Kin.set_as_mur( 0.118 );
	// Check if the phase-space point is broken or not
	if( !Kin.get_cuts() ){
		ff[0] = 0.0;
		return 0;
	}

	// New function ordered by channels (integers)
	dsigma_summed = dsigma_channels( Kin, channel );

	// abort();
	// Return the (integrand) differential cross-section in pb
	ff[0] = dsigma_summed * hbarc2;
	return 0;
}






// This is the main program, i.e. the one which is run by the executable
// we can supply it with some inputs at command line (lets just use integers)
int main(int argc, char *argv[])
{
	// Initialise the channels
	init_channels();
	if( argc < 2 ) print_usage();

  // Initialise program defaults
	init_default();
	// Some default values
	channel = 1;
	pdg_projectile = 12;
  // Initialise some process specific scales/variables
	scale_opt = -1;
	// Scale options
	mur_var = 1.;
	muf_var = 1.;	
  mu0 = mz;
	mu_loop = 100.;	// No results depend on mu_reg value
	// Set the collision environment (pp collisions at LHC 13 TeV)
	Ecms = 1e5;
	Ecms2 = pow(Ecms,2);

	int isetup(0); // If we perhaps want to select some specific processes

	// Print available channels
	print_channels();

	// Now read the specific set of command line arguments
	read_arguments(argc,argv,seed_cache,isetup,channel,pdg_projectile,active_virtual);	

	if( active_virtual and channel >= 10 ){
		cout << "Virtual should only be active for channels 1-9\n";
		abort();
	}

	// Update cuba dimensions according to the process
	update_process_dimensions();

	// Print main program settings
	print_settings();

	// Initialise recola settings
	init_recola();

	// Initialises the process registration and constructs <int,str> map
	init_recola_processes();

	////////////
	// Setups //
	////////////
	// Setup 0, most simple inclusive LHC setup
	if( isetup == 0 ){
		cout << "Inclusive cross-section at fixed Ecms = " << Ecms << endl;
	}
	else if( isetup == 1 ){
		cout << "Inclusive cross-section for an energy scan in Ecms\n";
	}

	///////////////////////////
	// Information for Vegas //
	///////////////////////////	
	double warmup_precision = 1e-2;
	double integration_precision = 1e-3;
	int grid_no = 3;
	grid_cache = grid_no;
	int max_evaluations = 2e7;

	// When evaluating virtual corrections, work with reduced precision for NLO coefficient
	if( active_virtual ){
		warmup_precision = 1e-1;
		integration_precision = 1e-2;
	}

	//////////////////////////////
	// Store cross-section data //
	//////////////////////////////
	const std::string s_analysis[3] = {"SigmaIncl","SigmaIncl_Ecms","dSigmaCosth13"};
  string outfile = s_analysis[isetup]+"_channel"+to_string(channel);
  // If active virtual, add NLO tag to file
  if( active_virtual ) outfile = outfile +"_virt";
  ofstream ofile_results;
  // open file
  ofile_results.open(outfile+"_s"+to_string(seed_cache)+".txt");
  // save general program settings
	write_settings(ofile_results,"");

	// Also write the scattering process
	ofile_results << endl << "# channel_id = " << channel << endl;
	ofile_results << "# active_virtual = " << active_virtual << endl;
	ofile_results << "# process  = " << process_map.at(channel) << endl;

	/////////////////
	// Start timer //
	/////////////////	
  struct timeval t0, t1;
  gettimeofday(&t0,NULL);


  // Fiducial results for isetup < 3
  if( isetup == 0 ){
		// Run a test integral, returns an array with two entries < integral, error >, accessed via test_setup[0] and test_setup[1] respectively.
		array<double,2> sigma_fiducial = {0.};
		//			
		sigma_fiducial = integration::vegasC4(Vegas_Interface, warmup_precision, integration_precision, cuba_dimensions, NULL, seed_cache, grid_cache, max_evaluations );

		// Output the results
		cout << "Integral = " << sigma_fiducial[0] << endl;
		cout << "Error = " << sigma_fiducial[1] << endl;	

		cout << sigma_Wl_incl( Ecms2, channel ) << endl;

		// cout << sigma_nu_incl( Ecms2, channel ) << endl;

		// Save this information to the file
		ofile_results << "# Sigma Fiducial: sigma\terror\tsigma_analytic\tratio" << endl;

		// Include a function which writes all relevant information to the text file
		ofile_results << sigma_fiducial[0] << "\t" << sigma_fiducial[1] << "\t" << sigma_nu_incl(Ecms2,channel) << "\t" << sigma_fiducial[0]/sigma_nu_incl(Ecms2,channel) << endl;
	}

	// Set up the energy scan
	if( isetup == 1 ){

		// Either derive Ecms for a varying Enu, or directly fix Ecms

		// Gaetano look at 10^8 eV^2 > 10^24 eV^2
		// Corresponds to 10^4 eV > 10^12 eV: 
		double Ecms_low = 1e-7;
		double Ecms_high = 1e5;
		// Number of bins to consider
		int n_bins = 500;
		vector<double> Ecms_values = linspace( log(Ecms_low),log(Ecms_high), n_bins);
		// Note, the cross-sections plateau above 5x10^3 GeV
		// ^^^^^^^^^^^^^^^^^^^^^^^

		// Vector to store the computed cross-section
		vector< array<double,2> > sigma;
		for( int i(0); i < n_bins; i++ ){
			// Select Ecms
			Ecms = exp(Ecms_values[i]);
			Ecms2 = pow(Ecms,2);
			// Compute the cross-section
			array<double,2> sigma_incl = integration::vegasC4(Vegas_Interface, warmup_precision, integration_precision, cuba_dimensions, NULL, seed_cache, grid_cache, max_evaluations );
			// Save the results in the vector
			sigma.push_back( sigma_incl );
		}

		ofile_results << "# Ecms\tsigma[pb]\tsigma_error[pb]\tsigma_analytic[pb]\tratio\n";
		// Write the results to the file
		for( unsigned int i=0; i < sigma.size(); i++ ){
		// for( auto i: sigma ){
			array<double,2> sig = sigma[i];

			// For performing some analytic checks: channels (1,2,6,9)
			if( channel == 1 or channel == 2 or channel == 6 or channel == 9 ){
				Ecms2 =	pow( exp(Ecms_values[i]), 2);
				double sigma_analytic = sigma_nu_incl(Ecms2,channel);
				ofile_results << exp(Ecms_values[i]) << "\t"	<< sig[0] << "\t" << sig[1] << "\t" << sigma_analytic << "\t" << sig[0]/sigma_analytic << endl;
			}
			else if( channel >= 28 and channel <= 33 ){
				Ecms2 =	pow( exp(Ecms_values[i]), 2);
				double sigma_analytic = sigma_Wl_incl(Ecms2,channel);
				ofile_results << exp(Ecms_values[i]) << "\t"	<< sig[0] << "\t" << sig[1] << "\t" << sigma_analytic << "\t" << sig[0]/sigma_analytic << endl;

			}
			else{
				ofile_results << exp(Ecms_values[i]) << "\t"	<< sig[0] << "\t" << sig[1] << endl;
			}

		}		
	}



  // Differential cross-section in dcosth_13 at fixed-energy
  if( isetup == 2 ){

  	// Perform the differential cross-section in dcosth_13
		double costh13_low = -1;
		double costh13_up = +1.;
		// Number of bins to consider
		int n_bins = 200;
		vector<double> costh13_values = linspace( costh13_low, costh13_up, n_bins);


		active_costh13_min = true;
		active_costh13_max = true;
		for( int ibin = 0; ibin < (n_bins-1); ibin++){

			costh13_min = costh13_values[ibin];
			costh13_max = costh13_values[ibin+1];

			// continue;
			double costh13_cen = ( costh13_min + costh13_max ) / 2.;
			double bin_width = costh13_max - costh13_min;


			// Compute the differential cross-section within a bin
			array<double,2> dsigma_diff = integration::vegasC4(Vegas_Interface, warmup_precision, integration_precision, cuba_dimensions, NULL, seed_cache, 0, max_evaluations );

			// Compute analytically at the bin centre
			double dsigma_analytic = dsigma_Interface_2to2(Ecms2, "costh13_com", costh13_cen, channel );

			cout << "Ratio of the results\n\n\n" << endl;
			cout << (dsigma_diff[0] / bin_width) / dsigma_analytic << endl;
		}

	}





	ofile_results.close();

  gettimeofday(&t1,NULL);
  double ta=t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*0.000001;
  cout << "Total time: " << ta << "s" << endl;		

	return 0;
}


