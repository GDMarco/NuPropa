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
#include "var.hh"
#include "tools.hh"
// Process files
#include "dsigma.hh"
//Timing
#include <sys/time.h>

#include <filesystem>
namespace fs = std::filesystem;

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
        << " -E --Ecms     <E> center of mass energy of the interaction in GeV\n"
		<< " -s --seed			<s>	seed number for integration\n"
		<< " -a --iobs			<a>	(0 = in costh, 1 = in t variable j)\n"
		<< " -n --pdg_nu		<n>	pdg of (anti)neutrino projecile\n"
		<< " -c --chan			<c>	channel selection\n"
		<< " -v --virt		  <v> virtual option for channels 1-9\n"
		<< " -h --help			<h> print available channels "
		<< endl;
	exit(0);
	return;
}

// Introduce a more user friendly interface
void read_arguments(int argc, char* argv[], double &Ecms, int &seed, int &analysis, int &ichan, int &flav_nu ) {
	// provide it length 4 integer array for PDG codes (these must all be entered)
	const char* const short_options = "E:s:a:n:c:h:";
	const struct option long_options[] = { { "help", 0, NULL, 'h' },
           { "Ecms", 1, NULL,  'E' },
           { "seed", 1, NULL,  's' },
		   { "analysis", 1, NULL,  'a' },
		   { "pdg_nu",1, NULL, 'n' },
		   { "chan", 1, NULL,  'c' },
		   { NULL, 0, NULL, 0 } };
	int next_option;
	do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
		switch (next_option) {
            case 'E':
                Ecms = stod(optarg, NULL);
                break;
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



// The following function acts as the main interface between all expressions for |M|^2 and dsigma / dX which were implented in the previous code
// The implementation of dsigma and |M|^2 relied on the construction of a KinematicData structure (four momentum of the scattering particles)
// This interface takes input values of s12, and some differential variable (e.g. t_hat or costh_CoM)
// From here it builds the kinematic data structure (explicit four momentum representation) and evaluates dsigma / dX

// After giving it the input variables (which fix the 2to2 scattering)
// It just needs to know the channel which is being evaluated "channel_number"
// The full list of channels are seen in "print_channels()"

// Rather than an interface to an external numerical integrator, we can also just supply the usual integration variable
double dsigma_Interface_2to2(double s12, std::string variable, double var_value, int channel_number ){

	// Appropriate 2to2 channel? should be < 100
	if( channel > 100 ){
		cerr << "dsigma_Interface_2to2: requested channel " << channel << endl;
		cerr << "not a supported 2to2 process currently\n";
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
	// 1) This initalises the particle masses (by channel)
	init_channel_information( channel, pdgs, masses );

	// Sanity check on threshold of final-state particles
	double mfinal(0.);
	for( int i(0); i < nfinal; i++)
		mfinal += masses[i];

	if( pow(mfinal,2) > s12 ){
		return 0.;
	}

	//////////////////////////////////////////////////////
	// 2) Construct p1 and p2, aligned along the z-axis //
	//////////////////////////////////////////////////////	

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

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 3) The last differential quantity "variable" and its value "var_value" then fixes the outoing p3, p4 momenta //
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Example a) using dcosth_CoM [-1,1] to fix the scattering
	// cos_th approach: Variable dependent implementation here
	double costh_13(0.), sinth_13(0.);
	if( variable == "costh13_com" ){
		costh_13 = var_value;
	}
	else{
		cerr << "unsupported variable " << variable << endl;
		abort();
	}
	sinth_13 = sqrt( (1.+costh_13)*(1.-costh_13) );

	// Example b) using the mandelstam "t_hat"
	// t_hat = - (p1 - p3)^2 - m1^2 - m3^3
	// One must check the value provided is physically consistent
	// t_hat = -s13 = - 2 ( E1 E3 - |p1||p3| costheta_13 )
	// Max and minimum allowed values depend on the particle masses (which is channel dependent)
	// Working in the CoM frame, we can say that:
	// Assuming massless incoming particles
	// E1 = sqrt(s12)/2 = Ecms/2
	// E3 = sqrt(s12)/2 ( 1 + m3^2 - m4^2 )
	// |p3| = |p4| = sqrt(s12)/2 * kallen(1,m3^2/s12,m4^2/s12)^(0.5)
	// Hence, max / min values of t_hat are:
	// t_hat = - s13 = -2 p1.p3
	// t_hat_{max/min} = - 2 E1 ( E3 - |p3| cos_theta ) -> - Ecms ( E3 [-/+] |p3| )
	// t_hat_{max/min} = Ecms^2/2 ( [1+m3^2/s12-m4^2/s12] -/+ Kallen(1,m3^2/s12,m4^2/s12)^(0.5) )

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
	double scale = 1.0;//KinematicScales( Kinematics, scale_opt );
	Kinematics.set_muf( scale * muf_var );
	Kinematics.set_mur( scale * mur_var );
	// Have used that lambda^{1/2}(1,s_3,s_4) = 2 |p3| / sqrt(s12) in partonic CoM frame, m12=sqrt(s12)
	double ps_factor = (2.*M_PI) / ( 2. * pow(4.*M_PI,2) ) * ( 2. * mod_p3 / m12 );
	Kinematics.set_weight( ps_factor );
	Kinematics.set_flux( 2.0 * s12 );
	// Apply analysis cuts (check if phase-space point passes cuts or not)
	// Check if the phase-space point passes any differential selections

	// If the phase-space point passes any restrictions, evaluate d sigma / d costh in pb
	return dsigma_channels( Kinematics, channel ) * hbarc2;
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
    
	int isetup(0); // If we perhaps want to select some specific processes

	// Print available channels
	print_channels();

	// Now read the specific set of command line arguments
	read_arguments(argc,argv,Ecms,seed_cache,isetup,channel,pdg_projectile);

    Ecms2 = pow(Ecms,2);
	// Print main program settings
	print_settings();
    
	////////////
	// Setups //
	////////////
	// Setup 0, most simple inclusive LHC setup
	if( isetup == 0 ){
		cout << "Implement a dcostheta scan at fixed ECMS\n";
	}
	else if( isetup == 1 ){
		cout << "Implement a dt_hat scan at fixed ECMS\n";

	}

	//////////////////////////////
	// Store cross-section data //
	//////////////////////////////
    std::string output_dir = "./dataDifferentialXS/channel"+to_string(channel)+"/";
    
    if (!fs::exists(output_dir)) {     // probably with fs::filesystem in OS
        fs::create_directories(output_dir);
    }
    
	const std::string s_analysis[3] = {"dsigdcosth", "dsigdt"};
	string outfile = s_analysis[isetup]+"_channel"+to_string(channel);
	ofstream ofile_results;
	
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(5) << Ecms2;
    std::string energyCms2 = "_EcmsSq" + ss.str();
    
    // open file
	ofile_results.open(output_dir+outfile+energyCms2+"_s"+to_string(seed_cache)+".txt");
	// save general program settings
	write_settings(ofile_results,"");

    /**
    if (!ofile_results) {
        std::cerr << "Error: could not open file for writing: " << full_filename << std::endl;
    }
    */
    
	// Also write the scattering process
	ofile_results << endl << "# channel_id = " << channel << endl;
	ofile_results << "# process  = " << process_map.at(channel) << endl;

	/////////////////
	// Start timer //
	/////////////////	
	struct timeval t0, t1;
	gettimeofday(&t0,NULL);
    
	// Differential cross-section in dcosth_13 at fixed-energy
	if( isetup == 0 ){

		// Perform the differential cross-section in dcosth_13
		double costh13_low = -1;
		double costh13_up = +1.;
		// Number of bins to consider
		int n_bins = 100;
		vector<double> costh13_values = linspace( costh13_low, costh13_up, n_bins+1);

		for( int ibin = 0; ibin < n_bins; ibin++){
			double costh13_min = costh13_values[ibin];
			double costh13_max = costh13_values[ibin+1];
			// bin centre
			double costh13_cen = ( costh13_min + costh13_max ) / 2.;
			// Compute analytically at the bin centre
			double dsigma_analytic = dsigma_Interface_2to2(Ecms2, "costh13_com", costh13_cen, channel );
			ofile_results << costh13_cen << "\t" << dsigma_analytic << endl;

		}
	}





  ofile_results.close();

  gettimeofday(&t1,NULL);
  double ta=t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*1e-6;
  cout << "Total time: " << ta << "s" << endl;		

	return 0;
}


