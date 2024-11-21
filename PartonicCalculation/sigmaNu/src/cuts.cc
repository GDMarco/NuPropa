#include <cmath>
#include <iostream>
#include "var.hh"
#include "tools.hh"
#include "cuts.hh"

using namespace Variables;

// Some particle specific cuts
bool pass_cuts_neutrino( p4vec &pnu ){
        return true;

        // Currently not set up cuts
}



// General function for applying cuts (set true if all cuts are passed)
void apply_cuts( KinematicData &Kin ){

        // How many particles in the Kinematic data structure
        // int nparticles = Kin.length();

        // Currently just test the costh13_CoM cross-section for 2to2 scattering
        if( active_costh13_min or active_costh13_max ){
                // Derive costh13
                double costh13 = Kin.p(3).pz() / Kin.p(3).modp();

                if( active_costh13_min ){
                        if( costh13 < costh13_min ) return;
                }

                if( active_costh13_max ){
                        if( costh13 > costh13_max ) return;
                }
        }

        // Otherwise these cuts are passed, and set True!
        Kin.set_cuts(true);
        return;
}