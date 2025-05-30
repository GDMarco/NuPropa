#include "tools.hh"
#include "ME2_Analytic.hh"
#include "var.hh"
#include "dsigma.hh"

using namespace std;
using namespace Variables;



// The previous routines were becoming cumbersome, so they are now updated with a channel dependent integer
double dsigma_channels( KinematicData &Kin, int channel_id ){

        // The summary of all channels is stored via the global map<int,str> process_map
        // The general form of the differential cross-section is always:
        // d sigma = \sum |M|^2 / flux * dphi_n

        // The phase space dphi_n is generated according to whether we are in 2to2 or 2to3, accessed via Kin.weight() which contains all jacobian factors
        // The flux is always 1 / (4p1.p2) for massless particles, accessed via Kin.flux()

        // Thus, here we evaluate |M|^2 given a channel id
        double ME2(0.);

        // All equal flavour results obtained from |M|^2 for nu1 + nu1 > nu1 + nu1
        // 1) nu1 + nu1 > nu1 + nu1
        if( channel_id == 1 )  ME2 = ME2_Analytic::nu1nu1_nu1nu1(1,2,3,4, Kin, pdg_projectile, 12);
        // 2) nu1 + nu1bar > nu1 + nu1bar
        if( channel_id == 2  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(1,4,3,2, Kin, pdg_projectile, 12);
        // 3) nu1bar + nu1bar > nu1bar + nu1bar
        if( channel_id == 3  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(3,4,1,2, Kin, pdg_projectile, 12);
        // 4) nu1bar + nu1 > nu1bar + nu1
        if( channel_id == 4  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(3,2,1,4, Kin, pdg_projectile, 12);
        // Mixed flavour results obtained from |M|^2 for nu1 + nu2 > nu1 + nu2
        // 5) nu1 + nu2 > nu1 + nu2        
        if( channel_id == 5  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(1,2,3,4, Kin, pdg_projectile, 12);
        // 6) nu1 + nu2bar > nu1 + nu2bar        
        if( channel_id == 6  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(1,4,3,2, Kin, pdg_projectile, 12);
        // 7) nu1bar + nu2 > nu1bar + nu2        
        if( channel_id == 7  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(3,2,1,4, Kin, pdg_projectile, 12);
        // 8) nu1bar + nu2bar > nu1bar + nu2bar        
        if( channel_id == 8  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(3,4,1,2, Kin, pdg_projectile, 12);
        // neutrino annihilation channels (3 channels)        
        // 9) nu1 + nu1bar > nu2 + nu2bar
        if( channel_id == 9  ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 14);
        // 10-12) nu1 + nu1bar > l1 + l1bar (provide outgoing fermion PDG)
        if( channel_id == 10 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, 11);
        if( channel_id == 11 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, 13);
        if( channel_id == 12 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, 15);                
        // 13-18) nu1 + nu2bar > l1 + l2bar [CC process, no Z-coupling info needed]
        if( channel_id > 12 and channel_id < 19 )  ME2 = ME2_Analytic::nu1nu2bar_l1l2bar(1,2,3,4, Kin );
        // 19-27 nu nubar > f fbar
        if( channel_id == 19 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 1);
        if( channel_id == 20 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 2);
        if( channel_id == 21 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 3);
        if( channel_id == 22 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 4);
        if( channel_id == 23 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 5);
        if( channel_id == 24 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 6);
        if( channel_id == 25 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 11);
        if( channel_id == 26 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 13);
        if( channel_id == 27 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 15);
        // 28-33
        if( channel_id > 27 and channel_id < 34 ) ME2 = ME2_Analytic::nugumma_Wl(1,2,3,4, Kin);

        // testing channel, for charm production
        if( channel_id == 99 ) ME2 = ME2_Analytic::eeB0g0NCM(1,2,3,4, Kin, 11, 4);

        // Channels 101+
        // 101-109) nu_1 gamma -> nu_1 l_2 l_2~
        if( channel_id == 101 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 11);
        if( channel_id == 102 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 13);
        if( channel_id == 103 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 15);
        if( channel_id == 104 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 1);
        if( channel_id == 105 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 2); // Here there is some problem
        if( channel_id == 106 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 3);
        if( channel_id == 107 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 4);
        if( channel_id == 108 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 5);
        if( channel_id == 109 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, 12, 6);
        // 110-115) nu1  + gamma > l1 + nu2 + l2~
        if( channel_id > 109 and channel_id < 116 ) ME2 = ME2_Analytic::nu1gamma_l1nu2l2x(1, 2, 3, 4, 5, Kin);
        // 116-118) nu1  + gamma > l1 + nu1 + l1~ [same flavour]
        if( channel_id > 115 and channel_id < 119 ){
                // Currently only one I have not validated/cross-checked, just use Recola for now
                // ME2 = 2.0 * ME2_Analytic::compute_process_recola(Kin,channel,0);
                cerr << "Check calculation again\n";
                exit(0);
        }
        // 119-121, nu1  + gamma > l1 + u d~ (massless quarks)        
        if( channel_id > 118 and channel_id < 122 ) ME2 = ME2_Analytic::nu1gamma_l1qqbar(1, 2, 3, 4, 5, Kin);

        // Implement more 2to3 channels here...
        if( channel_id > 122 ){
                cout << "Computation still in progress, use recola for now\n";
                abort();
        }

        return ME2 * Kin.weight() / Kin.flux();
}



// Implementation of inclusive result for channels 1,2,6 (corresponding to four neutrino scattering)
// Follows the results in https://cds.cern.ch/record/244041/files/PhysRevD.47.5247.pdf
// Note I analytically expand these results to be stable in the limit s12 << mz^2

// Eq. 2.1
double F0(double s,int pdg){

        if( abs(pdg) != 12 ){
                cerr << "FO: function implemented only for neutrinos currently\n";
                abort();
        }
        double t3 = 1./2.;
        complex<double> prop = 1. / ( s - MZ2C );
        complex<double> spz = prop * conj(prop) * pow(mz,4) * s;
        // 
        double sigma = 2. * gf * gf * spz.real() * pow(t3 ,2) / ( 3. * pi );
        return sigma * hbarc2;   
}

// Eq 2.2
double F1(double s){
        double prefac = pow(gf,2) * s / (2. * pi);
        double y = s / pow(mz,2);
        double F1 = ( pow(y,2) + 2. * y - 2. * (1.+y)* log1p(y) ) / pow(y,3);
        if( y < 1e-4 ){
                F1 = 1. / 3. - y / 6. + y*y/10. - pow(y,3)/15.; // O(y^4)
        }
        double sigma = prefac * F1;
        return sigma * hbarc2;
}

// Eq 2.3
double F2(double s){
        double prefac = pow(gf,2) / (2. * pi);
        double y = s / pow(mz,2);
        double F2 = ( 3. * pow(y,2) + 2. * y - 2. * pow( 1.+y, 2 )* log1p(y) ) / pow(y,3);
        if( y < 1e-4 ){
                F2 = -2. / 3. + y / 6. - y*y/15. + pow(y,3)/30.; // O(y^4)
        }
        // Factor of PZ ( s-mz^2 )
        complex<double> prop = 1. / ( s - MZ2C );
        complex<double> pz = prop * conj(prop) * pow(mz,4);        
        double sigma = prefac * pz.real() * ( s - pow(mz,2) ) * y * F2;
        return sigma * hbarc2;
}

// Eq 2.9
double nui_nui(double s){
        double prefac = pow(gf,2) / (2. * pi);
        // Full expression
        double mz2 = pow(mz,2);
        //
        double expr =  mz2 * (s / (s + mz2) + 2. * mz2 / ( 2. *  mz2 + s ) * log1p(s/mz2) );
        if( s/mz2 < 1e-4 ){
                expr = 2. * s - 2. * pow(s,2)/mz2;
        }
        double sigma = prefac * expr;
        return sigma * hbarc2;

}

// Implemented for channels: 1, 2, 6
double sigma_nu_incl(double shat, int chan){

        if( chan == 1 ) return nui_nui(shat);
        if( chan == 2 ) return F0(shat,12) + F1(shat) + F2(shat);
        if( chan == 6 ) return F1(shat);
        if( chan == 9 ) return F0(shat,12);
        return 0.0;
}


// The neutrino + photon > W + l, cross-section Seckl 9709290 - Eq. (2)

double sigma_Wl_incl(double shat, int chan){

        double msql(0.);

        if( chan == 28 or chan == 31 ){
                msql = pow(me,2);
        }
        else if( chan == 29 or chan == 32 ){
                msql = pow(mm,2);
        }
        else if( chan == 30 or chan == 33 ){
                msql = pow(mtau,2);
        }
        else{
                cerr << "sigma_Wl_incl: unsupported channel " << chan << endl;
                abort();
        }

        double shat_min = pow( mw + sqrt(msql), 2 );
        if( shat < shat_min ) return 0.0;

        double y = shat / pow(mw,2);

        double sigma = ( 2. * ( 1. - 1./y ) * ( 1. + 2. / pow(y,2) - 1. / pow(y,2) * log(y) )
                + 1. / y * ( 1. - 2. / y + 2. / pow(y,2) ) * log( pow(mw,2) / msql * pow( y - 1.,2) / y ) );

        return sqrt(2.) * ALPHA.real() * gf * sigma * hbarc2;
}


