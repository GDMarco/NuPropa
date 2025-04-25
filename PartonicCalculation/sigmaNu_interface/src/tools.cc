#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include "tools.hh"

using namespace std;

p4vec::p4vec() {
	p = {0.,0.,0.,0.};
}

p4vec::p4vec(double E, double px, double py, double pz){
	p = {E, px, py, pz};
}

// General boost
void p4vec::boost(double bx, double by, double bz){
    // beta_v = {bx,by,pz}
    double b2 = bx*bx + by*by + bz*bz;
    double gam  = 1./sqrt(1.-b2);
    if( gam < 0. or isnan(gam) ){
        cout << gam << " " << b2 << endl;
        cout << bx*bx + by*by + bz*bz << endl;
    }
    double bp = p[1]*bx + p[2]*by + p[3]*bz;
    double gam2 = (gam-1.)/b2;
    if( b2 <= 0 ) gam2 = 0.; // For numerical stability when beta = 0.
    p[1] = p[1] + gam2*bp*bx + gam*bx*p[0];
    p[2] = p[2] + gam2*bp*by + gam*by*p[0];
    p[3] = p[3] + gam2*bp*bz + gam*bz*p[0];
    p[0] = gam*( p[0] + bp );
}

// General rotation around axis = i, by alpha (given by cos[alpha])
void p4vec::rot(double c_a, int i){
    p[0] = p[0]; // Energy component un-touched
    array< double, 3 > temp = {p[1],p[2],p[3]};
    array< array<double, 3>, 3 > rotation;
    double s_a = sqrt(1.-pow(c_a,2));
    if( i == 1 ){
        rotation[0] = { { 1.0, 0.0, 0.0} };
        rotation[1] = { { 0.0, c_a, s_a} };
        rotation[2] = { { 0.0,-s_a, c_a} };        
    }
    else if( i == 2 ){
        rotation[0] = { { c_a, 0.0,-s_a} };
        rotation[1] = { { 0.0, 1.0, 0.0} };
        rotation[2] = { { s_a, 0.0, c_a} };         
    }
    else if( i == 3 ){
        rotation[0] = { { c_a, s_a, 0.0} };
        rotation[1] = { {-s_a, c_a, 0.0} };
        rotation[2] = { { 0.0, 0.0, 1.0} };           
    }
    else{
        cout << "p4vec::rot unsuported axis rotation " << i << endl;
        abort();
    }
    for( int j = 0; j < 3; j++){
        p[j+1] = 0.;
        for( int k = 0; k < 3; k++){
            p[j+1] += rotation[j][k] * temp[k];
        }
    }
}

double p4vec::rap(){
    return (0.5*log( (p[0]+p[3])/(p[0]-p[3]) ) );
    // For better numerical stability
}

double p4vec::modp(){
    return sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
}

double p4vec::eta(){
    double modp = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
    return (0.5*log( (modp+p[3])/(modp-p[3]) ) );
}

double p4vec::pT(){
    return sqrt(p[1]*p[1]+p[2]*p[2]);
}

double p4vec::ET(){
    double ET2 = (p[0] + p[3]) * (p[0] - p[3]);
    // For numerical instabilities with extreme kinematics    
    if( ET2 < 0. ){ // Should never be encountered
        cout << "p4vec::ET repaired\n";
        cout << "new ET2 = " << abs(ET2) << endl;
        return sqrt( abs(ET2) );
    }
    return sqrt( ET2 );
}

double p4vec::m2(){
    return p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
}

double p4vec::m(){
    return sqrt(p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3]);
}

double p4vec::E(){
    return p[0];
}
double p4vec::px(){
    return p[1];
}
double p4vec::py(){
    return p[2];
}
double p4vec::pz(){
    return p[3];
}
double p4vec::plus(){// (E + Pz) / sqrt(2)
        return (p[0] + p[3])/sqrt(2.);
}
double p4vec::minus(){// (E - Pz) / sqrt(2)
        return (p[0] - p[3])/sqrt(2.);
}
// Access individual elements i = 0 (E), 1 (px), 2 (py), 3 (pz) of four-momentum
double p4vec::pi(int i){
    if( i==0 or i==4 ) return p[0];
    else if( i>0 and i<4 )
        return p[i];
    else{
        std::cout << "p4vec::pi trying to access element i = " << i << std::endl;
        return 0.0;        
    }
}

p4vec p4vec::operator+(const p4vec& vec){
    p4vec p4out;
    p4out.p[1] = this->p[1] + vec.p[1];
    p4out.p[2] = this->p[2] + vec.p[2];
    p4out.p[3] = this->p[3] + vec.p[3];
    p4out.p[0] = this->p[0] + vec.p[0];
    return p4out;
}

p4vec p4vec::operator+=(const p4vec& vec){
    p4vec p4out;
    this->p[1] += vec.p[1];
    this->p[2] += vec.p[2];
    this->p[3] += vec.p[3];
    this->p[0] += vec.p[0];
    return p4out;
}

p4vec p4vec::operator-(const p4vec& vec){
    p4vec p4out;
    p4out.p[1] = this->p[1] - vec.p[1];
    p4out.p[2] = this->p[2] - vec.p[2];
    p4out.p[3] = this->p[3] - vec.p[3];
    p4out.p[0] = this->p[0] - vec.p[0];
    return p4out;
}

p4vec p4vec::operator-=(const p4vec& vec){
    p4vec p4out;
    this->p[1] -= vec.p[1];
    this->p[2] -= vec.p[2];
    this->p[3] -= vec.p[3];
    this->p[0] -= vec.p[0];
    return p4out;
}

p4vec p4vec::operator*(double x){
    p4vec p4out;
    p4out.p[1] = this->p[1] * x;
    p4out.p[2] = this->p[2] * x;
    p4out.p[3] = this->p[3] * x;
    p4out.p[0] = this->p[0] * x;
    return p4out;
}

p4vec p4vec::operator/(double x){
    p4vec p4out;
    p4out.p[1] = this->p[1] / x;
    p4out.p[2] = this->p[2] / x;
    p4out.p[3] = this->p[3] / x;
    p4out.p[0] = this->p[0] / x;
    return p4out;
}

// General kinematic functions given below

// Kallen function
double kallen(const double x, const double y, const double z){
    return pow(x,2) + pow(y,2) + pow(z,2) - 2.*x*y - 2.*x*z - 2.*y*z;
}

// Dot product of four vectors, pi.pj
double dotp4(p4vec &pi, p4vec &pj){
    return pi.E()*pj.E()-pi.px()*pj.px()-pi.py()*pj.py()-pi.pz()*pj.pz();
}

// dy_ij = y_i = y_j
double dy_ij(p4vec &pi, p4vec &pj){
    return pi.rap()-pj.rap();
}

// deta_ij = eta_i - eta_j
double deta_ij(p4vec &pi, p4vec &pj){
    return pi.eta()-pj.eta();
}

// dphi_ij = phi_i - phi_j
double dphi_ij(p4vec &pi, p4vec &pj){
    double phi1 = atan2( pi.py(), pi.px() );
    // [0,2Pi]
    if( phi1 < 0 ) phi1 += 2*M_PI;
    double phi2 = atan2( pj.py(), pj.px() );
    if( phi2 < 0 ) phi2 += 2*M_PI;
    double dphi = fabs(phi1-phi2);
    // [0,Pi]
    if( dphi > M_PI) dphi = 2*M_PI-dphi;
    return dphi;
}

// dR_ij = sqrt( dphi_ij^2 + dy_ij^2 )
double dR_ij(p4vec &pi, p4vec &pj){
    double dy2 = pow(pi.rap()-pj.rap(),2);
    double phi1 = atan2( pi.py(), pi.px() );
    // [0,2Pi]
    if( phi1 < 0 ) phi1 += 2*M_PI;
    double phi2 = atan2( pj.py(), pj.px() );
    if( phi2 < 0 ) phi2 += 2*M_PI;
    double dphi = fabs(phi1-phi2);
    // [0,Pi]
    if( dphi > M_PI) dphi = 2*M_PI-dphi;
    return sqrt(dy2+dphi*dphi);
}

// costheta_CS, i.e. in the Collins-Soper frame
double costheta_CS(p4vec &lm, p4vec &lp){
    // Construct V
    p4vec V = lm + lp;
    double mll = V.m();
    double ptll = V.pT();
    // + and - components
    double p1p = lm.E()+lm.pz();
    double p1m = lm.E()-lm.pz();
    double p2p = lp.E()+lp.pz();
    double p2m = lp.E()-lp.pz();
    // Fix sign of pz
    double sign_pz = 1.0;
    if( V.pz() < 0 ) sign_pz *= -1.0;
    //
    return sign_pz*((p1p*p2m)-(p1m*p2p))/mll/sqrt(mll*mll+ptll*ptll);
}

// phi_CS, i.e. in the Collins-Soper frame
// Also can compute costheta_CS this way
double phi_CS(p4vec &lm, p4vec &lp){
    // Construct V
    p4vec V = lm + lp;
    double QV2 = V.m2();
    double QV = sqrt(QV2);
    double PTV = V.pT();
    double PTV2 = pow(PTV,2);
    double XTV = sqrt(QV2+PTV2);
    //
    double sign_pz = 1.0;
    if( V.pz() < 0 ) sign_pz *= -1.0;
    // 
    double pl_x = 
        +PTV/QV/XTV * ( -V.E()*lm.E() +V.pz()*lm.pz() )
        +XTV/QV/PTV * ( +V.px()*lm.px() +V.py()*lm.py() );
    double pl_y = sign_pz/PTV * ( -V.py()*lm.px() +V.px()*lm.py() );
    // double pl_z = sign_pz/XTV * (
    //      -V.pz()*lm.E()
    //      +V.E()*lm.pz() );
    // Now evaluate cos_theta_cs, phi_cs
    // double cos_th_cs = pl_z / sqrt( pow(pl_x,2) + pow(pl_y,2) + pow(pl_z,2) );
    // double th_cs = atan2( sqrt( pow(pl_x,2) + pow(pl_y,2)), pl_z );
    double phi_cs = atan2( pl_y, pl_x );
    if( phi_cs < 0 ) phi_cs += 2.*M_PI;
    return phi_cs;
}

// phistar variable
double phistar(p4vec &lm, p4vec &lp){
    double dphi = dphi_ij( lm, lp );
    double deta = lm.rap() - lp.rap(); // eta_lm - eta_lp, massless leptons
    double cos_th = tanh(deta/2.0); 
    double sin_th = sqrt( (1.+cos_th)*(1.-cos_th) ); // independent of sign of deta
    return tan( (M_PI-dphi)/2.0 ) * sin_th;
}

// mT variable, Charged Current (see 1701.07420)
double MT(p4vec &nu, p4vec &l){
    double dphi = dphi_ij( nu, l );
    return sqrt(2. * nu.pT() * l.pT() * ( 1. - cos(dphi) ));
}

// Could also pass an array by reference to a void function?
// Function for defining energies and momenta in a particle decay of the form mij -> mi + mj in the mij CoM frame
// Forms the basis for generating particle momenta using sequences of 2-body phase spaces
std::array<double,3> CoM_1to2_kinematics(double mij, double mi, double mj){
    // Define dimensionless variables
    double xi = max(0.,pow(mi/mij,2));
    double xj = max(0.,pow(mj/mij,2));
    // 3d array to contain Ei, Ej, |pi|
    std::array<double,3> output = {0.};
    // Ei = mij / 2. ( 1 + xi - xj )
    output[0] = mij / 2. * ( 1. + xi - xj );
    // Ej = mij / 2. ( 1 - xi + xj )
    output[1] = mij / 2. * ( 1. - xi + xj );
    // |pi| = |pj| = mij/2. * Kallen^{1/2}(1.,xi,xj)
    double vel = sqrt( kallen(1.,xi,xj) );    
    if( xi == 0.0 ){
        vel = ( 1. - xj );
    }
    if( xj == 0.0 ){
        vel = ( 1. - xi );
    }
    output[2] = mij / 2. * vel;    
    // output[2] = mij / 2. * sqrt( kallen(1.,xi,xj) );
    return output; 
}

std::array<double,3> CoM_1to2_kinematics_sq(double msq_ij, double msq_i, double msq_j){
    // Define dimensionless variables
    double xi = max(0.,msq_i/msq_ij);
    double xj = max(0.,msq_j/msq_ij);
    // 3d array to contain Ei, Ej, |pi|
    std::array<double,3> output = {0.};
    double mij = sqrt(msq_ij);
    // Ei = mij / 2. ( 1 + xi - xj )
    output[0] = mij / 2. * ( 1. + xi - xj );
    // Ej = mij / 2. ( 1 - xi + xj )
    output[1] = mij / 2. * ( 1. - xi + xj );
    // |pi| = |pj| = mij/2. * Kallen^{1/2}(1.,xi,xj)
    // sqrt( 1 -)
    double vel = sqrt( kallen(1.,xi,xj) );
    // What about the approximation?
    // if( xi < 1e-8 and xj < 1e-8 ){
    //     vel = 1. - xi - xj;
    // }
    if( xi == 0.0 ){
        vel = ( 1. - xj );
    }
    if( xj == 0.0 ){
        vel = ( 1. - xi );
    }
    output[2] = mij / 2. * vel;
    // Perform expansion?
    return output; 
}



////////////////////////////////////////////////////////////
// Information concerning the kinematic information class //
////////////////////////////////////////////////////////////


// The constructor for kinematic data
// Create the variable length objects accordingly
KinematicData::KinematicData(int npar) : npar(npar), pset(new p4vec[npar]), pdgs(new int[npar]) {
}

// Kinematic functions
// sij (currently didn't implement the SH functions)

// pi.pj
double KinematicData::pij(int i, int j) {
	if( i > npar or j > npar ) {
		cout << "KinematicData::pij, Error, p(" << max(i, j) << ") not initalised" << endl;
		cout << npar << endl;
		abort();
	}
	return dotp4(pset[i - 1], pset[j - 1]);
}

// pijext(int, p4vec)
double KinematicData::pijext(int i, p4vec j){
    if( i > npar ) {
        cout << "KinematicData::pijext, Error, p(" << i << ") not initalised" << endl;
        cout << npar << endl;
        abort();     
    }
    return dotp4(pset[i-1], j);
}

double KinematicData::sij(int i, int j) {
	if( i > npar or j > npar or min(i, j) < 1 ) {
		cout << "KinematicData::Za, p_i(" << i << ") not initalised" << endl;
		cout << "KinematicData::Za, p_j(" << j << ") not initalised" << endl;
		abort();
	}
    return 2. * pij(i,j);
}

// Return function to access the privately stored momentum information
p4vec KinematicData::p(int i) {
	if( i > npar ) {
		cout << "KinematicData::p, Error, p(" << i << ") not initalised" << endl;
		abort();
	}
	return pset[i - 1];
}

// Access Lab momenta (for cuts)
p4vec KinematicData::p_Lab(int i, double b_x, double b_y, double b_z) {
    if( i > npar ) {
        cout << "KinematicData::p, Error, p(" << i << ") not initalised" << endl;
        abort();
    }
    // Perform the boost to the Lab frame
    p4vec p_lab = pset[i - 1]; p_lab.boost(b_x,b_y,b_z);
    return p_lab;
}

// Function to activate storage of n random doubles [0,1]
void KinematicData::activate_ran( int nran ) {
    // Check if it was already activated
    if( r_i ){
        cout << "r_i has previously been activated" << endl;
    }

    // Otherwise create the nran length array of double
    r_i.reset( new double[nran] );
    // Also adjust the n_ran variable
    n_ran = nran;
    return;
}

// Access random variables (if set)
double KinematicData::r(int i) {
    if( i > n_ran ) {
        cout << "KinematicData::r, Error, r(" << i << ") not initalised" << endl;
        abort();
    }
    return r_i[i - 1];
}



