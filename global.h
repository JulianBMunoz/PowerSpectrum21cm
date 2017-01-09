////////////////////////////////////////////////////////////////////
//
//    Some utility constants
//
////////////////////////////////////////////////////////////////////


const int MAX_LINE=200;

////////////////////////////////////////////////////////////////////
//
//    Mathematical constants
//
////////////////////////////////////////////////////////////////////

const double PI=3.1415926535;

////////////////////////////////////////////////////////////////////
//
//    Conversion to natural units [ eV^# ]
//
////////////////////////////////////////////////////////////////////

const double METER=1./1.97e-7;
const double MPC=3.0857e22*METER;
const double SEC=1./6.58e-16;
const double KELVIN=1./1.16e4;
const double KILOGRAM=1./1.78e-36;

////////////////////////////////////////////////////////////////////
//
//    Physical constants
//
////////////////////////////////////////////////////////////////////

const double c=299792.458;   // speed of light [ km/s ]
const double alphaem=1./137;  // electromagnetic fine structure constant
const double echarge=sqrt(4.*PI*alphaem); // unit electric charge
const double SigT=0.6652e-28/pow(MPC/METER,2); // Thomson scattering cross section [ Mpc^2 ]
const double me=0.511e6*MPC;  // electron mass [ Mpc^-1 ]
const double mp=938.3e6*MPC;  // proton/hydrogen mass [ Mpc^-1 ]
const double mHe=mp/1.008*4.003; // helium mass [ Mpc^-1 ]

// energies in Mpc^-1
const double B1=13.6*MPC;    // hydrogen 1st ionization
const double chi0=23.72*MPC; // helium 1st ionization
const double chi1=52.5*MPC;  // helium 2nd ionization

// rates in Mpc^-1
const double LAM2s1s=8.227/SEC*MPC;  // hydrogen 2-photon decay rate from 2s

////////////////////////////////////////////////////////////////////
//
//    Model constants
//
////////////////////////////////////////////////////////////////////

const double Mpl=2.435e18*1e9/1.97e-7*3.0857e22;  // Planck mass ( Mpl^2 = 1/(8 PI G)) [ Mpc^-1 ]
const double h=0.678;    // reduced Hubble constant
const double H0=100*h/c; // Hubble constant [ Mpc^-1 ]
const double rhocri=3.*Mpl*Mpl*H0*H0;    // critical density at present time
const double T0=2.725/1.16e4/1.97e-7*3.0857e22;   // present CMB temperature [ Mpc^-1 ]
const double Omegac=0.227; // cold dark matter
const double Omegab=0.0456; // baryon
const double OmegaM=Omegac+Omegab;
const double Omegaga=PI*PI/15.*pow(T0,4)/3./Mpl/Mpl/H0/H0; // photons
const double Omeganu=21./8.*pow(4./11.,4./3.)*Omegaga; // neutrinos (3 flavors, left-handed, massless)
const double OmegaR=Omegaga+Omeganu;
const double OmegaL=1.-OmegaM-OmegaR; // dark energy
const double YHe=0.23;  // helium mass fraction
const double aCMB=1./1100; // scale factor at recombination

// primordial curvature power spectrum
const double DeltaZeta=exp(3.089)*1e-10;
const double ns=1-0.0397;
const double kpivot=0.05; // Planck pivot scale [ Mpc^-1 ]

// primordial tensor power spectrum
const double Tensor2Scalar=0.2;   // tensor-to-scalar ratio
const double nt=-Tensor2Scalar/8.;    // according to single-field slow-roll consistency relation
const double Kpivot=0.05; // sams as scalar pivot scale

// initial time
const double etaini=5e-4;

////////////////////////////////////////////////////////////////////
//
//    Computation of Wigner symbols
//
////////////////////////////////////////////////////////////////////

//const int LMAX=205;                    // maximum BiPoSH multipole
//const int lMAX=1010;                  // maximum CMB multipole
//double wig3j[LMAX][lMAX][LMAX];       // array for Wigner 3j symbols

////////////////////////////////////////////////////////////////////
//
//    Tabulation of spherical Bessel functions
//
////////////////////////////////////////////////////////////////////



const int lmax=1000;//15010;                  // maximum multipole
const int nls=1000;                  // how many ls we actually want
const double xmax=25000;               // maximum argument ( k_max * r_max )
const long jlsample=xmax/PI*10;               // no. of sampling (xmax/PI*10)
const double xinc=xmax/jlsample;      // increment for argument 
double sbesj[nls][jlsample+1][3];    // array for spherical Bessel functions
int llist[nls];			//list of the ls we want
const int ncentral=30;


// 	int llinear=500; //The cut in l where it starts increasing logarithmically//
// 	double lstep=(log(ltop)-log(llinear))/(npoints-llinear-1);
// 	
// 	
// 	for(j1=0;j1<llinear;++j1){
// 		larray[j1]=j1;	
// 	}
// 	
// 	//The first l we separate to avoid overcounting//
// 	
// 	for(j1=llinear;j1<npoints;++j1){
// 		larray[j1]=llinear*exp(lstep*(j1-llinear));
// //		printf("%d \n",larray[llinear+j1]);
// 	}

////////////////////////////////////////////////////////////////////
//
//    Cosmography
//
//////////////////////////////////////////////////////////////////// 

// range of scale factor (logarithmic scale)
const double amin=1e-9;
const double amax=1.5;//aCMB*2;
const int aNo=1000;
double eta_grid[3][aNo+1]; // grid to pre-compute conformal time eta

////////////////////////////////////////////////////////////////////
//
//    Recombination ( helium + hydrogen )
//
//////////////////////////////////////////////////////////////////// 

const double aHei=0.00015; // start of helium recombination
const double aHe2H=0.00058; // completion of helium recombination
const double aHf=0.00175; // completion of hydrogen recombination
const int aHeNo=400;   // sampling helium recombination
const int aHNo=400;    // sampling hydrogen recombination
double ne_grid[2][aHeNo+aHNo+1]; // grid for free electron number density
const int tauSampleNo=400;
const double amin_tau=5e-4;
const double amax_tau=1e-2;
double tau_grid[2][tauSampleNo+1]; // grid for photon optical depth

// photon/neutrino damping (smear out due to free streaming)
const double DAMPNU=0.1;
const double DAMPGA=0.5;

// BAO switch
const double ifBAO=1;

////////////////////////////////////////////////////////////////////
//
//    Reionization
//
//////////////////////////////////////////////////////////////////// 

const double zreion=10.5;  // reionization redshift
const double Delz=0.5;     // reionization width

const int tauSampleNo_re=400;
const double amin_tau_re=1./18;
const double amax_tau_re=0.9999;
double tau_re_grid[2][tauSampleNo_re+1]; // grid for reionization optical depth

////////////////////////////////////////////////////////////////////
//
//    (Linear) matter transfer function (at various source redshifts)
//
////////////////////////////////////////////////////////////////////

// source redshifts
const int zsNo=6;
const double zs[zsNo]={200,100,10,2,1,0.5};

// range of scalar wave number k [ Mpc^-1 ]
const double ksmin=0.001;
const double ksmax=1;
const int ksNo=200;
// tabulate transfer function
// 1st row is k
double Tdelta[zsNo+1][ksNo+1];

// data file for matter transfer functions
char matter_transfer_filename[]="matter_transfer.dat";

////////////////////////////////////////////////////////////////////
//
//    (Linear) tensor transfer function
//
////////////////////////////////////////////////////////////////////

const int lmax_ten=20;  // cut off the Boltzmann hierarchy 

// sample of tensor wave number K
const double Kmin=0.001;
const double Kmax=0.1;
const int KSampleNo=2400; // number of K's sampled

// CMB polarization angular modes
const int Jmin=2;
const int Jmax=20;
const int JNo=Jmax-Jmin+1;

//*************** recombination contribution **********************//

// range of conformal time for line of sight integral 
const double etai_ten=5e-4;
const double etaf_ten=1350;
const int etaNo_ten=500;

// tensor source function stored here (for each K)
// the 1st row is the conformal time
double source_ten[KSampleNo+1][etaNo_ten+1];

// transfer function for B-mode polarization
// the 1st row is the tensor wave number K
double TBgrid[JNo+1][KSampleNo];
double CBB_grid[JNo];   // polarization B-mode power spectrum

//*************** reionization contribution **********************//

// range of conformal time for line of sight integral 
const double etai_ten_re=5e-4;
const double etaf_ten_re=1.421714e4;
const int etaNo_ten_re=500;

// tensor source function stored here (for each K)
// the 1st row is the conformal time
double source_ten_re[KSampleNo+1][etaNo_ten+1];

// transfer function for B-mode polarization
// the 1st row is the tensor wave number K
double TBgrid_re[JNo+1][KSampleNo];
double CBB_re_grid[JNo];   // polarization B-mode power spectrum

////////////////////////////////////////////////////////////////////
//
//    Some toy cases
//
////////////////////////////////////////////////////////////////////

const double OmegaMRD=0.1;
const double eta0=1;

////////////////////////////////////////////////////////////////////
//
//    GW fossil signal
//
////////////////////////////////////////////////////////////////////

// sample number for wave number integration
const int kintegNo=10000;
const int KintegNo=10000;

// range for matter scales [Mpc^-1]
const double klow=ksmin;
const double kup=ksmax;

// range for tensor scales [Mpc^-1]
const double Klow=Kmin;
const double Kup=Kmax;

// matter angular modes
const int ellmin=200;
const int ellmax=1000;
const int ellNo=ellmax-ellmin+1;
// matter angular power spectra
double Cdd_grid[ellNo][zsNo]; // call the (ell - ellmin)th element

// auxillary functions: multipole, redshift, type (1 or 2)
double alpha_grid[ellNo][zsNo][2];
double beta_grid[JNo][zsNo][2];


