////////////////////////////////////////////////////////////////////////////////////////
//
//	Code to integrate two spherical bessel functions to get low ell Cells in dark ages
//	not necessarily l1=l2, to get off-diagonal components.
//
//
////////////////////////////////////////////////////////////////////////////////////////

const int velocityswitch=1; //Switch to activate(1) or deactivate (0) the relative velocity effect.
#define npoints 100 //number of ks we probe from 0 to ktop, logarithmically spaced.
#define nr 50 //number of rs we probe from 1/2 to 1 or from 1/r2-1 to 1 for r2 and r3 respectively.
//definitions for global.h to work.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h" /*Initialization and constants */
#include "auxiliar.h" /*definition of spherical bessel functions there */


#define nmax nls  //maximum low-ell w e use. Defined in global.h


int main(){

	FILE *fp;
	double z=30; // Redshift at which we calculate everything 
	double a=1./(z+1); //scale factor
	double Deltanu=1.0; //Bandwidth in MHz.
	
	
	
	int noiseswitch=0; // switch for noise. 0 is no noise. 1 is SKA. 2 is optimistic SKA (as defined in 1603.01206).
						
	int fkswitch=0; // switch for integrating with non standard fk. 0 is fk=1, 1 is k, 2 is k^2.
					// 3 is 1-k, 4 is (1-k)^2, 5 is 1/k and 6 is 1/k^2. All k in Mpc^(-1).
	
	int lowellswitch=1; //whether we include the lower ell. Code becomes slower.
	
	
	
	
	

//  Calculates the radial distance for a certain redshift.	
	
	double rcenter;	
	double rc0=15173.984582718711;
	double rc1=16855.594931908756;
	
	rcenter=rc0-rc1/sqrt(1.+z); //very good fit during Matter Domination with r0=15173.984582718711 and r1=16855.594931908756/





//we read a different file for the high-ell coefficients. Technically could be unified.
	long nelems=100; //Number of rows in the file.
	
	double zlist[nelems]; //Scale factors considered to integrate//
	double alphalist[nelems]; //Scale factors considered to integrate//
	double ttolist[nelems]; //Scale factors considered to integrate//
	
	double alpha, tto; //dT21/ddelta, dT21/dv .
	
	fp=fopen("coeffs_21cm.dat","r");	//File with {z, T21bar, alpha, beta, gamma} (tto=T21bar) from z=200 to z=20 in DESCENDING order.
	
	 //We will use ja to count the number of elements//
	
	int ja;
	
	double tempa; //useess temp variable.
	
	for(ja=0;fscanf(fp,"%le %le %le %le %le",&zlist[ja],&ttolist[ja],&alphalist[ja],&tempa,&tempa)==5;++ja){
 	}

	fclose(fp);	
	
		
	reverse(alphalist,nelems); 
	reverse(ttolist,nelems);
	reverse(zlist,nelems);//to make it ascend in order
				

	alpha =interpol(alphalist,zlist,nelems,z)*1000.;	//in mK	
	tto=interpol(ttolist,zlist,nelems,z)*1000.;	//
	
	printf("alpha=%le, tto=%le \n",alpha, tto);










//We first assign to the nls elements of llist the values of l that we want from 0 to 110 //



	extern int llist[nls];				//list of the ls we will calculate


	int g,g2,j;

	for(g=0;g<nls;g+=1){
		llist[g]=g;	
	}
	
	long j1,j2,i2;






	
	
	double fpre=2./PI*DeltaZeta*2.*PI*PI*pow(kpivot,1-ns); //z-independent prefactor	
	
	double fprez;//z-dependent prefactor, growth^2
	


//Now we get the 21-cm coefficients during the dark ages

	
	
	
	
//We get the growth factor.
//1/N spaced.
	int length2=100; //Number of elements growth function output
	double atab2[length2],Dptab[length2]; //from z=10 to 1000.
	
	fp=fopen("growth.dat","r"); //file with (a, Dp)
	
	for(j=0;fscanf(fp,"%le %le ",&atab2[j],&Dptab[j])==2;++j)
		;
	
	fclose(fp);		
	
	double Dp;

	
	Dp=interpol(Dptab, atab2, length2, a);	
	



	fprez=Dp*Dp; //T in mK^2 def before.



//First we do the tabulation of r.


	double datapoints=10000;//20000.; //number of k-points we integrate over.
	int rdatapoints=200; //how many r-points we integrate over.
	
	
	double *rintegrand; //internal integrand
	rintegrand=create_1D_array(rdatapoints);
	
	double *W1;// Window function.
	W1=create_1D_array(rdatapoints);
	double normW=0.;//normalization of window

	
	double r;
	double Deltar= Deltanu*sqrt((1.+z)/51.)*60.;//to convert from dnu to dr. dnu in MHz, dr in Mpc.
	int Nsigma=4; //number of sigmas away from center of bandwidth we integrate over.
	
	double rmin=rcenter-Nsigma*Deltar;
	double rmax=rcenter+Nsigma*Deltar;
	
//	printf("%le %le %le \n",rmin,rmax,rcenter);
	
	double *rinterpol;
	rinterpol=create_1D_array(rdatapoints);

	double dr=(rmax-rmin)/(rdatapoints-1.);
//	printf("%le \n",dr);
	
	
	for(j1=0;j1<rdatapoints;j1++){
		rinterpol[j1]=rmin+(rmax-rmin)/(rdatapoints-1.)*j1;
		W1[j1]=exp(-(rinterpol[j1]-rcenter)*(rinterpol[j1]-rcenter)/Deltar/Deltar/2.);
		normW+=W1[j1]; //so it integrates to 1.
//		printf("%ld, %le %le %le \n", j1, rinterpol[j1],W1[j1],normW);
	}
	

	double checksum=0;	
 	for(j1=0;j1<rdatapoints;j1++){
 		W1[j1]=W1[j1]/normW;
 		checksum+=W1[j1];
 	}
// 	printf("%le \n",checksum);






	
	
//Reads the Transfer function from CAMB and does an interpolation as a function of x.
	
	
	char filename[200];
	int lengthname;
		
				
	fp=fopen("cambz0.dat","r");
	
	int length=250; //Number of elements in the transfer function output, usually 302, but cut off to be k/h<0.3 Mpc^-1.
	double TF[length],koverh[length]; /*The baryon transfer function and k/h obtained
	from Lambda CAMB code */
	
	

	double temp[length];
	
	for(j=0;j<length;++j){
		fscanf(fp,"%le %le %le %le %le %le %le ",&koverh[j],temp ,&TF[j],temp,temp,temp,temp);
	}
//   	for(j=0;j<length;++j)
//   		printf("k/h=%le TF=%le \n",koverh[j],TF[j]);
	
	fclose(fp);

	
//	We need to give the right values to TF, evaluated over k and not k/h, we define kgrid for that //
	
	double kgrid[length];
	
	for(j=0;j<length;++j)
		kgrid[j] = h * koverh[j];
	
	double kmin=kgrid[0];
	double kmax=kgrid[length-1];
		

	
	double *kinterpol;
	kinterpol=create_1D_array((int)datapoints);	

//	This part does the interpolation. //
	
	double **fk; //Different than standard (f=1) function.
	int Nfk=7; //number of different functions
	fk=create_2D_array(Nfk,(int)datapoints);	
	
	
	
	
	for (j=0; j<datapoints;++j){
		kinterpol[j]=kmin*exp((log(kmax)-log(kmin))*j/(datapoints-1.));
//	 	printf("j= %d, k=%le \n",j,kinterpol[j]);
		fk[0][j]=1.;
		fk[1][j]=kinterpol[j];
		fk[2][j]=kinterpol[j]*kinterpol[j];
		fk[3][j]=1.-kinterpol[j];
		fk[4][j]=(1.-kinterpol[j])*(1.-kinterpol[j]);
		fk[5][j]=1./kinterpol[j];
		fk[6][j]=1./kinterpol[j]/kinterpol[j];
	}





	double bes1,*integrand;
	integrand=create_1D_array((int)datapoints);
	
	
	
	
//	We tabulate the spherical bessel functions up to l=lmax. //		
 	int sbesselj_ini(double sbesselj[nls][jlsample+1][3]);	
	if(lowellswitch==1){	
 		sbesselj_ini(sbesj);
 	}
 	//Remember to initiate sbesj variable in global.h




	  	
//Now we integrate dr W(r)jl(kx)			//
//		as a function of k					//
//											//
//											//
//											//


	


	double ki, r1,r2; //kintegration  and rintegration
	double dk; //differential of k
	double dlogk=log(kinterpol[3]/kinterpol[2]);

	double power; //matter power spectrum

	
	double tf;// transfer function
	
		

	
	
	int lpoints=60; //number of ls we integrate over.
	long llistint[lpoints]; //list of the ls we integrate over.		
	int lcutoff=lpoints/2; //where we go from linear to log spacing.
	double lstep=exp((log((nls-2.)/lcutoff))/(1.*(lpoints-lcutoff-1.))); //logstep from end of linear to final

	long gg;
	
	if(nls<=lpoints){
		for(gg=0;gg<nls;gg+=1){
			llistint[gg]=gg;	
		}
	}
	else{
		for(gg=0;gg<lcutoff;gg+=1){
			llistint[gg]=gg;	
		}
		for(gg=lcutoff;gg<lpoints;gg+=1){
			llistint[gg]=(long)(llistint[lcutoff-1]*pow(lstep,1.*(gg-lcutoff))); //
		}
	}
	
	
for(gg=1;gg<lpoints;gg+=1){ //to make sure they are in ascending order and not repeated
		if(llistint[gg]<=llistint[gg-1])
		{	
			llistint[gg]=llistint[gg-1]+1;
		}
//		printf(" %ld \n", llistint[gg]);
	}
	
 

	double **tl,**tl1,**tl2; //21-cm transfer function (l,k) at fixed z, for Deltal=0, 1 and 2.
	double **tl3,**tl4;
	tl=create_2D_array(lpoints,(int)datapoints);
	tl1=create_2D_array(lpoints,(int)datapoints);
	tl2=create_2D_array(lpoints,(int)datapoints);
	tl3=create_2D_array(lpoints,(int)datapoints);
	tl4=create_2D_array(lpoints,(int)datapoints);	
		
		
	double *Cl0,*Cl1,*Cl2; //Atual Cells we want to measure, 0, +-1 or +-2 away. For low ell
	double *Cl3,*Cl4;
	Cl0 = (double*)calloc(lpoints, sizeof(double));
	Cl1 = (double*)calloc(lpoints, sizeof(double));
	Cl2 = (double*)calloc(lpoints, sizeof(double));
	Cl3 = (double*)calloc(lpoints, sizeof(double));
	Cl4 = (double*)calloc(lpoints, sizeof(double));

		
	long l1;
	double a1[3];
				
if(lowellswitch==1){				
	for (j1=2;j1<lpoints;j1++){
		l1=llistint[j1];
		a1[0]=-l1*(l1-1.)/(4.*l1*l1-1)*tto*velocityswitch;
		a1[1]=(2.*l1*l1+2.*l1-1)/(4.*l1*l1+4*l1-3)*tto+alpha;
		a1[2]=-(l1+2.)*(l1+1)/(2.*l1+1)/(2.*l1+3)*tto*velocityswitch;
//		printf("%ld, %ld , %le , %le , %le \n",j1, l1,a1[0],a1[1],a1[2]);	
		for (i2=0;i2<datapoints;i2++){
			ki=kinterpol[i2];	
			tl[j1][i2]=0;
			for(g=0;g<rdatapoints;g++){
				r1=rinterpol[g];
				bes1=a1[0]*get_sbesj(sbesj,l1-2,ki*r1)+a1[1]*get_sbesj(sbesj,l1,ki*r1)+a1[2]*get_sbesj(sbesj,l1+2,ki*r1);
				rintegrand[g]=bes1*W1[g]; // no dr since we modified W1 to integrate to 1 w/o it.
//				printf("%d, %le, %le \n",l1,get_sbesj(sbesj,l1-2,ki*r1),rintegrand[g]);
				tl[j1][i2]+=rintegrand[g];	
				}
		}
//		printf("%ld, %ld , %le  \n",j1, l1,tl[j1][(int)datapoints-2]);	
		
		l1+=1;
		a1[0]=-l1*(l1-1.)/(4.*l1*l1-1)*tto*velocityswitch;
		a1[1]=(2.*l1*l1+2.*l1-1)/(4.*l1*l1+4*l1-3)*tto+alpha;
		a1[2]=-(l1+2.)*(l1+1)/(2.*l1+1)/(2.*l1+3)*tto*velocityswitch;
//		printf("%ld, %ld \n",j1, l1);	
		for (i2=0;i2<datapoints;i2++){
			ki=kinterpol[i2];	
			tl1[j1][i2]=0;
			for(g=0;g<rdatapoints;g++){
				r1=rinterpol[g];
				bes1=a1[0]*get_sbesj(sbesj,l1-2,ki*r1)+a1[1]*get_sbesj(sbesj,l1,ki*r1)+a1[2]*get_sbesj(sbesj,l1+2,ki*r1);
				rintegrand[g]=bes1*W1[g]; // no dr since we modified W1 to integrate to 1 w/o it.
				tl1[j1][i2]+=rintegrand[g];	
				}
		}
//		printf("%ld, %le \n",j1,tl[j1][(int)datapoints-3]);

		
		l1+=1;
		a1[0]=-l1*(l1-1.)/(4.*l1*l1-1)*tto*velocityswitch;
		a1[1]=(2.*l1*l1+2.*l1-1)/(4.*l1*l1+4*l1-3)*tto+alpha;
		a1[2]=-(l1+2.)*(l1+1)/(2.*l1+1)/(2.*l1+3)*tto*velocityswitch;
//		printf("%ld, %ld \n",j1, l1);	
		for (i2=0;i2<datapoints;i2++){
			ki=kinterpol[i2];	
			tl2[j1][i2]=0;
			for(g=0;g<rdatapoints;g++){
				r1=rinterpol[g];
				bes1=a1[0]*get_sbesj(sbesj,l1-2,ki*r1)+a1[1]*get_sbesj(sbesj,l1,ki*r1)+a1[2]*get_sbesj(sbesj,l1+2,ki*r1);
				rintegrand[g]=bes1*W1[g]; // no dr since we modified W1 to integrate to 1 w/o it.
				tl2[j1][i2]+=rintegrand[g];	
				}
		}
		
		
		l1+=1;
		a1[0]=-l1*(l1-1.)/(4.*l1*l1-1)*tto*velocityswitch;
		a1[1]=(2.*l1*l1+2.*l1-1)/(4.*l1*l1+4*l1-3)*tto+alpha;
		a1[2]=-(l1+2.)*(l1+1)/(2.*l1+1)/(2.*l1+3)*tto*velocityswitch;
//		printf("%ld, %ld \n",j1, l1);	
		for (i2=0;i2<datapoints;i2++){
			ki=kinterpol[i2];	
			tl3[j1][i2]=0;
			for(g=0;g<rdatapoints;g++){
				r1=rinterpol[g];
				bes1=a1[0]*get_sbesj(sbesj,l1-2,ki*r1)+a1[1]*get_sbesj(sbesj,l1,ki*r1)+a1[2]*get_sbesj(sbesj,l1+2,ki*r1);
				rintegrand[g]=bes1*W1[g]; // no dr since we modified W1 to integrate to 1 w/o it.
				tl3[j1][i2]+=rintegrand[g];	
				}
		}
		
		
		
		l1+=1;
		a1[0]=-l1*(l1-1.)/(4.*l1*l1-1)*tto*velocityswitch;
		a1[1]=(2.*l1*l1+2.*l1-1)/(4.*l1*l1+4*l1-3)*tto+alpha;
		a1[2]=-(l1+2.)*(l1+1)/(2.*l1+1)/(2.*l1+3)*tto*velocityswitch;
//		printf("%ld, %ld \n",j1, l1);	
		for (i2=0;i2<datapoints;i2++){
			ki=kinterpol[i2];	
			tl4[j1][i2]=0;
			for(g=0;g<rdatapoints;g++){
				r1=rinterpol[g];
				bes1=a1[0]*get_sbesj(sbesj,l1-2,ki*r1)+a1[1]*get_sbesj(sbesj,l1,ki*r1)+a1[2]*get_sbesj(sbesj,l1+2,ki*r1);
				rintegrand[g]=bes1*W1[g]; // no dr since we modified W1 to integrate to 1 w/o it.
				tl4[j1][i2]+=rintegrand[g];	
				}
		}
		
	}
}
 

	
	
	
	
	
	
	double fkk=0.; //evaluation of fk(k)
	
	

if(lowellswitch==1){		
	for (j1=2;j1<lpoints;j1++){
		Cl0[j1]=0.;
		Cl1[j1]=0.;
		Cl2[j1]=0.;
		Cl3[j1]=0.;
		Cl4[j1]=0.;
		for (i2=0;i2<datapoints;i2++){
			ki=kinterpol[i2];
			dk=dlogk*ki;
			fkk=fk[fkswitch][i2]; //non-standard is fk=1.
			tf=interpol(TF, kgrid, length,ki)*ki*ki;
			power=1.*pow(ki,-4+ns)*tf*tf; //Power spectum with unity amplitude	
			integrand[i2]=ki*ki*tl[j1][i2]*tl[j1+0][i2]*power*dk*fkk;		
			Cl0[j1]+=integrand[i2];
			integrand[i2]=ki*ki*tl[j1][i2]*tl1[j1][i2]*power*dk*fkk;		
			Cl1[j1]+=integrand[i2];
			integrand[i2]=ki*ki*tl[j1][i2]*(tl2[j1][i2])*power*dk*fkk;		
			Cl2[j1]+=integrand[i2];
			integrand[i2]=ki*ki*tl[j1][i2]*(tl3[j1][i2])*power*dk*fkk;		
			Cl3[j1]+=integrand[i2];
			integrand[i2]=ki*ki*tl[j1][i2]*(tl4[j1][i2])*power*dk*fkk;		
			Cl4[j1]+=integrand[i2];
		}
		Cl0[j1]*=fpre*fprez;
		Cl1[j1]*=fpre*fprez;
		Cl2[j1]*=fpre*fprez;
		Cl3[j1]*=fpre*fprez;
		Cl4[j1]*=fpre*fprez;		
//		printf("%ld, %le \n",j1,Cl0[j1]);
	}
}	
	
	double noise=0.; //CLNoise, constant in ell.
	
	double Dbase; //in km
	double fcover;
	double Nyears; 
	double t0; //in MHz^-1
	double fsky;
	
	double lcover; 
	double Tsys= 180.*pow(45./180.*31./(1.+z),-2.6)*1000.;	 // in mK
	
	if(noiseswitch==1){//SKA-like noise
		fsky=0.75;
		Dbase=6.; //in km
		fcover=0.02;
		Nyears=5.;
		t0=Nyears*365.25*86400.*pow(10.,6.); //T of observation in MHz^-1
		lcover= (2*PI)*Dbase*1000./(0.21*(1.+z));
		noise=pow(2*PI,3.)*Tsys*Tsys/(Deltanu*t0)*1./fcover/fcover/lcover/lcover;
	}
	if(noiseswitch==2){//futuristic (earth-based) noise
		fsky=0.75;	
		Dbase=100.; //in km
		fcover=0.2;
		Nyears=10.;
		t0=Nyears*365.25*86400.*pow(10.,6.); //T of observation in MHz^-1
		lcover= (2*PI)*Dbase*1000./(0.21*(1.+z));
		noise=pow(2*PI,3.)*Tsys*Tsys/(Deltanu*t0)*1./fcover/fcover/lcover/lcover;
	}
	
	
	for (j1=2;j1<lpoints;j1++){
		Cl0[j1]+=noise;
		Cl1[j1]+=noise;
		Cl2[j1]+=noise;
		Cl3[j1]+=noise;
		Cl4[j1]+=noise;
	}
	
	
	
	
	
// 	 fp=fopen("lowellCl.dat","w");
// 
//  	 for(j1=2;j1<lpoints;++j1){
//  	 	fprintf(fp,"%ld \t\t\t %le \n",llistint[j1], Cl0[j1]); //in mK^2
// 	 }
// 
// 	fclose(fp);
//UNCOMMENT IF YOU WANT TO SAVE THE LOW ELL ONLY.
	
	
	
	
	//Now we calculate high-ell in flat-sky approximation//	
	//													 //
	//													 //
	//													 //
	//													 //
	//													 //



	
		
	fp=fopen("cambz0.dat","r");
	
	int length3=300; //Number of elements in the transfer function output
	double TF3[length3],koverh3[length3]; //The baryon transfer function and k/h obtained from Lambda CAMB code
		

	double temp3[length3];
	
	for(j=0;fscanf(fp,"%le %le %le %le %le %le %le ",&koverh3[j],temp3 ,&TF3[j],temp3,temp3,temp3,temp3)==7;++j)
		;
//   	for(j=0;j<length3;++j)
//   		printf("k/h=%le TF=%le \n",koverh[j],TF[j]);
	
	fclose(fp);

	
//	We need to give the right values to TF, evaluated over k and not k/h, we define kgrid for that //
	
	double kgrid3[length3];
	
	for(j=0;j<length3;++j)
		kgrid3[j] = h * koverh3[j];
	
	double kmin3=kgrid3[0];
	double kmax3=kgrid3[length3-1];
	
	
	 		
	
	


	double datapoints3=120.; //Number of datapoints that we will integrate over. 100 seems enough.

		
	
	double fk3[Nfk]; //Different than standard (f=1) function. For high-ell now	
	
	
//	long npoints=100; DEFINED IN THE HEADER.
	

	double *Pow; //Power spectrum delta delta
	double *Powv; //Power spectrum delta v
	double *Powvv; //Power spectrum v v
	
	Pow=(double*)calloc(npoints, sizeof(double));	
	Powv=(double*)calloc(npoints, sizeof(double));	
	Powvv=(double*)calloc(npoints, sizeof(double));	
	
	
	
	double *karray,*karrayz;
	karray= (double*)calloc(npoints, sizeof(double));	
	karrayz= (double*)calloc(datapoints3, sizeof(double));	

	double ktop=10; //Maximum k [Mpc-1], corresponds to l~5*10^5.
	double kbot=llistint[lpoints-1]/rcenter; //Minimum k [Mpc-1], below k=0.001 the N would be too small for our approximations.

	double logstep=log(ktop/kbot)/(npoints-1.);

	for(j1=0;j1<npoints;j1++){
		karray[j1]=kbot*exp(j1*logstep);//The k list for the orthogonal directions (radially)
//		printf("%le \n",karray[j1]);
	}
	

	


	printf("Max k perpendicular=%le, if bigger than %le danger!! \n",ktop,kmax3/2); //Should always be fine.



	double *Cl0h,*Cl1h,*Cl2h,*lhigh; //Atual Cells and ells we want to measure, 0, +-1 or +-2 away. For hgh ell

	lhigh= (double*)calloc(npoints, sizeof(double));
	Cl0h = (double*)calloc(npoints, sizeof(double));
	Cl1h = (double*)calloc(npoints, sizeof(double));
	Cl2h = (double*)calloc(npoints, sizeof(double));



    
   
	for(j1=0;j1<datapoints3;j1++){
		karrayz[j1]=kmin3*exp(1.*j1/datapoints3*(log(kmax3/3)-log(kmin3))); //The k list for the z direction (positive and negative)
		//The maximum k is kmax/3 so that we do not overflow the interpol(x) function.
	}		
   
	int n1,n1z;
   
	double k1,tf1,dk1z;
   
   
	double k1z;
	
	double win1; //window function in k-space.
	double sigma=60.*sqrt(1.+z)/sqrt(51.)*Deltanu; //width in k of the window function.
	   
    	   
	   
	   
	   
//We calculate the Power Spectrum.


	for (n1=0;n1<npoints;n1++){	
		Pow[n1]=0.;
		for (n1z=1;n1z<datapoints3;n1z++){
			k1=sqrt(karray[n1]*karray[n1]+karrayz[n1z]*karrayz[n1z]);
			fk3[0]=1.;
			fk3[1]=k1;
			fk3[2]=k1*k1;
			fk3[3]=1.-k1;
			fk3[4]=(1.-k1)*(1.-k1);
			fk3[5]=1./k1;
			fk3[6]=1./k1/k1;
			fkk=fk3[fkswitch]; //non-standard is fk=1.
			win1=exp(-karrayz[n1z]*karrayz[n1z]*sigma*sigma); //it's actually the window squared, so it does not have 1/2. in exponent.
			tf1=interpol(TF3, kgrid3, length3,k1)*k1*k1;
			dk1z=(karrayz[n1z]-karrayz[n1z-1])*win1*fkk; 
			Pow[n1]+=(tf1*tf1/pow(k1,4-ns))*dk1z; //Not including pow(L,-2)*pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*
			Powv[n1]+=(tf1*tf1/pow(k1,4-ns))*dk1z *karrayz[n1z]*karrayz[n1z]/k1/k1 ; //Not including pow(L,-2)*pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*
			Powvv[n1]+=(tf1*tf1/pow(k1,4-ns))*dk1z *karrayz[n1z]*karrayz[n1z]/k1/k1 *karrayz[n1z]*karrayz[n1z]/k1/k1; //Not including pow(L,-2)*pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*
		} 
		Pow[n1]*=pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*2;//Last *2 to account for kz>0 and <0.	 //pow(L,-2)*	 
		Powv[n1]*=pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*2;//Last *2 to account for kz>0 and <0.	 //pow(L,-2)*	 
		Powvv[n1]*=pow(2*PI,-1)*(2*PI*PI*DeltaZeta*pow(kpivot,1-ns))*2;//Last *2 to account for kz>0 and <0.	 //pow(L,-2)*	 
	 } 	






	
// 	lengthname=sprintf(filename,"PowDiscA-%.1f.dat",Deltanu); //We reuse the same filename variable name. A for anisotropies
// 	fp=fopen(filename,"w");		
// 
// 
// 	 for(j1=0;j1<npoints;++j1){
// 			fprintf(fp,"%le \t\t\t %le \n",karray[j1],Pow[j1]);
// 		}
// 	fclose(fp);	
// 	
// 	lengthname=sprintf(filename,"PowDiscAv-%.1f.dat",Deltanu); //We reuse the same filename variable name.
// 	fp=fopen(filename,"w");	
// 
// 
// 	 for(j1=0;j1<npoints;++j1){
// 			fprintf(fp,"%le \t\t\t %le \n",karray[j1],Powv[j1]);
// 		}
// 	fclose(fp);	
// 	
// 	lengthname=sprintf(filename,"PowDiscAvv-%.1f.dat",Deltanu); //We reuse the same filename variable name.
// 	fp=fopen(filename,"w");	
// 
// 
// 	 for(j1=0;j1<npoints;++j1){
// 			fprintf(fp,"%le \t\t\t %le \n",karray[j1],Powvv[j1]);
// 		}
// 	fclose(fp);			
//UNCOMMENT IF YOU WANNA SAVE THE POWER SPECTRUM INSTEAF OF CELLS.

















	for (n1=0;n1<npoints;n1++){	
		Cl0h[n1]=(alpha*alpha*Pow[n1]+
		2.*alpha*tto*Powv[n1]+ //Note we need the factor of 2, since we did not calculate both permutations (<delta v> and <v delta>). This isn't true for Bispectrum.
		tto*tto*Powvv[n1]
		)*pow(Dp,2)/rcenter/rcenter;
		Cl1h[n1]=Cl0h[n1]; //it is the same than ll
		Cl2h[n1]=Cl0h[n1]; // it is the same too
		lhigh[n1]=rcenter*karray[n1]; 
	 } 
	 
	 
	 
	 //We add noise now, which we defined a while ago. No need to redefine here.
	 for(n1=0;n1<npoints;n1++){
		Cl0h[n1]+=noise;
		Cl1h[n1]+=noise;
		Cl2h[n1]+=noise;
	}





	//Now we combine high and low-ell //	
	//								  //
	//								  //
	//								  //
	//								  //



	//int ntotal=nls+npoints;
	

  	lengthname=sprintf(filename,"Cllz%.0f-bw%.1f-n%d-fk%d.dat",z,Deltanu,noiseswitch,fkswitch); //We reuse the same variable name.

	fp=fopen(filename,"w"); //file with (l, Cell) for whatever noise, fk, z and bandwidth we specified.
	
	for(j=4;j<lpoints;j++){
		fprintf(fp,"%ld %le \n",llistint[j],Cl0[j]);
	}
	
	for(n1=0;n1<npoints;n1++){
		fprintf(fp,"%le %le \n",lhigh[n1],Cl0h[n1]);
	}
	
	fclose(fp);	
	
	  
	  
	lengthname=sprintf(filename,"Cll+1z%.0f-bw%.1f-n%d-fk%d.dat",z,Deltanu,noiseswitch,fkswitch); //We reuse the same variable name.

	fp=fopen(filename,"w"); //file with (l, Cell) for whatever noise, fk, z and bandwidth we specified.
	
	for(j=4;j<lpoints;j++){
		fprintf(fp,"%ld %le \n",llistint[j],Cl1[j]);
	}
	
	for(n1=0;n1<npoints;n1++){
		fprintf(fp,"%le %le \n",lhigh[n1],Cl1h[n1]);
	}
	
	fclose(fp);	
	
	lengthname=sprintf(filename,"Cll+2z%.0f-bw%.1f-n%d-fk%d.dat",z,Deltanu,noiseswitch,fkswitch); //We reuse the same variable name.

	fp=fopen(filename,"w"); //file with (l, Cell) for whatever noise, fk, z and bandwidth we specified.
	
	for(j=4;j<lpoints;j++){
		fprintf(fp,"%ld %le \n",llistint[j],Cl2[j]);
	}
	
	for(n1=0;n1<npoints;n1++){
		fprintf(fp,"%le %le \n",lhigh[n1],Cl2h[n1]);
	}
	
	fclose(fp);	
	
	
	lengthname=sprintf(filename,"Cll+3z%.0f-bw%.1f-n%d-fk%d.dat",z,Deltanu,noiseswitch,fkswitch); //We reuse the same variable name.

	fp=fopen(filename,"w"); //file with (l, Cell) for whatever noise, fk, z and bandwidth we specified.
	
	for(j=4;j<lpoints;j++){
		fprintf(fp,"%ld %le \n",llistint[j],Cl3[j]);
	}
	
	for(n1=0;n1<npoints;n1++){
		fprintf(fp,"%le %le \n",lhigh[n1],Cl0h[n1]);
	}
	
	fclose(fp);	
	
	lengthname=sprintf(filename,"Cll+4z%.0f-bw%.1f-n%d-fk%d.dat",z,Deltanu,noiseswitch,fkswitch); //We reuse the same variable name.

	fp=fopen(filename,"w"); //file with (l, Cell) for whatever noise, fk, z and bandwidth we specified.
	
	for(j=4;j<lpoints;j++){
		fprintf(fp,"%ld %le \n",llistint[j],Cl4[j]);
	}
	
	for(n1=0;n1<npoints;n1++){
		fprintf(fp,"%le %le \n",lhigh[n1],Cl0h[n1]);
	}
	
	fclose(fp);		




	
	
	
	//Finally, we free memory.
	
		
	free(Cl0);
	free(Cl1);
 	free(Cl2);
	free(Cl3);
 	free(Cl4); 	
 			
	free(Cl0h);
	free(Cl1h);
 	free(Cl2h);
	free(lhigh);		
		
	free(W1);
	free(rinterpol);
	free(rintegrand);
	
	
	free(Pow);
	free(Powv);
	free(Powvv);
	
	free(kinterpol);
	free(integrand);
  	
	free_2D_array(fk,Nfk);	
	free_2D_array(tl, lpoints);
	free_2D_array(tl1, lpoints);
	free_2D_array(tl2, lpoints);
	free_2D_array(tl3, lpoints);
	free_2D_array(tl4, lpoints);	
	
	
	free(karray);
	free(karrayz);	
	
	
	
}






