double *create_1D_array(unsigned n1){
   double *matrix = (double *) calloc(n1, sizeof(double));
   if (matrix == NULL) {
      fprintf(stderr, "Error in create_1D_array: unable to allocate memory\n");
      exit(1);
    }
   return matrix;
}



double **create_2D_array(unsigned n1, unsigned n2){
   unsigned i;
   double **matrix = (double **) calloc(n1, sizeof(double *));
   if (matrix == NULL){
      fprintf(stderr, "Error in create_2D_array: unable to allocate memory\n");
      exit(1);
    }
   for (i = 0; i < n1; i++)  matrix[i] = create_1D_array(n2);
   
   return matrix;
}

void free_2D_array(double **matrix, unsigned n1){
    unsigned i;
    for (i = 0; i < n1; i++) free(matrix[i]);
    free(matrix);
}


double ***create_3D_array(unsigned n1, unsigned n2, unsigned n3){
   unsigned i;
   double ***matrix = (double ***) calloc(n1, sizeof(double **));
   if (matrix == NULL){
      fprintf(stderr, "Error in create_3D_array: unable to allocate memory\n");
      exit(1);
    }
   for (i = 0; i < n1; i++)  matrix[i] = create_2D_array(n2,n3);
   
   return matrix;
}


void free_3D_array(double ***matrix, unsigned n1,unsigned n2){
    unsigned i;
    
    for (i = 0; i < n1; i++){    
    	free_2D_array(matrix[i],n2);
    }
    
    free(matrix);
}



//////////////////////////////////////////////////////////////////////////////////////
//
//    This code tabulates spherical Bessel functions of the first kind:  j_l (x)
//                                                and first derivative:  dj_l/dx (x)
//                                               and second derivative:  d^2j_l/dx^2 (x)
//
//    truncation: 0 <= l <= lmax
//                0 <= x <= xmax
//          step: xinc
//      sampling: jlsample=1+xmax/xinc
//
//    data storage:  sbesselj[l][i][0] == j_l (i*xinc)
//                   sbesselj[l][i][1] == dj_l/dx (i*xinc)
//                   sbesselj[l][i][2] == d^2 j_l/dx ^2 (i*xinc)
//FROM LIANG
///////////////////////////////////////////////////////////////////////////////////////


int sbesselj_ini(double sbesselj[][jlsample+1][3])
{
    long double x;
    long i;
    long l;
    double d=5e-8;
    int j; //Positions of the ls we will store. nls of them.
    
    printf("Start jl tabulation: lmax=%d, xmax=%f, xinc=%f\n",lmax,xmax,xinc);
  
    sbesselj[0][0][0]=1.;
    sbesselj[0][0][1]=0.;
    sbesselj[0][0][2]=-1./3.;
    sbesselj[1][0][0]=0.;
    sbesselj[1][0][1]=1./3.;
    sbesselj[1][0][2]=0.;
    
    double *ptr1,*ptr2,*ptr3;
    double *d1,*d2,*d3;
    double *dd1,*dd2,*dd3;
    
    ptr1 = (double*)calloc(jlsample, sizeof(double)); //Pointers for the values of Bessel functions
    ptr2 = (double*)calloc(jlsample, sizeof(double));
    ptr3 = (double*)calloc(jlsample, sizeof(double));
    d1 = (double*)calloc(jlsample, sizeof(double));
    d2 = (double*)calloc(jlsample, sizeof(double));
    d3 = (double*)calloc(jlsample, sizeof(double));
    dd1 = (double*)calloc(jlsample, sizeof(double));
    dd2 = (double*)calloc(jlsample, sizeof(double));
    dd3 = (double*)calloc(jlsample, sizeof(double));

    
    //Starting at l=0
    
     ptr1[0]=1.,ptr2[0]=0.,ptr3[0]=-1./3;
     d1[0]=0.,d2[0]=1./2,d3[0]=0.;
     dd1[0]=0.,dd2[0]=0.,dd3[0]=2./15;
    
 
     for(l=2; l<nls; l++) {
        sbesselj[l][0][0]=0.;
    	sbesselj[l][0][1]=0.;
        sbesselj[l][0][2]=(l==2)?(2./15.):0.;
    }
    
    for(i=1; i<jlsample; i++) {
             
        x=i*xinc;
        ptr1[i]=sbesselj[0][i][0]=sin(x)/x;
		ptr2[i]=sbesselj[1][i][0]=sin(x)/x/x-cos(x)/x;
		d1[i]= sbesselj[0][i][1]=-sbesselj[1][i][0];
		d2[i]= sbesselj[1][i][1]=sbesselj[0][i][0]-2*sbesselj[1][i][0]/x;
		dd1[i]= sbesselj[0][i][2]=-2.*cos(x)/x/x+sin(x)/x*(2./x/x-1.);
		dd2[i]= sbesselj[1][i][2]=sbesselj[0][i][1]/3.+sbesselj[1][i][0]*(6./x/x-2./3.)-2./x*sbesselj[0][i][0];
       
        for(l=j=2; l<=2*x;) {
            
 //           if((fabs(ptr2[i])<d)&&(x<l)) {l++; j++; break;} //This is a precautionary measure, to make sure there are no weird bumps.

        
        	ptr3[i] = (2*l-1)*ptr2[i]/x-ptr1[i];
        	d3[i]=ptr2[i]-(l+1)*ptr3[i]/x;
        	dd3[i]=d2[i]*l/(2*l+1.)+ptr3[i]*((l+1.)*(l+2.)/x/x-(l+1.)/(2*l+1.))
            -ptr2[i]/x*(l+1.)*(l+2.)/(2*l+1.);
            
			if(l==llist[j]){
            	sbesselj[j][i][0]=ptr3[i];
				sbesselj[j][i][1]=d3[i];
				sbesselj[j][i][2]=dd3[i];   
				j++	;	
            }
            
            l++;
            
        	ptr1[i] = (2*l-1)*ptr3[i]/x-ptr2[i];
        	d1[i]=ptr3[i]-(l+1)*ptr1[i]/x;
        	dd1[i]=d3[i]*l/(2*l+1.)+ptr1[i]*((l+1.)*(l+2.)/x/x-(l+1.)/(2*l+1.))
            -ptr3[i]/x*(l+1.)*(l+2.)/(2*l+1.);
            
			if(l==llist[j]){
            	sbesselj[j][i][0]=ptr1[i];
				sbesselj[j][i][1]=d1[i];
				sbesselj[j][i][2]=dd1[i];   
				j++;	
            }
            
            
            l++;
            
        	ptr2[i] = (2*l-1)*ptr1[i]/x-ptr3[i];
        	d2[i]=ptr1[i]-(l+1)*ptr2[i]/x;
        	dd2[i]=d1[i]*l/(2*l+1.)+ptr2[i]*((l+1.)*(l+2.)/x/x-(l+1.)/(2*l+1.))
            -ptr1[i]/x*(l+1.)*(l+2.)/(2*l+1.);
            
            
			if(l==llist[j]){
            	sbesselj[j][i][0]=ptr2[i];
				sbesselj[j][i][1]=d2[i];
				sbesselj[j][i][2]=dd2[i];   
				j++;		
            }
            
            l++;

            if((fabs(ptr2[i])<d)&&(x<l)) {l++; j++; break;} //This would mean it's 0, and hence all bigger l will be 0 too.
        }
        while(j<nls) {
            sbesselj[j][i][0]=0.;
            sbesselj[j][i][1]=0.;
            sbesselj[j][i][2]=0.;
            j++;
        }
//        printf("i=%d \n",i);
   }
//     
    printf("jl tabulation done!\n\n");
    
    free(ptr1);
    free(ptr2);
    free(ptr3);
    free(d1);
    free(d2);
    free(d3);
    free(dd1);
    free(dd2);
    free(dd3);    
    
    
    
    
    
    return 1;
}



/////////////////////////////////////////////////////////////////////////
//
//    The following code evaluate spherical Bessel function j_l (x)
//    from tabulation:
//
//      0 <= l <= lmax       0 <= x <= xmax
//
/////////////////////////////////////////////////////////////////////////

double get_sbesj(double sbesselj[][jlsample+1][3], int l, double x)
{
    if((l>lmax)||(x>xmax)) {
        printf("Spherical Besssel function out of range: l=%d, x=%f\n",l,x);
        return 0;
    }
    
    
    long double x1,x2;
    double v1,v2,d1,d2,res;
    int i=int(x/xinc);
    x1=i*xinc;
    x2=(i+1)*xinc;
    v1=sbesselj[l][i][0];
    v2=sbesselj[l][i+1][0];
    d1=sbesselj[l][i][1];
    d2=sbesselj[l][i+1][1];

    
    res=(-(v1*(x-x2)*(x-x2)*(2*x - 3*x1 + x2)) +
         (x - x1)*(v2*(x - x1)*(2*x + x1 - 3*x2) + (x - x2)*(x1 - x2)*(d1*x + d2*x - d2*x1 - d1*x2)))/
    (x1-x2)/(x1-x2)/(x1-x2);
    
    return res;
    
}




double interpol(double data[], double xtab[],int length, double x){
  //Linear interpolation algorithm
  
    double d1,res,v1,v2,x1,xmin,xmax;
    xmin=xtab[0];
    xmax=xtab[length-1];
    
    
    if((x<xmin)) {
        printf("interpol(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmin;
    }
    
    if((x>xmax)) {
        printf("interpol(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmax;
    }    
    
    

	int i=length/2,count; 
  
    for(count=2; x>=xtab[i+1] || x<xtab[i];++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol \n"); 
    }
    
    
    double step = xtab[i+1]-xtab[i]; 
    

//     printf("%e,%e, %f, %d \n",v1,v2,step,length);

    x1=xtab[i];
    v1=data[i];
    v2=data[i+1];
    d1=(v2-v1)/step;
    
    res= v1 + d1*(x-x1);
    
//     printf("x=%e,v1=%e, res=%e \n",x1,v1,res);
    
    return res;
    
}


double interpol2d(double data[][nls], double xtab[],int length, long x, int s){
// Special case where we have a diagonal part of the matrix and we want to calculate +-s away.
  
    double d1,res,v1,v2,x1,xmin,xmax;
    xmin=xtab[0];
    xmax=xtab[length-1];
    
    if((x<xmin)||(x>xmax)) {
        printf("interpol2d(x) is out of range. Range is x=(%f, %f) \n",xmin,xmax);
        return 0;
    }
    
  if ((s!=0) && (s!=1) && (s!=2))
  	printf("ERROR on interpol2d, what is s? \n");    

	int i=length/2,count; 
  
    for(count=2; x>=xtab[i+1] || x<xtab[i];++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol2d \n"); 
    }
    	
    	
    double step = xtab[i+1]-xtab[i]; 
      
    x1=xtab[i];
	v1=data[i][i+s-1];
    v2=data[i+1][i+s]; //The -1 is to account for the fact that s starts at 0.
    d1=(v2-v1)/step;  
  

    
    
	res= v1 + d1*(x-x1);    
    
    return res;      
         
    
}


double interpol3d(double data[][ncentral][ncentral], double xtab[],int length, int x,int y,int z){
// Bispectrum case, if l1+l2+l3 is odd cancels out, but we want to average it over.
  
    double res,v1,x1,y1,z1,xmin,xmax;
    double v2x,v2y,v2z,d1z,d1y,d1x;
    xmin=xtab[0];
    xmax=xtab[length-1];
    
    
    
    if((x<xmin)||(x>xmax)||(y<xmin)||(y>xmax)||(z<xmin)||(z>xmax)) {
        printf("interpol3d(x) is out of range. Range is x=(%f %f) \n",xmin,xmax);
        return 0;
    }
 

	int i=length/2,j=length/2,k=length/2,count; 
  
    for(count=2; x>=xtab[i+1] || x<xtab[i];++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol3d \n"); 
    }
    
    for(count=2; y>=xtab[j+1] || y<xtab[j];++count){
    	if (y>=xtab[j+1])
    		j+=1;
//    		i=i+length/(count);
    	else if (y<xtab[j])
			j-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol3d \n"); 
    } 
       
    for(count=2; z>=xtab[k+1] || z<xtab[k];++count){
    	if (z>=xtab[k+1])
    		k+=1;
//    		i=i+length/(count);
    	else if (z<xtab[k])
			k-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol3d \n"); 
    }    
    	
    	
 
      
    x1=xtab[i],y1=xtab[j],z1=xtab[k];

    double xstep = xtab[i+1]-xtab[i]; 
    double ystep = xtab[j+1]-xtab[j]; 
    double zstep = xtab[k+1]-xtab[k];
	v1=data[i][j][k];
    v2x=data[i+1][j][k],v2y=data[i][j+1][k],v2z=data[i][j][k+1];
    d1x=(v2x-v1)/xstep, d1y=(v2y-v1)/ystep, d1z=(v2z-v1)/zstep;  

  

    
    
	res= v1 + d1x*(x-x1)+ d1y*(y-y1)+ d1z*(z-z1);    
    
    return res;      
         
    
}






double nintegrate(double data[],double xtab[], int length){

	double sum = 0;
	long double step; //To avoid possible overflows when two values are too close.

	long i;
	
    for (i = 0; i < length-2; ++i) {
    	step = xtab[i+1]-xtab[i];
        sum +=  0.5*(data[i] + data[i+1]) * step;
//        printf("i= %d, sum=%lf \n",i,sum);
    }

return sum;

}



double Hfhigh(double q, long l1, long l2, double alpha, double r){
//Integral of 2 spherical bessel functions (j_l1 and j_l2) with a power k^(q+ns-1-alpha)

 	double gammas = lgamma(2.+alpha-q-ns) + lgamma((l1+l2+q+ns-alpha)/2.) - lgamma((l1+l2+4-q-ns+alpha)/2.)- lgamma((l1-l2 + 3-q-ns+alpha)/2.)- lgamma((-l1+l2+3-q-ns+alpha)/2.);

//	double gammas = lgamma(2+alpha-q-ns) - lgamma((l1-l2 + 3-q-ns+alpha)/2)- lgamma((-l1+l2+3-q-ns+alpha)/2) + (q+ns-alpha-2)*(log(l1+l2+4+alpha-q-ns)-log(2));

	double rest = log(PI) +(alpha-q-ns) * log(r) - (3.+alpha-q-ns)*log(2.);

	double result=exp(rest+gammas);

//	printf("gammas=%le, rest=%le, result=%le \n", gammas, rest, result);

	return result;
	



}

double Thetahigh(double q, long l1, long l2, double r,double S,double T21bar,double Tsbar,double DeltaZeta,double Tcmb,double z){
	
	
	double alpha,A;
//	double t0= S*Tcmb/Tsbar + 2*(3.*l1*l1+3*l1-2)/(4.*l1*l1+4*l1-3)*(1-Tcmb/Tsbar);
//	double tp=((l1+2.)*(l1+1))/((2.*l1+1)*(2*l1+3))*(1-Tcmb/Tsbar);
//	double tm=((l1*(l1-1.))/(4.*l1*l1-1))*(1-Tcmb/Tsbar);

	double	t0= S*Tcmb/(Tsbar-Tcmb) + 2*(3.*l1*l1+3*l1-2)/(4.*l1*l1+4*l1-3);
	double	tp=velocityswitch*((l1+2.)*(l1+1))/((2.*l1+1)*(2*l1+3)); //Change in the rest of codes.
	double	tm=velocityswitch*((l1*(l1-1.))/(4.*l1*l1-1));


//		alpha=0.6,A=250000.*pow(100/z,2); //within 10% OLD, redshift depedent.
	
	if (l1>6900)
		alpha=0.6,A=4.5*pow(10,4); //within 10% at redshift 0.
	else if(l1<=6900)
		printf("Error, l too low to use H!! \n");	
	else
		printf("Error, wtf is l? \n");
		
	double H1= (Hfhigh( q,  l1,  l2,  alpha,  r)-0.18*Hfhigh( q-1,  l1,  l2,  alpha,  r));
	double H2= (Hfhigh( q,  l1,  l2+2,  alpha,  r)-0.18*Hfhigh( q-1,  l1,  l2+2,  alpha,  r));
	double H3= (Hfhigh( q,  l1,  l2-2,  alpha,  r)-0.18*Hfhigh( q-1,  l1,  l2-2,  alpha,  r));
	
//	double res=(DeltaZeta*T21bar)*(H1*(S*Tcmb/Tsbar + 2*(3.*l1*l1+3*l1-2)/(4.*l1*l1+4*l1-3)*(1-Tcmb/Tsbar))-(((l1+2.)*(l1+1))/((2.*l1+1)*(2*l1+3))*(1-Tcmb/Tsbar))*H2-H3*(((l1*(l1-1.))/(4.*l1*l1-1))*(1-Tcmb/Tsbar)));

			
	double res=(2*PI*PI*DeltaZeta*T21bar)*A*A*(H1*t0-tp*H2-H3*tm);
	return res;
	
	
//Form of the H function at high l, with different values of alpha at different ls, and different amplitudes//

		
}

double Clhigh(double q, long l1, double r,double S,double T21bar,double Tsbar,double DeltaZeta,double Tcmb,double z){
	
	
	double alpha,A; //Power spectrum parametrization.
	double t0,tp,tm; //Aux variables

//		alpha=0.6,A=250000.*pow(100/z,2); //within 10% OLD, z dependent.
	
	if (l1>6900)
		alpha=0.6,A=4.5*pow(10,4);////within 10%
	else if(l1<=6900)
		printf("Error, l too low to use H!! \n");	
	else
		printf("Error, wtf is l? \n");
		
	double H1= (Hfhigh( q,  l1,  l1,  alpha,  r)-0.18*Hfhigh( q-1,  l1,  l1,  alpha,  r));
	double H2=(Hfhigh( q,  l1+2,  l1+2,  alpha,  r)-0.18*Hfhigh( q-1,  l1+2,  l1+2,  alpha,  r));
	double H3= (Hfhigh( q,  l1-2,  l1-2,  alpha,  r)-0.18*Hfhigh( q-1,  l1-2,  l1-2,  alpha,  r));
	double H4= (Hfhigh( q,  l1,  l1+2,  alpha,  r)-0.18*Hfhigh( q-1,  l1,  l1+2,  alpha,  r));
	double H5=(Hfhigh( q,  l1,  l1-2,  alpha,  r)-0.18*Hfhigh( q-1,  l1,  l1-2,  alpha,  r));
	double H6=(Hfhigh( q,  l1+2,  l1-2,  alpha,  r)-0.18*Hfhigh( q-1,  l1+2,  l1-2,  alpha,  r));	
	
		t0= S*Tcmb/(Tsbar-Tcmb) + 2*(3.*l1*l1+3*l1-2)/(4.*l1*l1+4*l1-3);
		tp=velocityswitch*((l1+2.)*(l1+1))/((2.*l1+1)*(2*l1+3)); //Change in the rest of codes.
		tm=velocityswitch*((l1*(l1-1.))/(4.*l1*l1-1));
					
	double	res=A*A*2./PI*(2*PI*PI*pow(kpivot,ns-1.)*DeltaZeta*T21bar*T21bar)*(H1*t0*t0+tp*tp*H2+H3*tm*tm-2*H4*t0*tp-2*H5*t0*tm-2*H6*tm*tp);

	return res;

		
}



double Threej(double l1,double l2,double l3){
	
	if(lround(l1+l2+l3)%2 !=0){
		return 0; //To make sure we only get the even permutations, that are not 0. 
		}
		//NOTE: We are not checking for triangle identities.
		
		
	long double log1;
	
	int sign=lround((l1+l2+l3)/2)%2;	
		
//	if(l1<100000){
		log1=lgamma((1+l1-l2+l3)/2.)+lgamma((1-l1+l2+l3)/2.)+lgamma((2+l1+l2+l3)/2.)-lgamma(l1+l2-l3+1.)-lgamma((l1-l2+l3)/2.+1)-lgamma((-l1+l2+l3)/2.+1.)-lgamma((3+l1+l2+l3)/2.) - 2.*lgamma((1-l1-l2+l3)/2.) ;
		//}
// 	else{	
// 		log1=-1/2.*log((2+l1-l2+l3)/2.)-1/2.*log((2-l1+l2+l3)/2.)-1/2.*log((3+l1+l2+l3)/2.)-(l1+l2-l3)*log(l1+l2-l3)+(l1+l2-l3)*log((l1+l2-l3)/2.);
// 	}
	
	log1=log1+2.*(-1.+l1+l2-l3)/2.*log(2);
	
	double res=exp(log1/2.)*pow(-1,sign)*pow(PI,1/4.); //*pow(2.,(-1.+l1+l2-l3)/2.) included right above
	
	return res;
	
		
}

int triangle(double l1,double l2,double l3){
	//Checking for triangle identities. Returns 0 if not valid and 1 if valid.
	int res=1;
	
	if((l1+l2-l3<0) || (l1-l2+l3<0) || (-l1+l2+l3<0)){
		res=0;	
		}
		
	return res;
}


int delta(double r2, double r3){
	//Supposedly we have k1, k2=r2*k1 and k3=k2*r3.
	//First we check for triangle identities. Returns 0 if not valid and runs the function if valid.
	int res;
	
	
	if((r2<0.5)||(r2>1.)||(r3>1.)||(r3<1./r2-1.)){
		res=0;	
		printf("Error on delta(), r2 and r3 do not satisfy triangle identities \n");
	}
			
	
	if((r2!=1)&&(r3!=1)){
		res=6; //Three sides different.
	}
	else if((r3==1)||(r2==1)){
		res=3;//Two sides are equal.
	}
	else if((r3==1)&&(r2==1)){
		res=1;//Three sides the same.
	}	


		
	return res;
}


double Sixj1(double l1,double l2,double l3,int s1, int s2){
	
// This is the 6 j symbol for l4=1, so s1=l5-l1 and s2=l6-l2 have to be +-1.	
	
	double res=0;
	
	double l=l1+l2+l3;
	
	if(s2==-1){
		if(s1<0)
			res=pow(-1,l)*sqrt(l*(l+1.)*(l-2*l3-1)*(l-2*l3)/((2*l1-1)*2*l1*(2*l1+1)*(2*l2-1)*2*l2*(2*l2+1)));
		else
			res=pow(-1,l)*sqrt((l-2*l1-1.)*(l-2*l1)*(l-2*l2+1.)*(l-2*l2+2)/((2*l1+1)*(2*l1+2)*(2*l1+3)*(2*l2-1)*2*l2*(2*l2+1)));			
	}
	else if(s2==1){
		l2+=1;
		l1+=s1;
		l=l1+l2+l3;
		if(s1>0) //It's inverted now!
			res=pow(-1,l)*sqrt(l*(l+1.)*(l-2*l3-1)*(l-2*l3)/((2*l1-1)*2*l1*(2*l1+1)*(2*l2-1)*2*l2*(2*l2+1)));
		else
			res=pow(-1,l)*sqrt((l-2*l1-1.)*(l-2*l1)*(l-2*l2+1.)*(l-2*l2+2)/((2*l1+1)*(2*l1+2)*(2*l1+3)*(2*l2-1)*2*l2*(2*l2+1)));			
	}
		
		

	return res;
	
		
}


double Sixj2(double l1,double l2,double l3,int s1, int s2){
	
// This is the 6 j symbol for l4=2, so s1=l5-l1 and s2=l6-l2 have to be +-2 or 0. From tables.
	
	double res=0;
	
	long l=l1+l2+l3;
	
	int sign=l%2;
	
	double x; //For the 0-0 case.
	
	if(s2==-2){
		if(s1<0)
			res=pow(-1,sign)*sqrt(l*(l-2.)*(l-1)*(l+1)*(l-2*l3-3)*(l-2.*l3-2)*(l-2*l3-1)*(l-2*l3)/((2.*l1-3)*(2*l1-2)*(2*l1-1)*(2*l1)*(2*l1+1)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1))); //s1=-2, s2=-2 
		else if (s1==0) 
			res=pow(-1,l)*sqrt(6.*l*(l+1.)*(l-2*l3-1)*(l-2*l3)*(l-2*l1-1)*(l-2*l1)*(l-2*l2+1)*(l-2*l2+2)/((2.*l1+3)*(2*l1+2)*(2*l1+1)*(2*l1)*(2*l1-1)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1)));// 0 -2
		else
			res=pow(-1,l)*sqrt((l-2.*l1-3)*(l-2*l1-2)*(l-2*l1-1)*(l-2*l1)*(l-2*l2+1)*(l-2*l2+2)*(l-2*l2+3)*(l-2*l2+4)/((2.*l1+3)*(2*l1+2)*(2*l1+1)*(2*l1+4)*(2*l1+5)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1))); //2 -2
	}
	else if (s2==0){ 
		if(s1<0){
			x=l1; //temp storage//
			l1=l2;
			l2=x;
			res=pow(-1,l)*sqrt(6.*l*(l+1.)*(l-2*l3-1)*(l-2*l3)*(l-2*l1-1)*(l-2*l1)*(l-2*l2+1)*(l-2*l2+2)/((2.*l1+3)*(2*l1+2)*(2*l1+1)*(2*l1)*(2*l1-1)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1)));// -2 0 (with l1 and l2 swapped)
			}
		else if (s1==0){ 
			x=l1*(l1+1.)+l2*(l2+1.)-l3*(l3+1.);
			res=pow(-1,l)*2.*(3.*x*(x-1.)-4*l1*(l1+1)*l2*(l2+1))/sqrt((2.*l1-1)*(2.*l1)*(2*l1+1)*(2*l1+2)*(2*l1+3)*(2*l2-1)*(2.*l2)*(2*l2+1)*(2*l2+2)*(2*l2+3));//0 0
			}
		else{
			l1+=s1; //to account for the inversion//
			l=l1+l2+l3;
			x=l1; //temp storage//
			l1=l2;
			l2=x;
			res=pow(-1,l)*sqrt(6.*l*(l+1.)*(l-2*l3-1)*(l-2*l3)*(l-2*l1-1)*(l-2*l1)*(l-2*l2+1)*(l-2*l2+2)/((2.*l1+3)*(2*l1+2)*(2*l1+1)*(2*l1)*(2*l1-1)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1)));// 2 0 (with l1 and l2 swapped)
			}
	}
	else if(s2==2){ 
		l2+=s2;
		l1+=s1;		
		l=l1+l2+l3;
		sign=l%2;
		if(s1>0)
			res=pow(-1,sign)*sqrt(l*(l-2.)*(l-1)*(l+1)*(l-2*l3-3)*(l-2.*l3-2)*(l-2*l3-1)*(l-2*l3)/((2.*l1-3)*(2*l1-2)*(2*l1-1)*(2*l1)*(2*l1+1)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1))); //s1=-2, s2=2 
		else if (s1==0) 
			res=pow(-1,l)*sqrt(6.*l*(l+1.)*(l-2*l3-1)*(l-2*l3)*(l-2*l1-1)*(l-2*l1)*(l-2*l2+1)*(l-2*l2+2)/((2.*l1+3)*(2*l1+2)*(2*l1+1)*(2*l1)*(2*l1-1)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1)));// 0 2
		else
			res=pow(-1,l)*sqrt((l-2.*l1-3)*(l-2*l1-2)*(l-2*l1-1)*(l-2*l1)*(l-2*l2+1)*(l-2*l2+2)*(l-2*l2+3)*(l-2*l2+4)/((2.*l1+3)*(2*l1+2)*(2*l1+1)*(2*l1+4)*(2*l1+5)*(2*l2-3)*(2*l2-2)*(2*l2-1)*(2*l2)*(2*l2+1))); //2 2
	}
	
		
		

	return res;
	
		
}


double interpol_cubic(double x0, double dx, double *ytab, unsigned int Nx, double x) {
  //cubic interpolation routine, assumes x equally spaced.
  long ix;
  double frac;

  if (Nx < 4) {
    fprintf(stderr, "Error: interpol_cubic: Table needs to be of dimension 4 at least\n");
    exit(1);
  }

  // Check if in range
  if (  (dx > 0 && (x<x0 || x>x0+dx*(Nx-1))) 
      ||(dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) ) {
	fprintf(stderr, "Error: interpol_cubic: x-value out of range. Range is (%le,%le) and x=%le. \n",x0,x0+dx*(Nx-1),x);
 //   exit(1);
  }
 
  // Identify location to interpolate
  ix = (long)floor((x-x0)/dx);
  if (ix<1) ix=1; 
  if (ix>Nx-3) ix=Nx-3;
  frac = (x-x0)/dx-ix;
  ytab += ix-1;

  // Return value
  return(
    -ytab[0]*frac*(1.-frac)*(2.-frac)/6.
    +ytab[1]*(1.+frac)*(1.-frac)*(2.-frac)/2.
    +ytab[2]*(1.+frac)*frac*(2.-frac)/2.
    -ytab[3]*(1.+frac)*frac*(1.-frac)/6.
  );
}


double interpol_2d_cubic(double x0, double dx, double data[npoints][npoints], int Nx, double x, double y) {
  
  long ix,iy;
  
  if(npoints!=Nx) printf("Error on interpol_2d_cubic, mismatching dimensions \n"); 

  if (Nx < 4) {
    fprintf(stderr, "Error: interpol_2d_cubic: Table needs to be of dimension 4 at least\n");
    exit(1);
  }

// Check if in range
  if (  (dx > 0 && (x<x0 || x>x0+dx*(Nx-1))) 
      ||(dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_2d_cubic: x-value out of range in interpolation routine.\n");
    exit(1);
  }
  if (  (dx > 0 && (y<x0 || y>x0+dx*(Nx-1))) 
      ||(dx < 0 && (y>x0 || y<x0+dx*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_2d_cubic: y-value out of range in interpolation routine.\n");
    exit(1);
  }  
 
// Identify location to interpolate
  ix = (long)floor((x-x0)/dx);
	if (ix<1) ix=0; 
	if (ix>Nx-3) ix=Nx-3;
  iy = (long)floor((y-x0)/dx);
	if (iy<1) iy=0; 
	if (iy>Nx-3) iy=Nx-3;  
 // 
  
  int i;
  double vx[4];
  
  for(i=0;i<4;i++){
	vx[i]=interpol_cubic(x0,dx,data[ix+i], Nx, y);
  }
  
  double res=interpol_cubic(x0+(ix)*dx,dx,vx,4,x);
  
   return res;
  
}


double interpol_2d_aux(double x0,double y0, double dx,double dy, double **data, double x, double y) {
  //2d interpolation to use with unknown-dimension arrays.
  int Nx=4; //To use in the 3d case.
  
  long ix,iy;
  

// Check if in range
  if (  (dx > 0 && (x<x0 || x>x0+dx*(Nx-1))) 
      ||(dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_2d_aux: x-value out of range. Range is x=(%f, %f) and x=%le \n",x0,x0+dx*(Nx-1),x);
//    exit(1);
  }
  if (  (dy > 0 && (y<y0 || y>y0+dy*(Nx-1))) 
      ||(dy < 0 && (y>y0 || y<y0+dy*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_2d_aux: y-value out of range. Range is x=(%f, %f) and x=%le \n",y0,y0+dy*(Nx-1),y);
//    exit(1);
  }  
 
// Identify location to interpolate
  ix = (long)floor((x-x0)/dx);
//  printf("ix=%ld \t\t\t",ix);
	if (ix<1) ix=0; 
	if (ix>Nx-3) ix=Nx-4; 
  iy = (long)floor((y-y0)/dy);
//  printf("iy=%ld \n",iy);  
	if (iy<1) iy=0; 
	if (iy>Nx-3) iy=Nx-4; 
 
 
// printf("ix=%ld and iy= %ld \n",ix,iy);
  
  int i;
  double vx[4];
  
  for(i=0;i<4;i++){
	vx[i]=interpol_cubic(y0,dy,data[i],4,y);
  }
  
	double res=interpol_cubic(x0,dx,vx,4,x);
  
	return res;
  
}


double interpol_3d_cubic(double x0, double dx, double data[npoints][npoints][npoints], int Nx, double x, double y, double z) {
  //given a npoints^3 array n k-space.
  long ix,iy,iz;
  
  if(npoints!=Nx) printf("Error on interpol_3d_cubic, mismatching dimensions \n"); 

  if (Nx < 4) {
    fprintf(stderr, "Error: interpol_3d_cubic: Table needs to be of dimension 4 at least\n");
    exit(1);
  }

  // Check if in range
  if (  (dx > 0 && (x<x0 || x>x0+dx*(Nx-1))) 
      ||(dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_3d_cubic: x-value out of range in interpolation routine.");
    printf("limit is (%le,%le)  and you are inputting %le. \n",x0,x0+dx*(Nx-1),x);
    exit(1);
  }
  if (  (dx > 0 && (y<x0 || y>x0+dx*(Nx-1))) 
      ||(dx < 0 && (y>x0 || y<x0+dx*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_3d_cubic: y-value out of range in interpolation routine.\n");
    printf("limit is (%le,%le)  and you are inputting %le. \n",x0,x0+dx*(Nx-1),y);
    exit(1);
  }
  if (  (dx > 0 && (z<x0 || z>x0+dx*(Nx-1))) 
      ||(dx < 0 && (z>x0 || z<x0+dx*(Nx-1))) ) {
    fprintf(stderr, "Error: interpol_3d_cubic: y-value out of range in interpolation routine.\n");
    printf("limit is (%le,%le)  and you are inputting %le. \n",x0,x0+dx*(Nx-1),z);
    exit(1);
  }
  
      
 
  // Identify location to interpolate
  ix = (long)floor((x-x0)/dx);
	if (ix<1) ix=0; 
	if (ix>Nx-3) ix=Nx-3;
  iy = (long)floor((y-x0)/dx);
	if (iy<1) iy=0; 
	if (iy>Nx-3) iy=Nx-3;    
  iz = (long)floor((z-x0)/dx);
	if (iz<1) iz=0; 
	if (iz>Nx-3) iz=Nx-3;  
 
  
  int i,j,k;
  double vx[4];
  double dat[4][4];
  

	if(ix>Nx-4) ix=Nx-4;
	for(i=0;i<4;i++){
	  for(j=0;j<4;j++){  	
		  for(k=0;k<4;k++){  	
  	 		 dat[j][k]=data[ix+i][iy+j][iz+k];
  	  	}
  	  }
//	vx[i]=interpol_2d_aux(x0+iy*dx,x0+iz*dx,dx,dx,dat,y,z);
  }
  
  double res=interpol_cubic(x0+ix*dx,dx,vx,4,x);
  
  return res;
  
}


double interpol_3d_triang(double xtab[],double data[][npoints][npoints],int length, double x,double y,double z){
  //checks for triangle identities.
    double res,v1,x1,y1,z1,xmin,xmax;
    double v2x,v2y,v2z,d1z,d1y,d1x;
    xmin=xtab[0];
    xmax=xtab[length-1];
    
    
    
    if((x<xmin)||(x>xmax)||(y<xmin)||(y>xmax)||(z<xmin)||(z>xmax)) {
        printf("interpol_3d_triang(x) is out of range. Range is x=(%f, %f) \n",xmin,xmax);
        return 0;
    }
 

	long i=length/2,j=length/2,k=length/2,count; 
  
    for(count=2; x>=xtab[i+1] || x<xtab[i];++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol_3d_triang \n"); 
    }
    
    for(count=2; y>=xtab[j+1] || y<xtab[j];++count){
    	if (y>=xtab[j+1])
    		j+=1;
//    		i=i+length/(count);
    	else if (y<xtab[j])
			j-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol_3d_triang \n"); 
    } 
       
    for(count=2; z>=xtab[k+1] || z<xtab[k];++count){
    	if (z>=xtab[k+1])
    		k+=1;
//    		i=i+length/(count);
    	else if (z<xtab[k])
			k-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol_3d_triang \n"); 
    }   
    
    
    x1=xtab[i],y1=xtab[j],z1=xtab[k];     
	v1=data[i][j][k];    
    int di=0,dj=0,dk=0;	
    // Now we check for the next non-zero point.
    
    if(triangle(xtab[i+1],y1,z1)){
    	di=1;
    }
    if(triangle(x1,xtab[j+1],z1)){
    	dj=1;
    }
    if(triangle(x1,y1,xtab[k+1])){
    	dk=1;
    }    	
 
    double xstep;

	if(di){
    xstep = (xtab[i+di]-xtab[i]); 
    v2x=data[i+di][j][k];
    d1x=(v2x-v1)/xstep;
    }
    else {
    d1x=0;
    }
    
    if(dj){
    xstep = (xtab[j+dj]-xtab[j]); 
    v2y=data[i][j+dj][k];
    d1y=(v2y-v1)/xstep;
    }
    else {
    d1y=0;
    }
    
    if(dk){
    xstep = (xtab[k+dk]-xtab[k]); 
    v2z=data[i][j][k+dk];
    d1z=(v2z-v1)/xstep;
    }
    else {
    d1z=0;
    }
    
    
    
	res= v1 + d1x*(x-x1)+ d1y*(y-y1)+ d1z*(z-z1);    

    
    return res;      
         
    
}


double interpol_3d_ratio(double xtab[],long double r20, long double r30,double ***data,long lx, long lr, double x,long double r2,long double r3){
  //k2=r2*k and k3=r3*k2.
  //r20 and r30 are the initial r2 and r3 (0.5 and 1/r2-1. respectively).
  //We will take k ((xtab[])) in log space.
    double res,v1,x1,xmin,xmax;
    double d1z,d1y,d1x;
    xmin=xtab[0];
    xmax=xtab[lx-1];
    
    
    
    if((x<xmin)||(x>xmax)) {
        printf("interpol_3d_ratio(x) is out of range. Range is x=(%f, %f) and x=%le \n",xmin,xmax,x);
        return 0;
    }
    if((r2<r20)||(r2>1.)) {
        printf("interpol_3d_ratio(x) is out of range. Range is r2=(%f, %f) and r2=%.2Lf \n",.5,1.,r2);
        return 0;
    }
    if((r3<r30)||(r3>1.)) {
        printf("interpol_3d_ratio(x) is out of range. Range is r3=(%Lf, %f) and r3=%.2Lf \n",1./r2-1.,1.,r3);
        return 0;
    }        
    
 

	long ix=lx/2;
	long ir2, ir3;
  
    for(; x>=xtab[ix+1] || x<xtab[ix];){
    	if (x>=xtab[ix+1])
    		ix+=1;
    	else if (x<xtab[ix])
			ix-=1;
    	else
    		printf("ERROR on interpol_3d_ratio \n"); 
    }
    
    //We assume that the rlist increase linearly.
    long double stepr2=(1.-r20)/(lr-1.);
    long double stepr3=(1.-r30)/(lr-1.);
    double stepx=(xmax-xmin)/(lx-1.);
    
    ir2= (long)floor((r2-r20)/stepr2);
    if (ir2<1) ir2=0; 
	if (ir2>lr-4) ir2=lr-4;
    ir3= (long)floor((r3-r30)/stepr3);   
//    if(r3<(r30+ir3*stepr3)) r3=(r30+ir3*stepr3); 
    if (ir3<1) ir3=0; 
	if (ir3>lr-4) ir3=lr-4;
	
    if (ix<1) ix=0; 
	if (ix>lx-4) ix=lx-4;  	


	int j,k,m;
	double **datanew;
	
	// datanew=(double**)calloc(4, sizeof(double *));			
// 
// 	for(j=0;j<4;++j){
// 		datanew[j]=(double*)calloc(4, sizeof(double));
// 	}
	
	datanew=create_2D_array(4,4);
	
	double v2k[4]; //To interpolate in k space.
	
	for(m=0;m<4;m++){
		for(j=0;j<4;j++){
			for(k=0;k<4;k++){
//				if((ix+m>lx)||(j+ir2>lr)||(k+ir3>lr)) printf("error 1 \n");
//				if((ix+m<0)||(j+ir2<0)||(k+ir3<0)) printf("error 2 \n");
				datanew[j][k]=data[ix+m][j+ir2][k+ir3];
			}
		}	
		v2k[m]=interpol_2d_aux(r20+stepr2*ir2,r30+ir3*stepr3,stepr2,stepr3, datanew, r2, r3);		
	}
	
	free_2D_array(datanew,4);
	  

	
	res = interpol_cubic(xtab[ix], stepx, v2k, 4, x);//v1+ dvk*dk;// + dvr2*dr2 + dvr3*dr3;
	
	return res;      
         
}



void solve_syst(double  **A, double  *X, double  *B, int N) {
  if (N == 1) {
     if (A[0][0] == 0) {
       fprintf(stderr, "Error: solve_syst: singular system (N = 1).\n");
       exit(1);
     }
     X[0] = B[0]/A[0][0];
     return;
  }
  
  double  Apivot = A[N-1][N-1];
  int ipivot = N-1;

  /* Find a row for which A[i, N-1] is not zero */
  while (Apivot == 0. && ipivot > 0) {
    ipivot--;
    Apivot = A[ipivot][N-1];
  }
  if (Apivot == 0. && ipivot == 0 ) {
    fprintf(stderr, "Error: solve_syst: singular system (N = %i).\n", N);
      //exit(1);
  }

  /* Assign a new N-1 by N-1 matrix with all columns except N-1 and 
     all rows except ipivot */ 
  double  **newA = create_2D_array(N-1, N-1);
  double  *newB  =  create_1D_array(N-1);
  int i, j;

  for (j = 0; j < N-1; j++) {
     for (i = 0; i < ipivot; i++)     newA[i][j] = A[i][j] - A[i][N-1]*A[ipivot][j]/Apivot;
     for (i = ipivot+1; i < N; i++) newA[i-1][j] = A[i][j] - A[i][N-1]*A[ipivot][j]/Apivot;
     newB[j] = B[j] - A[j][N-1] * B[ipivot]/Apivot;
  }
  solve_syst(newA, X, newB, N-1);

  X[N-1] = B[ipivot]/Apivot;
  for (j = 0; j < N-1; j++) X[N-1] -= A[ipivot][j]*X[j]/Apivot;
  
  free_2D_array(newA, N-1);
  free(newB);
}

void reverse(double *list,int N){
	double aux;
	int j;
	for(j=0;j<=N/2;j++){
		aux=list[j];
		list[j]=list[N-j];
		list[N-j]=aux;
	}

}

double det(double **M, int N){
//Computes the determinant of a NxN matrix M
	if(N==1) return M[0][0];
		
	double res;
	double **Maux;
	Maux=create_2D_array(N-1,N-1);
	
	int i,j,k;
	
	int jaux; //0 if j<k, 1 if j>=k, to not copy the k column.
	
	for(k=0;k<N;k++){
		for(i=0;i<N-1;i++){	
			for(j=0;j<N-1;j++){
			(j<k)?(jaux=0):(jaux=1);
			Maux[i][j] = M[i+1][j+jaux];
			}
		}
	res+=pow(-1,k)*M[k][0]*det(Maux,N-1);
	}

	free_2D_array(Maux,N-1);
	
	return res;
}



