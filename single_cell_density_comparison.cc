// g++ single_cell_density_comparison.cc -o single_cell_density_comparison -I/usr/local/include/voro++/ -L/usr/local/lib -lvoro++

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

#include "voro++.hh"
using namespace voro;


// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

const double h=0.5, pi=M_PI;
const double Nc=100.0;

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=50;

int id_list[particles];
double x_list[particles], y_list[particles], z_list[particles];

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}


///////////////////////////

double F1_3d(double phi, double r0, double R_0, double B1) {

  double integral;
  double mu, a, logs, invtan, u;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h2 = pow(r0/h,2);
  r0h3 = pow(r0/h,3);
  r0h_2 = pow(r0/h,-2);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  I_2 = phi +   a2 * tan(phi);
  I_4 = phi + 2*a2 * tan(phi) + 1./3.*pow(a,4) * tan(phi)*(2. + 1./cosp2);

  //u = sqrt(1.-(1.+a2)*pow(mu,2));
  u = sin(phi)*sqrt(1-mu*mu);
  logs = log(1+u) - log(1-u);
  invtan = atan(u/a);
  I1  = invtan;

  I_1 = a/2.*logs + invtan;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-pow(u,2)) + logs);
  I_5 = I_3 + a*pow((1.+a2),2)/16. *( (10*u - 6*pow(u,3))/pow(1-pow(u,2),2) + 3.*logs);

  integral =  r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);

  /*
  printf("B1: %g\n", B1);
  printf("phi: %g\n", phi);
  printf("I_2: %g\n", I_2);
  printf("I_4: %g\n", I_4);
  printf("I_5: %g\n", I_5);
  printf("I0: %g\n\n", I0);
  */

  return integral;
}


double F2_3d(double phi, double r0, double R_0, double B2) {

  double integral;
  double mu, a, logs, invtan, u;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h2 = pow(r0/h,2);
  r0h3 = pow(r0/h,3);
  r0h_2 = pow(r0/h,-2);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  I_2 = phi +   a2 * tan(phi);
  I_4 = phi + 2*a2 * tan(phi) + 1./3.*pow(a,4) * tan(phi)*(2. + 1./cosp2);

  //u = sqrt(1.-(1.+a2)*pow(mu,2));  The two expressions for u are equivallent. The bottom one is chosen for numerical reasons.
  u = sin(phi)*sqrt(1-mu*mu);
  logs = log(1+u) - log(1-u);
  invtan = atan(u/a);
  I1  = invtan;

  I_1 = a/2.*logs + invtan;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-pow(u,2)) + logs);
  I_5 = I_3 + a*pow((1.+a2),2)/16. *( (10*u - 6*pow(u,3))/pow(1-pow(u,2),2) + 3.*logs);

  integral =  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
					  1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0);

  /*
  printf("B2: %g\n", B2);
  printf("phi: %g\n", phi);
  printf("I_2: %g\n", I_2);
  printf("I_3: %g\n", I_3);
  printf("I_4: %g\n", I_4);
  printf("I_5: %g\n", I_5);
  printf("I0: %g\n", I0);
  printf("I1: %g\n\n", I1);
  */

  return integral;
}

double F3_3d(double phi, double r0, double R_0, double B3) {

  double integral;
  double I0, I1;
  double a, a2, cosp2, r03, r0h3, r0h_3, u, invtan, mu;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h3 = pow(r0/h,3);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  u = sin(phi)*sqrt(1-mu*mu);
  invtan = atan(u/a);
  I1  = invtan;

  integral = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0);

  /*
  printf("B3: %g\n", B3);
  printf("phi: %g\n", phi);
  printf("I0: %g\n", I0);
  printf("I1: %g\n\n", I1);
  */

  return integral;
}


double full_integral(double phi, double r0, double R_0) {

  double B1, B2, B3, mu, a, logs, invtan, u;
  double full_int;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;
  double r, R, linedist, phi1, phi2;

  if(r0==0.0) return 0.0;
  if(R_0==0.0) return 0.0;
  if(phi==0.0) return 0.0;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h2 = pow(r0/h,2);
  r0h3 = pow(r0/h,3);
  r0h_2 = pow(r0/h,-2);
  r0h_3 = pow(r0/h,-3);


  if(r0 >= 2.0*h) {
    B3 = pow(h,3) /4.;
  }
  else if(r0 > h) {
    B3 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2);
    B2 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3);
  }
  else {
    B3 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2);
    B2 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2);
    B1 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3);
  }


  /*
  if(r0/mu < 1.0*h) {
    full_int = F1_3d(phi,r0,R_0,B1);
  }
  else if(r0/mu < 2.0*h) {
    full_int = F2_3d(phi,r0,R_0,B2);
  }
  else {
    full_int = F3_3d(phi,r0,R_0,B3);
  }

  */

  linedist = sqrt(r0*r0 + R_0*R_0);
  R = R_0/cos(phi);
  r = sqrt(r0*r0 + R*R);


    full_int = 0.0;

    if(linedist <= 1.0*h) {
      phi1 = acos(R_0/sqrt(h*h-r0*r0));
      phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));

      if(r <= 1.0*h) {
        full_int = full_int + F1_3d(phi, r0, R_0, B1) - F1_3d(0.0, r0, R_0, B1);
      }
      else if(r <= 2.0*h) {
        full_int = full_int + F1_3d(phi1, r0, R_0, B1) - F1_3d(0.0, r0, R_0, B1);
        full_int = full_int + F2_3d(phi, r0, R_0, B2) - F2_3d(phi1, r0, R_0, B2);
      }
      else { 
        full_int = full_int + F1_3d(phi1, r0, R_0, B1) - F1_3d(0.0, r0, R_0, B1);
        full_int = full_int + F2_3d(phi2, r0, R_0, B2) - F2_3d(phi1, r0, R_0, B2);
        full_int = full_int + F3_3d(phi, r0, R_0, B3) - F3_3d(phi2, r0, R_0, B3);
      }
    }


    if((linedist>1.0*h) && (linedist<=2.0*h)) {
      phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));

      if(r <= 2.0*h) {
        full_int = full_int + F2_3d(phi, r0, R_0, B2) - F2_3d(0.0, r0, R_0, B2);
      }
      else {
        full_int = full_int + F2_3d(phi2, r0, R_0, B2) - F2_3d(0.0, r0, R_0, B2);
        full_int = full_int + F3_3d(phi, r0, R_0, B3) - F3_3d(phi2, r0, R_0, B3);
      }
    }



    if(linedist > 2.0*h) {
      full_int = full_int + F3_3d(phi, r0, R_0, B3) - F3_3d(0.0, r0, R_0, B3);
    }

  return full_int;
}






double full_integral_new(double phi, double r0, double R_0) {

  double B1, B2, B3, mu, a, logs, u;
  double full_int;
  double a2, cosp, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3, tanp;
  double r, R, linedist, phi1, phi2;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double D1, D2, D3;

  if(r0==0.0) return 0.0;
  if(R_0==0.0) return 0.0;
  if(phi==0.0) return 0.0;

  r03 = r0*r0*r0;
  r0h2 = r0/h*r0/h;
  r0h3 = r0h2*r0/h;
  r0h_2 = h/r0*h/r0;
  r0h_3 = r0h_2*h/r0;


  if(r0 >= 2.0*h) {
    B3 = h*h*h /4.;
  }
  else if(r0 > h) {
    B3 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2);
    B2 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3);
  }
  else {
    B3 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2);
    B2 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2);
    B1 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3);
  }


  a = R_0/r0;
  a2 = a*a;

  linedist = sqrt(r0*r0 + R_0*R_0);
  R = R_0/cos(phi);
  r = sqrt(r0*r0 + R*R);


  full_int = 0.0;
  D2 = 0.0;
  D3 = 0.0;

  if(linedist <= 1.0*h) {
    ////// phi1 business /////
    phi1 = acos(R_0/sqrt(h*h-r0*r0));
    
    cosp = cos(phi1);
    cosp2 = cosp*cosp;
    mu = cosp/a / sqrt(1. + cosp2/a2);
    
    tanp = tan(phi1);
    
    I0  = phi1;
    I_2 = phi1 +   a2 * tanp;
    I_4 = phi1 + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);
    
    u = sin(phi1)*sqrt(1-mu*mu);
    logs = log(1+u) - log(1-u);
    I1 = atan(u/a);
    
    I_1 = a/2.*logs + I1;
    I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
    I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);
    
    D2 = -1./6.*I_2 + 0.25*(r0/h) *I_3 - 0.15*r0h2 *I_4 + 1./30.*r0h3 *I_5 - 1./60. *r0h_3 *I1 + (B1-B2)/r03 *I0;
    
    
    ////// phi2 business /////
    phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));
    
    cosp = cos(phi2);
    cosp2 = cosp*cosp;
    mu = cosp/a / sqrt(1. + cosp2/a2);
    
    tanp = tan(phi2);
    
    I0  = phi2;
    I_2 = phi2 +   a2 * tanp;
    I_4 = phi2 + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);
    
    u = sin(phi2)*sqrt(1-mu*mu);
    logs = log(1+u) - log(1-u);
    I1 = atan(u/a);
    
    I_1 = a/2.*logs + I1;
    I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
    I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);
    
    D3 = 1./3.*I_2 - 0.25*(r0/h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2;
  }
  else if(linedist <= 2.0*h) {
    ////// phi2 business /////
    phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));
    
    cosp = cos(phi2);
    cosp2 = cosp*cosp;
    mu = cosp/a / sqrt(1. + cosp2/a2);
    
    tanp = tan(phi2);
    
    I0  = phi2;
    I_2 = phi2 +   a2 * tanp;
    I_4 = phi2 + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);
    
    u = sin(phi2)*sqrt(1-mu*mu);
    logs = log(1+u) - log(1-u);
    I1 = atan(u/a);
    
    I_1 = a/2.*logs + I1;
    I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
    I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);
    
    D3 = 1./3.*I_2 - 0.25*(r0/h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2;
  }
  
  
  printf("phi1, phi2: %g %g\n", phi1, phi2);
  printf("D2, D3: %g %g\n", D2, D3);
  printf("linedist/h, r/h: %g %g\n", linedist/h, r/h);
  printf("%g %g %g %g %g %g %g\n", I0, I1, I_1, I_2, I_3, I_4, I_5);


  //////////////////////////////


  cosp = cos(phi);
  cosp2 = cosp*cosp;
  mu = cosp/a / sqrt(1. + cosp2/a2);

  tanp = tan(phi);

  I0  = phi;
  I_2 = phi +   a2 * tanp;
  I_4 = phi + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);

  u = sin(phi)*sqrt(1-mu*mu);
  logs = log(1+u) - log(1-u);
  I1 = atan(u/a);

  I_1 = a/2.*logs + I1;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
  I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);

  
  //if(r0/mu < 1.0*h) {
  if(r < 1.0*h) {
    //full_int = F1_3d(phi,r0,R_0,B1);
    full_int = r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);
  }
  //else if(r0/mu < 2.0*h) {
  else {
    if(r < 2.0*h) {
      //full_int = F2_3d(phi,r0,R_0,B2);
      full_int=  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
				     1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2);
    }
    else {
      //full_int = F3_3d(phi,r0,R_0,B3);
      full_int = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3);
    }
  }  

  //printf("mu: %g, r0/mu: %g\n",mu,r0/mu);
  //printf("I1: %g %g\n", I1, pi/2.- asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) ));

  return full_int;
}


///////////////////////////////////
/// Numerical solution ///////////
/////////////////////////////////

double kernel(double y, double x, double z) {
  double r, W;

  r = sqrt(x*x+y*y+z*z);
  W = 1.0/h/h/h/pi;

  if(r < h) {
    W = W * (1 - 1.5*(r/h)*(r/h) + 0.75*(r/h)*(r/h)*(r/h) );
  }
  else {
    if(r < 2*h) {
      W = W * 0.25*(2.0-r/h)*(2.0-r/h)*(2.0-r/h);
    }
    else {
      W = 0.0;
    }
  }

  return W;
}

double line(double x, double z, double phi) {
  int i,N;
  double dy, y1, y2, ym, f1, f2, fm, side;

  N = int(Nc*x*tan(phi)/h)+1;
  dy = x*tan(phi)/double(N);
  y1 = 0.0;
  f1 = kernel(y1,x,z);
  side = 0.0;

  for(i=1;i<=N;i++){
    y2 = double(i)*dy;
    ym = (y1+y2)/2.0;
    f2 = kernel(y2,x,z);
    fm = kernel(ym,x,z);
    side = side + (f1+4.0*fm+f2)/6.0*(y2-y1);
    y1 = y2;
    f1 = f2;
  }

  return side;
}

double area(double z, double r0, double R_0, double phi) {
  int i, N;
  double dx, x1, x2, xm, f1, f2, fm, R, S;

  R = z/r0 * R_0;
  N = int(Nc*R/h)+1;
  dx = R/double(N);
  x1 = 0.0;
  f1 = line(x1,z,phi);
  S = 0.0;

  for(i=1;i<=N;i++){
    x2 = double(i)*dx;
    xm = (x1+x2)/2.0;
    f2 = line(x2,z,phi);
    fm = line(xm,z,phi);
    S = S + (f1+4.0*fm+f2)/6.0*(x2-x1);
    x1 = x2;
    f1 = f2;
  }

  return S;
}

double volume(double r0, double R_0, double phi) {
  int i, N;
  double dz, z1, z2, zm, f1, f2, fm, vol;

  N = int(Nc*r0/h)+1;
  dz = r0/double(N);
  z1 = 0.0;
  f1 = area(r0-z1,r0,R_0,phi);
  vol = 0.0;

  for(i=1;i<=N;i++){
    z2 = double(i)*dz;
    zm = (z1+z2)/2.0;
    f2 = area(r0-z2,r0,R_0,phi);
    fm = area(r0-zm,r0,R_0,phi);
    vol = vol + (f1+4.0*fm+f2)/6.0*(z2-z1);
    z1 = z2;
    f1 = f2;
  }

  return vol;
}



int main() {
        int i,j,k;
	double x,y,z;
	voronoicell_neighbor c;
	vector<int> f_vert, n_list;
        vector<double> f_coord;
        vector<int>::iterator it, iv;  // declare an iterator to a vector of strings
        vector<double>::iterator ic;  // declare an iterator to a vector of strings

        double A, B, C, D, x1, x2, x3, y1, y2, y3, z1, z2, z3, r0, R_0, norm, xp, yp, zp, r23, r12, r13, cosa, phi1, phi2, ar0;
        double M, Msum, M1, Msum1, s1, s2, x0, y0, z0, Mwall, vol, density, r, R;
        int out_vert, N;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
        //particle_order po;

        ofstream myfile, myfile1;
        myfile.open ("for_3d_plotting.dat");
        myfile1.open ("mass_residuals.dat");

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}


        x0=-0.0;
        y0=-0.0;
        z0=-0.0;

        N = 5;
        k=0;

        c_loop_all clo(con);
        voronoicell_neighbor c_list[particles];

        i=0;
	if(clo.start()) do if(con.compute_cell(c,clo)) {
	      c_list[i]=c;
	      id_list[i]=clo.pid();
	      clo.pos(x,y,z);
	      x_list[i]=x;
	      y_list[i]=y;
	      z_list[i]=z;
	      //printf("%g\n", c_list[i].volume());
              i++;
	} while (clo.inc());




	for(k=0; k<particles; k++) {

		// Get the position of the current particle under consideration
	        //clo.pos(x,y,z);
	        c = c_list[k];
	        x = x_list[k];
	        y = y_list[k];
	        z = z_list[k];

	        c.face_vertices(f_vert);
	        c.vertices(f_coord);
                vol=c.volume();
                c.neighbors(n_list);

                //printf("%d -> %d\n", k, clo.pid());
                //printf("%d -> %d\n", k, id_list[k]);
 
                for(iv=n_list.begin() ; iv < n_list.end(); iv++) {
		  //printf("%d,  ", *iv); 
                }
		//printf("\n"); 

                //printf("%g\n", vol);
                //printf("%g %g %g\n", x, y, z);

                i=0;
                myfile << k+1;
                myfile << "\n";

                for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {
		  //if(k==6) printf("%d -> %g %g %g \n", i, *ic +x, *(ic+1) +y, *(ic+2) +z);  // prints d.
                  myfile << *ic +x;
                  myfile << " ";
                  i++;
                }
                myfile << "\n";

                for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {
                  myfile << *(ic+1) +y;
                  myfile << " ";
                  i++;
                }
                myfile << "\n";

                for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {
                  myfile << *(ic+2) +z;
                  myfile << " ";
                  i++;
                }
                myfile << "\n";


                //printf("---\n");

                Msum = 0;
                Msum1 = 0;

                it=f_vert.begin();
                do {
		  //printf("%d, ", *it);  // prints d.
                  //myfile << *it;
                  //myfile << "\n";

                  iv=it;
                  Mwall=0.0;
                  out_vert=0;
		  for(j=0; j<*iv; j++) {
                    it++;
                    if(it >= f_vert.end()) break;

                    if(j==0) {     // Calculating the distance from (x,y,z) to each face of a cell
                                   // http://mathinsight.org/distance_point_plane
                                   // http://mathinsight.org/forming_planes
                      x1 = f_coord[3*(*it)] +x;
                      y1 = f_coord[3*(*it)+1] +y;
                      z1 = f_coord[3*(*it)+2] +z;

	              x2 = f_coord[3*(*(it+1))] +x;
                      y2 = f_coord[3*(*(it+1))+1] +y;
                      z2 = f_coord[3*(*(it+1))+2] +z;

                      x3 = f_coord[3*(*(it+2))] +x;
		      y3 = f_coord[3*(*(it+2))+1] +y;
		      z3 = f_coord[3*(*(it+2))+2] +z;

                      A = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
                      B = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1);
                      C = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
                      D = -A*x1 -B*y1 -C*z1;

                      norm = sqrt(A*A+B*B+C*C);
                      r0 = (A*x0+B*y0+C*z0+D)/norm;
                      ar0 = abs(r0);

                      xp = x0 - r0*A/norm;
                      yp = y0 - r0*B/norm;
                      zp = z0 - r0*C/norm;

                      s2 = x1*(y2*z3-z2*y3) + y1*(z2*x3-x2*z3) + z1*(x2*y3-y2*x3);
                      //printf("face orientation: %g\n",s2);

                      //printf("The distance is: %g \n", abs(r0));
                      //if(k==6) printf("The projection is at: (%g,%g,%g) \n", xp,yp,zp);
                      //if(k==6) printf("Normal projection? -> %g \n", (x0-xp)*(x1-xp)+(y0-yp)*(y1-yp)+(z0-zp)*(z1-zp));
                      //if(k==6) printf("Normal projection? -> %g \n", (x0-xp)*(x2-xp)+(y0-yp)*(y2-yp)+(z0-zp)*(z2-zp));
                      //if(k==6) printf("Normal projection? -> %g \n", (x0-xp)*(x3-xp)+(y0-yp)*(y3-yp)+(z0-zp)*(z3-zp));
                      //printf("This should be zero: %g \n", A*xp+B*yp+C*zp+D);
                    }

                    x2 = f_coord[3*(*it)] +x;
                    y2 = f_coord[3*(*it)+1] +y;
                    z2 = f_coord[3*(*it)+2] +z;

                    if(j==(*iv-1)) {
	              x3 = x1;
                      y3 = y1;
                      z3 = z1;
                    }
                    else {
	              x3 = f_coord[3*(*(it+1))] +x;
                      y3 = f_coord[3*(*(it+1))+1] +y;
                      z3 = f_coord[3*(*(it+1))+2] +z;
		    }

                    r23 = sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2));
                    r12 = sqrt((xp-x2)*(xp-x2) + (yp-y2)*(yp-y2) + (zp-z2)*(zp-z2));
                    r13 = sqrt((xp-x3)*(xp-x3) + (yp-y3)*(yp-y3) + (zp-z3)*(zp-z3));
                    cosa = ((x3-x2)*(xp-x2) + (y3-y2)*(yp-y2) + (z3-z2)*(zp-z2)) /r12/r23;

                    R_0 = r12 * sqrt(1-cosa*cosa);
                    //printf("R0: %g\n", R_0);
                    s1 = xp*(y2*z3-z2*y3) + yp*(z2*x3-x2*z3) + zp*(x2*y3-y2*x3);
                    //if(k==6) printf("s1: %g\n", s1);

                    phi1 = acos(R_0/r12);
                    phi2 = acos(R_0/r13);

                    R = R_0/cos(phi1);
                    r = sqrt(r0*r0 + R*R);
                    if(r>=2*h) out_vert++;

                    if(s1*s2*r0 <= 0) {
                      M = -1;
                      M1 = -1;
                    }
                    else {
                      M = 1;
                      M1 = 1;
                    }


		    if((r12*sin(phi1) >= r23) || (r13*sin(phi2) >= r23)) {
		      if(phi1 >= phi2) {
                        M = M*(full_integral(phi1, ar0, R_0) - full_integral(phi2, ar0, R_0));
                        M1 = M1*(volume(ar0,R_0,phi1) - volume(ar0,R_0,phi2));
		      }
		      else {
                        M = M*(full_integral(phi2, ar0, R_0) - full_integral(phi1, ar0, R_0));
                        M1 = M1*(volume(ar0,R_0,phi2) - volume(ar0,R_0,phi1));
		      }
		    }
		    else {
                      M = M*(full_integral(phi1, ar0, R_0) + full_integral(phi2, ar0, R_0));
                      M1 = M1*(volume(ar0,R_0,phi1) + volume(ar0,R_0,phi2));
	            }

                    Msum = Msum+M;
                    Msum1 = Msum1+M1;
                    Mwall = Mwall+M;

 		    //if(k==6) printf("%d, ", *it);  // prints d.
 		    //printf("(%g, %g) ", full_integral(phi1, ar0, R_0), full_integral(phi2, ar0, R_0));  
                    //if(k==6) printf("(%g %g) ", r0,M);
                    myfile << *it;
                    myfile << " ";
		  }
                 
                  //if(out_vert==*iv) printf("(%d is outside) ", out_vert);
		  it++;
		  //if(k==6) printf("-> %g\n\n", Mwall);
		  //printf("\n");
                  myfile << "\n";
                } while(it < f_vert.end());

		// printf("\n\n");
		//printf("Msum: %g \n", Msum);  
		printf("Msum: %g %g\n", Msum, Msum1);  
                myfile << abs(Msum)/vol/5.;
                myfile << "\n";

                myfile1 << Msum1 << " " << Msum;
                myfile1 << "\n";

                //if(k==1) break;

	}

        for(i=0;i<=10;i++){
          //printf("Test the integral: %g\n", full_integral(pi/4.,0.2*i*h, 0.2*i*h)*2.*4.*6.);
        }

        myfile << "0";
        myfile.close();
        myfile1.close();

}



/*
int main() {
  double phi, r0, R_0;
  double mu, a, invtan, u, I1;
  double a2, cosp2;
  int i, j, k;

  h=0.25;
  r0=0.1*h;
  R_0=1.0*h;
  phi=0.0;

  printf("(%g, %g)\n", full_integral(phi, r0, R_0), full_integral_new(phi, r0, R_0));

  for(i=0;i<10;i++) {
    for(j=0;j<10;j++) {
      for(k=0;k<10;k++) {
	r0 = (3.0/10.0)*i*h;
	R_0 = (3.0/10.0)*j*h;
	phi = (pi/2.0/10.0)*k;

	a = R_0/r0;
	mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

	a2 = pow(a,2);
	cosp2 = pow(cos(phi),2);
	
	I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );
	
	u = sin(phi)*sqrt(1-mu*mu);
	invtan = atan(u/a);

	//printf("(%g %g %g) -> %g %g\n",r0,R_0,phi,I1+pi/2.,invtan);
	//printf("(%g %g %g) -> %g %g\n",r0,R_0,phi,r0*R_0*R_0*tan(phi)/6.0, volume(r0,R_0,phi,10));
	//printf("(%g %g %g) -> %g \n",r0,R_0,phi, full_integral_new(phi, r0, R_0) - volume(r0,R_0,phi,30));
      }
    }
  }

  //phi = 1.41861;
  h=1.5;
  r0=0.3*h;
  R_0=0.3*h;

  for(i=1;i<=20;i++) {
    phi = pi/2./20.*double(i-1);
    printf("(%g %g %g) -> %g %g %g\n\n",r0/h,R_0/h,phi, full_integral(phi, r0, R_0), full_integral_new(phi, r0, R_0), volume(r0,R_0,phi,20));
  }
  return 0;
}


*/

