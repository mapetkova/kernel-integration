/**
 * A simple program to calculate the mass of an SPH particle enclosed by a cube that
 * is centered around it. 
 *
 * To compile use something like: g++ density_function.cc -o density_function
 * Author: Maya A. Petkova (map32@st-andrews.ac.uk)
*/
#include <iostream>
#include <math.h>

double full_integral(double phi, double r0, double R_0, double h);

int main() {
  double phi, r0, R_0, h;
  double cell_size;

  h=5.0;
  cell_size=4*h;
  r0=cell_size/2.0;
  R_0=cell_size/2.0;
  phi=M_PI/4.0;

  /**
   * full_integral(phi, r0, R_0, h) is the expression evaluated at each vertex due to symmetry.
   * We multiply it by 2 because each edge contains 2 vertices.
   * We multiply it further by 4 because each face has 4 edges.
   * We multiply it further by 6 because a cube has 6 faces.
   */
  printf("A cube with side %g times h, and with an SPH particle at its centre\n", cell_size/h);
  printf("encloses the following fraction of the mass of the particle: %g\n", 6*4*2*full_integral(phi, r0, R_0, h));

  return 0;
}


/**
 * Function that gives the 3D integral of the kernel of
 * a particle for a given vertex of a cell face.
 *
 * phi: Azimuthal angle of the vertex.
 * r0: Distance from the particle to the face of the cell.
 * R_0: Distance from the orthogonal projection of the particle
 * onto the face of the cell to a side of the face (containing the vertex).
 * h: The kernel smoothing length of the particle.
 */
double full_integral(double phi, double r0, double R_0, double h) {

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

  // Setting up the B1, B2, B3 constants of integration.

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

  if(linedist < 1.0*h) {
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
  else if(linedist < 2.0*h) {
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
  
  
  //////////////////////////////////
  // Calculating I_n expressions. //
  //////////////////////////////////


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

  // Calculating the integral expression.

  if(r < 1.0*h) {
    full_int = r0h3/M_PI  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);
  }
  else {
    if(r < 2.0*h) {
      full_int=  r0h3/M_PI  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
				     1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2);
    }
    else {
      full_int = r0h3/M_PI  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3);
    }
  }  

  return full_int;
}
