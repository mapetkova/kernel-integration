// g++ cell_density_BFS_data.cc -o cell_density_BFS_data -I/usr/local/include/voro++/ -L/usr/local/lib -lvoro++

// Author: Maya A. Petkova (map32@st-andrews.ac.uk)

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <ctime>
using namespace std;

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry: Cartesian box to hold all the particles
const double x_min=-12.06, x_max=15.72;
const double y_min=-9.45, y_max=9.05;
const double z_min=-8.37, z_max=11.45;

const double pi=M_PI;
double h;
const double Nc=10.0;

// Set up the number of blocks that the container is divided into
int n_x=100,n_y=100,n_z=100;

// Set the number of particles that are going to be randomly introduced
const int particles=400728;

int id_list[particles], reverse_id_list[particles], kernel_list[particles], next_list[particles];
double x_list[particles], y_list[particles], z_list[particles], h_list[particles];
double x_centroid[particles], y_centroid[particles], z_centroid[particles], density_centroid[particles];
voronoicell_neighbor c_list[particles];
int next_list_top, next_list_bottom, jg;


// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}



double full_integral_new(double phi, double r0, double R_0) {

  double B1, B2, B3, mu, a, logs, u;
  double full_int;
  double a2, cosp, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3, tanp;
  double r2, R, linedist2, phi1, phi2, h2;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double D2, D3;

  if(r0==0.0) return 0.0;
  if(R_0==0.0) return 0.0;
  if(phi==0.0) return 0.0;

  h2 = h*h;
  r03 = r0*r0*r0;
  r0h2 = r0/h*r0/h;
  r0h3 = r0h2*r0/h;
  r0h_2 = h/r0*h/r0;
  r0h_3 = r0h_2*h/r0;


  if(r0 >= 2.0*h) {
    B3 = h2*h /4.;
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

  linedist2 = r0*r0 + R_0*R_0;
  R = R_0/cos(phi);
  r2 = r0*r0 + R*R;


  full_int = 0.0;
  D2 = 0.0;
  D3 = 0.0;

  if(linedist2 <= h2) {
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
  else if(linedist2 <= 4.0*h2) {
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

  
  if(r2 < h2) {
    full_int = r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);
  }
  else if(r2 < 4.0*h2) {
    full_int=  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
					  1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2);
  }
  else {
    full_int = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3);
  }
  
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
  int i, N;
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






double compute_mass(int ind, int k, double x0, double y0, double z0){
  // ind      - current particle index (starting from 1)
  // k        - current cell index (starting from 0)
  // x0,y0,z0 - current particle coordinates

	double x,y,z;
	voronoicell_neighbor c;
	vector<int> f_vert, n_list;
        vector<double> f_coord;
        vector<int>::iterator it, iv;  
        vector<double>::iterator ic;  

        double A, B, C, D, x1, x2, x3, y1, y2, y3, z1, z2, z3, r0, R_0, norm, xp, yp, zp, r23, r12, r13, cosa, phi1, phi2, ar0;
        double M, Msum, s1, s2, Mwall, vol, density, r, R;
        int out_vert, j, i, cell_ind, r_cell_ind, N;


	c = c_list[k];
	x = x_list[k];
	y = y_list[k];
	z = z_list[k];

	c.face_vertices(f_vert);
	c.vertices(f_coord);
	vol=c.volume();
	c.neighbors(n_list);

        N=12;
	Msum = 0;
        i=0;

	it=f_vert.begin();

	do {
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

	    if(abs(cosa)<1.0){
	      R_0 = r12 * sqrt(1-cosa*cosa);
	    }
	    else {
	      if(abs(cosa)-1.0 < 0.00001){
		R_0=0.0;
	      }
	      else {
                printf("Error: cosa > 1: %g\n", cosa);
	      }
	    }

	    s1 = xp*(y2*z3-z2*y3) + yp*(z2*x3-x2*z3) + zp*(x2*y3-y2*x3);
	    //printf("s1: %g\n", s1);

	    if(R_0<r12){
              phi1 = acos(R_0/r12);
            }
            else {
              if(R_0-r12 < 0.00001){
                phi1=0.0;
              }
              else {
                printf("Error: R0 > r12: %g\n", R_0-r12);
              }
            }
	    if(R_0<r13){
  	      phi2 = acos(R_0/r13);
            }
            else {
              if(R_0-r13 < 0.00001){
                phi2=0.0;
              }
              else {
                printf("Error: R0 > r13: %g\n", R_0-r13);
              }
            }
	    
	    R = R_0/cos(phi1);
	    r = sqrt(r0*r0 + R*R);
	    if(r>=2*h) out_vert++;
	    
	    if(s1*s2*r0 <= 0) {
	      M = -1;
	    }
	    else {
	      M = 1;
	    }

	    if((r12*sin(phi1) >= r23) || (r13*sin(phi2) >= r23)) {
	      if(phi1 >= phi2) {
		M = M*(full_integral_new(phi1, ar0, R_0) - full_integral_new(phi2, ar0, R_0));
		//M = M*(volume(ar0, R_0, phi1) - volume(ar0, R_0, phi2));
	      }
	      else {
		M = M*(full_integral_new(phi2, ar0, R_0) - full_integral_new(phi1, ar0, R_0));
		//M = M*(volume(ar0, R_0, phi2) - volume(ar0, R_0, phi1));
	      }
	    }
	    else {
	      M = M*(full_integral_new(phi1, ar0, R_0) + full_integral_new(phi2, ar0, R_0));
	      //M = M*(volume(ar0, R_0, phi1) + volume(ar0, R_0, phi2));
	    }
	    
	    Msum = Msum+M;
	  }
                 
	  if(out_vert!=*iv) {
            cell_ind = *(n_list.begin()+i);
            r_cell_ind = reverse_id_list[cell_ind];
            if((cell_ind>=0) && (kernel_list[r_cell_ind]!=ind)) {
	      next_list[next_list_top] = r_cell_ind;
              kernel_list[r_cell_ind]=ind;
	      next_list_top++;
	    }
	  }
	  it++;
          i++;

	} while(it < f_vert.end());


	return Msum;
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
        double M, Msum, s1, x0, y0, z0, Mwall, vol, density, r, R, Mtot,Mtot_centroid;
        int out_vert, ind;
        double Mass[particles];

        clock_t begint = clock();

	// Set up pre-container and importing particle postions from a file
        pre_container pcon(x_min,x_max,y_min,y_max,z_min,z_max,false,false,false);
        pcon.import("sphtest6106_voronoi_grid_part-400728.dat");
        pcon.guess_optimal(n_x,n_y,n_z);
  
        // Set up the container class and import the particles from the
        // pre-container  --> I wish I knew why this was necessary!
        container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
        pcon.setup(con);
 
        printf("Not broken yet ;) 1\n");


        //std::ifstream infile("sphtest6106_voronoi_grid_h.dat");
        std::ifstream infile("sphtest6106_voronoi_grid_h-400728.dat");
        i=0;
        while (infile >> h)
        {
	  h_list[i]=h;
	  i++;
        }
        infile.close();

        c_loop_all clo(con);

        clock_t endt = clock();
        printf("Importing data took: %g sec\n", double(endt - begint) / CLOCKS_PER_SEC);

        i=0;

	if(clo.start()) do if(con.compute_cell(c,clo)) {
	      c_list[i]=c;
	      id_list[i]=clo.pid();
              reverse_id_list[id_list[i]]=i;

	      clo.pos(x,y,z);
	      x_list[i]=x;
	      y_list[i]=y;
	      z_list[i]=z;

              Mass[i]=0.0;
              kernel_list[i]=0;
	      density_centroid[i]=0.0;

	      c.centroid(x,y,z);
	      x_centroid[i]=x_list[i]+x;
	      y_centroid[i]=y_list[i]+y;
	      z_centroid[i]=z_list[i]+z;

              i++;
	} while (clo.inc());


        begint=endt;
        endt = clock();
        printf("Setup for densities took: %g sec\n", double(endt - begint) / CLOCKS_PER_SEC);

        printf("Not broken yet ;) 3\n");

        jg=1;

        for(ind=1;ind<=particles;ind++) {
          x0=x_list[ind-1];
          y0=y_list[ind-1];
          z0=z_list[ind-1];

          k=ind-1;
          h=h_list[id_list[ind-1]-1];

	  next_list_bottom = 0;
	  next_list[next_list_bottom] = k;
	  kernel_list[next_list[next_list_bottom]]=ind;
	  next_list_top = 1;

	  while(next_list_bottom < next_list_top) {
	    i = next_list[next_list_bottom];
	    Mass[i] = Mass[i] + compute_mass(ind,i,x0,y0,z0);

	    density_centroid[i] = density_centroid[i] + kernel(y_centroid[i]-y0, x_centroid[i]-x0, z_centroid[i]-z0);

            next_list_bottom++;
	  }
	  //printf("%d\n", next_list_top);
	}

        begint=endt;
        endt = clock();
        printf("Calculating densities took: %g sec\n", double(endt - begint) / CLOCKS_PER_SEC);
        printf("Not broken yet ;) 4\n");

	////////////// Output for plotting ///////////////

        ofstream myfile, statfile;
        myfile.open ("for_3d_plotting.dat");
        statfile.open ("density_comparison.dat");

        Mtot=0.0;
	Mtot_centroid=0.0;

        for(i=0;i<particles;i++) {
	  c = c_list[i];
	  x = x_list[i];
	  y = y_list[i];
	  z = z_list[i];

	  c.face_vertices(f_vert);
	  c.vertices(f_coord);
	  vol=c.volume();

          statfile << x;
          statfile << " ";
          statfile << y;
          statfile << " ";

	  myfile << i+1;
	  myfile << "\n";

	  for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {
	    myfile << *ic +x;
	    myfile << " ";
	  }
	  myfile << "\n";
	  
	  for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {
	    myfile << *(ic+1) +y;
	    myfile << " ";
	  }
	  myfile << "\n";
	  
	  for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {
	    myfile << *(ic+2) +z;
	    myfile << " ";
	  }
	  myfile << "\n";



	  it=f_vert.begin();
	  do {
	    iv=it;
	    for(j=0; j<*iv; j++) {
	      it++;
	      if(it >= f_vert.end()) break;
	      
	      myfile << *it;
	      myfile << " ";
	    }
	    it++;
	    myfile << "\n";
	  } while(it < f_vert.end());
	  
	  myfile << abs(Mass[i])/vol/2000000.;
	  myfile << "\n";
	  
          //statfile << (Mass[i]-1.0)*100.0;
          statfile << (Mass[i]/vol-density_centroid[i])/(Mass[i]/vol)*100.0;
	  statfile << "\n";

	  //printf("%g\n", Mass[i]);
	  //printf("%g\n", (Mass[i]-1.0)*100.0);
          Mtot=Mtot+Mass[i];
	  Mtot_centroid = Mtot_centroid + density_centroid[i]*vol;
	}
        printf("Mtot: %g\n", Mtot);
        printf("Mtot_centroid: %g\n", Mtot_centroid);

        myfile << "0";
        myfile.close();

        statfile.close();


        begint=endt;
        endt = clock();
        printf("Outputting data took: %g sec\n", double(endt - begint) / CLOCKS_PER_SEC);

	return 0;
}

