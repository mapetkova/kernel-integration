/**
 * A program to calculate the cell density of a Vaoronoi grid with
 * a single SPH particle mapped on it. 
 *
 * To compile use something like: g++ single_particle_voronoi_density.cc -o single_particle_voronoi_density -I/usr/local/include/voro++/ -L/usr/local/lib -lvoro++
 * The output can be displayed with density_structure_3d.py
 *
 * Author: Maya A. Petkova (map32@st-andrews.ac.uk)
*/
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

#include "voro++.hh"
using namespace voro;


double vertex_integral(double phi, double r0, double R_0, double h);

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
  // Set up constants for the container geometry
  const double x_min=-1,x_max=1;
  const double y_min=-1,y_max=1;
  const double z_min=-1,z_max=1;
  // Set up the number of blocks that the container is divided into
  const int n_x=6,n_y=6,n_z=6;
  // Set the number of cell generating sites that are going to be randomly introduced
  const int cell_num=50;

  // Set up smoothing length and position of the particle
  const double h=0.5;
  const double x0=0.0;
  const double y0=0.0;
  const double z0=0.0;

  int id_list[cell_num];
  double x_list[cell_num], y_list[cell_num], z_list[cell_num];
  voronoicell_neighbor c_list[cell_num];

  int i,j,k;
  double x,y,z;
  voronoicell_neighbor c;
  vector<int> f_vert, n_list;
  vector<double> f_coord;
  vector<int>::iterator it, iv; 
  vector<double>::iterator ic;  

  double A, B, C, D, x1, x2, x3, y1, y2, y3, z1, z2, z3, xp, yp, zp, r23, r12, r13, norm;
  double r0, R_0, cosa, phi1, phi2, ar0, r, R;
  double M, Msum, s1, s2, Mtot, vol, density;


  // Create a container with the geometry given above, and make it
  // non-periodic in each of the three coordinates. Allocate space for
  // eight particles within each computational block
  container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);

  // Randomly add cell generating sites into the container
  for(i=0;i<cell_num;i++) {
    x=x_min+rnd()*(x_max-x_min);
    y=y_min+rnd()*(y_max-y_min);
    z=z_min+rnd()*(z_max-z_min);
    con.put(i,x,y,z);
  }

  // Loop over the cell generating sites and record them in the order 
  // that they appear in the container (the container shuffles them)
  c_loop_all clo(con);
  i=0;
  if(clo.start()) do if(con.compute_cell(c,clo)) {
    c_list[i]=c;
    id_list[i]=clo.pid();
    clo.pos(x,y,z);
    x_list[i]=x;
    y_list[i]=y;
    z_list[i]=z;
    i++;
  } while (clo.inc());


  // Open the output file for 3D plotting
  ofstream myfile;
  myfile.open ("for_3d_plotting.dat");

  Mtot = 0;

  // Loop over all cells
  for(k=0; k<cell_num; k++) {
    // Get the position of the current cell under consideration
    c = c_list[k];
    x = x_list[k];
    y = y_list[k];
    z = z_list[k];

    c.face_vertices(f_vert);
    c.vertices(f_coord);
    vol=c.volume();
    c.neighbors(n_list);

    i=0;

    myfile << k+1 << "\n";
    for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {myfile << *ic +x << " ";}
    myfile << "\n";
    for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {myfile << *(ic+1) +y << " ";}
    myfile << "\n";
    for(ic=f_coord.begin() ; ic < f_coord.end(); ic+=3) {myfile << *(ic+2) +z << " ";}
    myfile << "\n";

    Msum = 0;

    // Loop over the vertices of each face
    it=f_vert.begin();
    do {
      iv=it;
      for(j=0; j<*iv; j++) {
        it++;
        if(it >= f_vert.end()) break;

        if(j==0) {     // Calculate the distance from (x,y,z) to each face of a cell (i.e. r0)
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

          // The orthogonal projection of the particle position on the face is at (xp,yp,zp) 
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
        } else {
          x3 = f_coord[3*(*(it+1))] +x;
          y3 = f_coord[3*(*(it+1))+1] +y;
          z3 = f_coord[3*(*(it+1))+2] +z;
        }

        // Calculate the distance from (xp,yp,zp) to each edge of the face (i.e. R_0)
        r23 = sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2));
        r12 = sqrt((xp-x2)*(xp-x2) + (yp-y2)*(yp-y2) + (zp-z2)*(zp-z2));
        r13 = sqrt((xp-x3)*(xp-x3) + (yp-y3)*(yp-y3) + (zp-z3)*(zp-z3));
        cosa = ((x3-x2)*(xp-x2) + (y3-y2)*(yp-y2) + (z3-z2)*(zp-z2)) /r12/r23;

        R_0 = r12 * sqrt(1-cosa*cosa);
        s1 = xp*(y2*z3-z2*y3) + yp*(z2*x3-x2*z3) + zp*(x2*y3-y2*x3);

        phi1 = acos(R_0/r12);
        phi2 = acos(R_0/r13);

        R = R_0/cos(phi1);
        r = sqrt(r0*r0 + R*R);

        // Establish the sign of the edge contribution to the mass of the cell
        if(s1*s2*r0 <= 0) {
          M = -1;
        } else {
          M = 1;
        }

        // Different vertex contributions depending on where (xp,yp,zp) is projected 
        // on the line of the edge 
        if((r12*sin(phi1) >= r23) || (r13*sin(phi2) >= r23)) {
          if(phi1 >= phi2) {
            M = M*(vertex_integral(phi1, ar0, R_0, h) - vertex_integral(phi2, ar0, R_0, h));
          } else {
            M = M*(vertex_integral(phi2, ar0, R_0, h) - vertex_integral(phi1, ar0, R_0, h));
          }
        } else {
          M = M*(vertex_integral(phi1, ar0, R_0, h) + vertex_integral(phi2, ar0, R_0, h));
        }

        Msum = Msum+M;
        myfile << *it << " ";
      }
      it++;
      myfile << "\n";
    } while(it < f_vert.end());

    if(abs(Msum)<0.0000001) Msum = 0;
    printf("Mass of cell %d: %g \n", k, Msum);  
    myfile << abs(Msum)/vol/5. << "\n";

    Mtot = Mtot+Msum;
  }

  printf("Total cell mass: %g (it should be 1 for h<=0.5)\n", Mtot);  

  myfile << "0";
  myfile.close();

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
double vertex_integral(double phi, double r0, double R_0, double h) {

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
