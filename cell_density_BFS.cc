// g++ cell_density_BFS.cc -o cell_density_BFS -I/usr/local/include/voro++/ -L/usr/local/lib -lvoro++

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

const double pi=M_PI;
double h;

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=2000;

int id_list[particles], reverse_id_list[particles], kernel_list[particles], next_list[particles];
double x_list[particles], y_list[particles], z_list[particles];
voronoicell_neighbor c_list[particles];
int next_list_top, next_list_bottom;

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}


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

  I_1 = a/2.*logs + invtan;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-pow(u,2)) + logs);
  I_5 = I_3 + a*pow((1.+a2),2)/16. *( (10*u - 6*pow(u,3))/pow(1-pow(u,2),2) + 3.*logs);

  integral =  r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);

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

  I_1 = a/2.*logs + invtan;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-pow(u,2)) + logs);
  I_5 = I_3 + a*pow((1.+a2),2)/16. *( (10*u - 6*pow(u,3))/pow(1-pow(u,2),2) + 3.*logs);

  integral =  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
					  1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0);

  return integral;
}

double F3_3d(double phi, double r0, double R_0, double B3) {

  double integral;
  double I0, I1;
  double a, a2, cosp2, r03, r0h3, r0h_3;

  a = R_0/r0;

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h3 = pow(r0/h,3);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  integral = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0);

  return integral;
}


double full_integral(double phi, double r0, double R_0) {

  double B1, B2, B3, mu, a, logs, invtan, u;
  double full_int;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;
  double r, R, linedist, phi1, phi2;

  if(r0==0.0) return 0.0;

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



  if(r0/mu < 1.0*h) {
    full_int = F1_3d(phi,r0,R_0,B1);
  }
  else if(r0/mu < 2.0*h) {
    full_int = F2_3d(phi,r0,R_0,B2);
  }
  else {
    full_int = F3_3d(phi,r0,R_0,B3);
  }



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



double compute_mass(int ind, int k, double x0, double y0, double z0){

	double x,y,z;
	voronoicell_neighbor c;
	vector<int> f_vert, n_list;
        vector<double> f_coord;
        vector<int>::iterator it, iv;  // declare an iterator to a vector of strings
        vector<double>::iterator ic;  // declare an iterator to a vector of strings

        double A, B, C, D, x1, x2, x3, y1, y2, y3, z1, z2, z3, r0, R_0, norm, xp, yp, zp, r23, r12, r13, cosa, phi1, phi2, ar0;
        double M, Msum, s1, s2, Mwall, vol, density, r, R;
        int out_vert, j, i, cell_ind, r_cell_ind;


	c = c_list[k];
	x = x_list[k];
	y = y_list[k];
	z = z_list[k];

	c.face_vertices(f_vert);
	c.vertices(f_coord);
	vol=c.volume();
	c.neighbors(n_list);

	//printf("%d -> %d\n", k, id_list[k]);
 
	for(iv=n_list.begin() ; iv < n_list.end(); iv++) {
	  //printf("%d,  ", *iv); 
	}
	//printf("\n"); 

	//printf("%g\n", vol);
	//printf("%g %g %g\n", x, y, z);

	Msum = 0;
        i=0;

	it=f_vert.begin();
	do {
	  //printf("%d, ", *it); 

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

	    R_0 = r12 * sqrt(1-cosa*cosa);
	    //printf("R0: %g\n", R_0);
	    s1 = xp*(y2*z3-z2*y3) + yp*(z2*x3-x2*z3) + zp*(x2*y3-y2*x3);
	    //printf("s1: %g\n", s1);

	    phi1 = acos(R_0/r12);
	    phi2 = acos(R_0/r13);
	    
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
		M = M*(full_integral(phi1, ar0, R_0) - full_integral(phi2, ar0, R_0));
	      }
	      else {
		M = M*(full_integral(phi2, ar0, R_0) - full_integral(phi1, ar0, R_0));
	      }
	    }
	    else {
	      M = M*(full_integral(phi1, ar0, R_0) + full_integral(phi2, ar0, R_0));
	    }
	    
	    Msum = Msum+M;
	    Mwall = Mwall+M;
	    
	    //printf("%d, ", *it);  // prints d.
	    //printf("(%g, %g) ", full_integral(phi1, ar0, R_0), full_integral(phi2, ar0, R_0));  
	    //printf("(%g %g) ", r0,M);
	  }
                 
	  if(out_vert==*iv) {
            //printf("(%d is outside) ", out_vert);
	  }
	  else {
            cell_ind = *(n_list.begin()+i);
            r_cell_ind = reverse_id_list[cell_ind];
            //printf("(%d) ", r_cell_ind);
            if((cell_ind>=0) && (kernel_list[r_cell_ind]!=ind)) {
	      next_list[next_list_top] = r_cell_ind;
              kernel_list[r_cell_ind]=ind;
	      next_list_top++;
	    }
	  }
	  it++;
          i++;

	} while(it < f_vert.end());

	//printf("\n\n");
	//printf("Msum: %g \n", Msum);  

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
        double M, Msum, s1, x0, y0, z0, Mwall, vol, density, r, R, Mtot;
        int out_vert, ind;
        double Mass[particles];

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}



        c_loop_all clo(con);

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
              i++;
	} while (clo.inc());



        for(ind=1;ind<=particles;ind++) {
          x0=x_list[ind-1];
          y0=y_list[ind-1];
          z0=z_list[ind-1];

	  //x0=-1.5 + 1.0*ind;
	  //y0=-1.5 + 1.0*ind;
	  //z0=-1.5 + 1.0*ind;

	  //x0=-0.1;
	  //y0=-0.1;
	  //z0=-0.1;

	  //k=15;
	  //ind=1;
          k=ind-1;

          //h=0.01;
	  //if((1-abs(x0))/2 < 0.01) {h=(1-abs(x0))/2; printf("%d, h = %g\n", ind, h);}
	  //if((1-abs(y0))/2 < 0.01) {h=(1-abs(y0))/2; printf("%d, h = %g\n", ind, h);}
	  //if((1-abs(z0))/2 < 0.01) {h=(1-abs(z0))/2; printf("%d, h = %g\n", ind, h);}

          h=rnd()/2;

	  next_list_bottom = 0;
	  next_list[next_list_bottom] = k;
	  kernel_list[next_list[next_list_bottom]]=ind;
	  next_list_top = 1;

	  while(next_list_bottom < next_list_top) {
	    i = next_list[next_list_bottom];
	    Mass[i] = Mass[i] + compute_mass(ind,i,x0,y0,z0);
	    //Mass[i] = compute_mass(ind,i,x0,y0,z0);
            next_list_bottom++;
	  }
	}


	////////////// Output for plotting ///////////////

        ofstream myfile;
        myfile.open ("for_3d_plotting.dat");

        Mtot=0.0;

        for(i=0;i<particles;i++) {
	  c = c_list[i];
	  x = x_list[i];
	  y = y_list[i];
	  z = z_list[i];

	  c.face_vertices(f_vert);
	  c.vertices(f_coord);
	  vol=c.volume();

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
	  
	  myfile << abs(Mass[i])/vol/100.;
	  myfile << "\n";
	  
	  
	  printf("%g\n", Mass[i]);
          Mtot=Mtot+Mass[i];
	}
        printf("Mtot: %g\n", Mtot);

        myfile << "0";
        myfile.close();

	return 0;
}

