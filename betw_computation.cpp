//To update this file and the others in the same folder on niko:
//rsync -rv -e 'ssh -p 8022' ~/cpp/code_Betweenness/*.* abaudena@localhost:~/cpp/code_Betweenness/

#include <iostream>
#include <iomanip>
#include <ctime>  
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <omp.h>

using namespace std;
#include "parameters.h" // Function to read parameters.dat
#include "read_velocity.h" // Function to read velocities 

/////////////////////////////////////////////////////////////////////////////////////////////////
/////// Prototype-functions used in the main code (written below, at the end of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

void print_usage(const string me);
int countlines(string namefile);
string leadingzeros(unsigned int ndigits, double number);
string leadingzeros_int(unsigned int ndigits, int number);

// boundary condition matching function (tracer long lat)
vector<double> boundary_cond(vector<double> grib_tr, int &pass_WE_YN);

vector<double> velocity_xy_interp(double cord_t, vector<double> cord_xy_vec);

// rk4 function
vector<double> rungeKutta4(double start_cord_t,
                           vector<double> start_cord_xy,
                           vector<double> velocity_xy_interp(double cord_t,
							     vector<double> cord_xy_vec));
// largest eigenvalue computation function
double eigenvalue(vector<double> Delta);	
// ftle computation function
double ftle_func(vector<double> x0, vector<double> xf, double t_ftle,double eigenvalue(vector<double> Delta));
/////////////////////////////////////////////////////////////////////////////////////////////////
// Declaring constant and global variables
/////////////////////////////////////////////////////////////////////////////////////////////////

const double pig = 3.14159265358; 
const double deg2rad = (pig/180.); // Conversion from degree to radians
const double day2sec = 60. * 60. * 24.; // Daily data -> velocities converted from seconds to days
const double Radius = 6372795.48;     // radius of the Earth (m)

const string input_dir = "/home/abaudena/cpp/code_Betweenness/Outputs/";
const string output_dir = "/home/abaudena/cpp/code_Betweenness/Outputs/";


// Tracer domain
/*AB_s double start_tracer_long;
double start_tracer_lat;
double end_tracer_long; 
double end_tracer_lat;  *///AB_e

/////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CODE 
/////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  // READ COMMAND LINE ARGUMENTS:
  // Usage in command line when executing: ./(path executable) (path to parameters file) 
  // (ex. ./grid_construction.out parameters.dat)
  if (2 != argc) 
    {
      // If it is not specified the path to parameters file
      cout << "Parameter file not specified" << endl; // Show a error message  
      string me=argv[0];
      print_usage(me); // and a usage reminder (me = name of executable) 
      return 1;
    }

  // READ PARAMETERS FROM FILE
  string fnameparameters = argv[1]; // Path to parameters file = Second argument in the command line 
  if(readparams(fnameparameters)==1)
    {
      // If it is not specified the path to parameters file
      cout << "Error in reading file parameters" << endl;// Show a error message 
      return 1;
    }

// File names postfixs
string rand_str= "_" + leadingzeros_int(8,rand_idx);
string postfixgrid = rand_str + ".grid";
string postfixtracer = rand_str + ".trac";
string postfixbetw = rand_str + ".betw";

  //////////////////// READ TRACERS INITIAL POSITIONS //////////////////

//AB  string sdomain_network;
//AB  string svel_product;
//AB  string slayer_index;
  string sparticle_latspacing;
  string sdelta_betw;
  string sstart_date;
  string stau;
  string sfile_initial_tracers;
  string sfile_betw;
  string sint_step;
  char buffer[50];

/*AB_s  sprintf(buffer, "networkdomain%d", domain_network);
  sdomain_network = buffer;
  sprintf(buffer, "_vel_product%d", vel_product);
  svel_product=buffer;*///AB_e
//AB  sprintf(buffer, "_depth%d", layer_index);
//AB  slayer_index=buffer;
  sparticle_latspacing = "_particlelatspacing";

  // Reconstruct name file
  sfile_initial_tracers = 
    input_dir + 
    prod_name +
//AB    slayer_index + 
//AB    sparticle_latspacing + leadingzeros(6,particle_latspacing) + 
    postfixtracer; 

  // Count number of lines = total number of lagrangian particles
  int num_itracers=countlines(sfile_initial_tracers);
  if(num_itracers < 0)
    {
      cout << "ERROR in counting lines" << endl;
      return 1;
    }

  // Print initial number of tracers to screen for debugging
  if(verbose == 1) 
    {
      cout <<" Initial number of tracers = "<< num_itracers << endl;
    }

  double *itracer_lon;
  double *itracer_lat;

  itracer_lon = new double [num_itracers];
  itracer_lat = new double [num_itracers];

  // open tracer file 
  ifstream file_initial_tracers(sfile_initial_tracers.c_str());// open tracer file 

  int i;

  // Read initial positions of particles 
  for(i=0; i<num_itracers; i++)
    {
      file_initial_tracers >> itracer_lat[i] >> itracer_lon[i];
    }
  file_initial_tracers.close();
  
  // READ layer numbered "layer_index" of the VELOCITY FIELD (from initial day "start_date" to "start_date + tau + 1")
  //AB_s
  int n_velfield2load;
  n_velfield2load=tau/vel_field_timestep+2;//+2 cos I add the starting time step and a final one, after tau, which is necessary for the temporal interpolation in velocity_xy_interp
//  cout<<"n_velfield2load "<<n_velfield2load<< endl;
  //AB_e
  //AB_s: computing betweenness backward of tau days is equivalent to computing it forward starting from start_date-tau. This cos betw is simmetric in time
  sprintf(buffer, "_tau%d", (int) tau*time_dir);
  stau=buffer;
  sstart_date="_startdate";
  sstart_date=sstart_date + leadingzeros_int(2,start_date.tm_mday) + leadingzeros_int(2,start_date.tm_mon+1) + leadingzeros_int(4,start_date.tm_year) + "_h" + leadingzeros_int(2,start_date.tm_hour)  ;
  if (time_dir==-1)
  {
  start_date.tm_mday= start_date.tm_mday-tau;
  mktime(&start_date);
  time_dir=1;
  }
  //AB_e
  if(read_vfield(start_date, n_velfield2load,time_dir)!=0)//AB (removed layer_index)
    {
      cout << "Error in reading velocities"<< endl;
      return 1;
    }
    
  
 //////////////////// BETW COMPUTATION //////////////////
  
  	double t,t_ftle;
  	vector<double> betw(num_itracers),ftle_forw((int)tau+1),ftle_back((int)tau+1);
  	for ( i=0;i<num_itracers;i++) { betw[i]=0;} // I initialise the betw vector
    int t0,j1,num_par=5;//how many particle I compute the Ftle with
	vector<int> n_pass_WE(num_par); //n_pass_WE: how many times each particle passed the WE line    
   // Parallelizing the code for computing trajectories 

/*   int number_proc = omp_get_num_procs();
   printf("There are %d processors\n", number_proc);
   int threads_used = 4; // Number of cores I use
   printf("We have %d threads\n", threads_used);*/

//#pragma omp parallel for default(shared) private(t) num_threads(threads_used)


   for(i=0; i<num_itracers; i++)
    {
     // Variable to overwrite the successive positions of each particle
     vector<double> x00(2),x0(10),xf(10);
	 double delta0_ftle_lon; 
	 delta0_ftle_lon=delta0_ftle/cos(x00[1]*deg2rad);

     x00[0] = itracer_lon[i];
     x00[1] = itracer_lat[i];
      
     x0[0]=x00[0];
     x0[1]=x00[1];
     n_pass_WE[0]=0;

     x0[2]=x00[0]+delta0_ftle_lon;
     x0[3]=x00[1];

     x0[4]=x00[0];
     x0[5]=x00[1]+delta0_ftle;
    
     x0[6]=x00[0]-delta0_ftle_lon;
     x0[7]=x00[1];

     x0[8]=x00[0];
     x0[9]=x00[1]-delta0_ftle;
     
     xf[0]=x0[0];//I don't know if these 2 lines are really useful
     xf[1]=x0[1];

	 // I start the loop on the days, from 0 to tau
	 for(t0=0; t0<=tau; t0++)
	 {

	 // I do the forward advection
	 	time_dir=1;
     	for(j1=1; j1<num_par; j1++)
     	{
     		x00[0] = x0[j1*2];
     		x00[1] = x0[j1*2+1];
     		n_pass_WE[j1]=0;
     		x00 = boundary_cond(x00,n_pass_WE[j1]);

     		for(t=t0; t<tau; t+=int_step)
	 		{	    
	  		x00 = rungeKutta4(t, x00, velocity_xy_interp);
	  		x00 = boundary_cond(x00,n_pass_WE[j1]);
      		}
      		xf[j1*2] = x00[0] + 360*n_pass_WE[j1];
     		xf[j1*2+1] = x00[1];	
    	}
	 	// I compute forward Ftle
	 	t_ftle=tau-(double)t0;
	 	ftle_forw[(int)t0]= ftle_func(x0,xf,t_ftle,eigenvalue);
	 	
	 	
	 	// I do the backward advection
	 	time_dir=-1;
     	for(j1=1; j1<num_par; j1++)
     	{
     		x00[0] = x0[j1*2];
     		x00[1] = x0[j1*2+1];
     		n_pass_WE[j1]=0;
     		x00 = boundary_cond(x00,n_pass_WE[j1]);

     		for(t=t0; t>0; t-=int_step)
	 		{	    
	  		x00 = rungeKutta4(t, x00, velocity_xy_interp);
	  		x00 = boundary_cond(x00,n_pass_WE[j1]);
      		}
      		xf[j1*2] = x00[0] + 360*n_pass_WE[j1];
     		xf[j1*2+1] = x00[1];	
    	}
	 	// I compute backward Ftle
	 	t_ftle=(double)t0;
		ftle_back[(int)t0]= ftle_func(x0,xf,t_ftle,eigenvalue);	 	
	 	time_dir=1;	 
	 
	 
	 
	 }
	 
 	 // I compute betweenness
    for(t0=0; t0<=tau; t0++)
    {
    betw[i]=betw[i]+exp(ftle_back[(int)t0]*double(t0))*exp(ftle_forw[(int)t0]*(tau-double(t0)))/tau;
    }

        
    
    }

  ////////////////////////////////// OUTPUT WRITING ///////////////////////
  sdelta_betw = "_delta0betw";
  sint_step = "_intstep";
    
  sfile_betw=
    output_dir + 
    prod_name + 
    sdelta_betw + leadingzeros(6,delta0_ftle) + 
    sstart_date +
    stau +
    sint_step + leadingzeros(3,int_step) +
    postfixbetw;

 ofstream file_betw(sfile_betw.c_str());
 
   for(i=0; i<num_itracers; i++)
    {
    file_betw << betw[i]<< endl;
    }
 
  file_betw.close();  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ////////////////////////////////// OUTPUT WRITING ///////////////////////

/*  string sfile_matrix;
  
  // Construct name file
  sfile_matrix=
    output_dir + 
    prod_name + 
//AB    slayer_index + 
    snode_size + leadingzeros(4,node_size) + 
    sparticle_latspacing + leadingzeros(6,particle_latspacing) + 
    sstart_date + leadingzeros_int(2,start_date.tm_mday) + leadingzeros_int(2,start_date.tm_mon+1) + leadingzeros_int(4,start_date.tm_year) +
    stau +
    sint_step + leadingzeros(3,int_step) +
    postfixmatrix;*/

 /* ofstream file_matrix(sfile_matrix.c_str());
  int sum_tracers = 0;
 
  // Write only non-null elements of the matrix in the output file
  for(int i=0; i<num_nodes; i++)
    {
      for(int j=0; j<num_nodes; j++)
	{
	  if(trans_matrix[i][j]>0)
	    {
	      file_matrix << i <<" "<< j <<" "<< trans_matrix[i][j] << endl;
	      sum_tracers+=trans_matrix[i][j];
	    }
	}
    }
  
  // Print to screen sum of all elements of the matrix
  if(verbose == 1) 
    {
      cout << " Total sum of matrix elements (should equal total number of tracers) = " << sum_tracers <<endl;
    }

  // Close file
  file_matrix.close();*/

  //DEALLOCATE MEMORY
  delete[] itracer_lat;
  delete[] itracer_lon;
//  delete ftle;
  // deallocate memory of velocity field
  deallocate_vfield(n_velfield2load); //AB

/*  delete[] ftracer_lat;
  delete[] ftracer_lon;
  delete[] land_ratio;
  delete[] delta_node_lon;
  delete[]  max_node_lon;
  delete[]  max_jindex_node;
  for (i=0; i<max_iindex_node; i++)
    {
      delete[] node_index[i];
    }
  delete[] node_index;
  for (int i=0; i<num_nodes; i++)
    {
      delete[] trans_matrix[i];
    }
  delete[] trans_matrix;*/
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////// functions used in the main code (declared at the beginning of the main code)
/////////////////////////////////////////////////////////////////////////////////////////////////

// This function prints a little command line manual if there are invalid parameters 
void print_usage(const string me) 
{
  cout << "Usage: " << me << " <file_parameters>" << endl << endl;
}

// This function counts the number of lines of a file
int countlines(string namefile)
{
  ifstream ifile(namefile.c_str());
  int numlines;
  string line;
  
  if (!ifile.is_open()) 
    {
      string me = "countlines()";
      cout << me <<": Skipping unreadable file \"" << namefile.c_str() << "\" "<<endl;
      return -1;
    }
  for ( numlines = 0; getline(ifile, line); ++numlines)
    ;

  ifile.close();
  return numlines;
}

/* 
 * it only works for numbers < 10.0 
 */
string leadingzeros(unsigned int ndigits, double number)
{
  ostringstream stream;// it needs to include  <sstream>
  double factor=pow(10,ndigits-1); // it needs to include <cmath>
 
  stream << setprecision(ndigits) << setw(ndigits) << setfill('0') << floor(number*factor);
  
  return stream.str();
}
/* function leadingzeros_int:
 * it only works for integer numbers
 */
string leadingzeros_int(unsigned int ndigits, int number)
{
  ostringstream stream;// it needs to include  <sstream>
 
  stream << setw(ndigits) << setfill('0') << number;
  
  return stream.str();
}

// FUNCTION to check specifically boundary condition (to maintain all particles in the domain)
vector<double> boundary_cond(vector<double> grib_tr, int &pass_WE_YN) 								
{
 if (periodic_conditionWE_YN==0)
 {
  if (((start_tracer_long < end_tracer_long) &&
       (grib_tr[0] < start_tracer_long || grib_tr[0] > end_tracer_long)) // if first lon of domain (velocity field) is lower than last lon of domain [e.g. from -150 to 150]
      || 
      ((start_tracer_long >= end_tracer_long)
       && (grib_tr[0] < start_tracer_long && grib_tr[0] > end_tracer_long))) // if first lon of domain (velocity field) is greater than last lon of domain [e.g. from 150 to -150]
    {
      if (pow((start_tracer_long - grib_tr[0]),2) < pow((end_tracer_long - grib_tr[0]),2))
	grib_tr[0] = start_tracer_long;
      else
	grib_tr[0] = end_tracer_long;
    }
  }
  else
  {
  if (grib_tr[0] < -180)
  	{grib_tr[0]=grib_tr[0]+360;
  	pass_WE_YN=pass_WE_YN-1;}
  if (grib_tr[0] >= 180)
  	{grib_tr[0]=grib_tr[0]-360;
  	pass_WE_YN=pass_WE_YN+1;}
  }
  if (grib_tr[1] < start_tracer_lat)
    grib_tr[1] = start_tracer_lat;
  
  if (grib_tr[1] > end_tracer_lat)
    grib_tr[1] = end_tracer_lat;
  
  return grib_tr;
}

// FUNCTION interpolating (in space and time) the velocity field for the point whose coordinates are (cord_xy_vec[0]=lon, cord_xy_vec[1]=lat, cord_t=t) 
// To obtain velocity field at the exact positions of particle along its trajectory
vector<double> velocity_xy_interp (double cord_t, vector<double> cord_xy_vec)			
{
  int i, j, k, jl, jr;//AB l and r suffix indicate, resp., right and left
  
  double ii, jj;
  double t, u, dfrac,cord_t2;//AB
  
  vector<double> v_xy(2);
  
  // transformation from coordinates to indices
  jj = (cord_xy_vec[0] - v_field_lon[0]) * degree_resolution_inv;//AB
  ii = (cord_xy_vec[1] - v_field_lat[0]) * degree_resolution_inv;//AB		
  
  i = (int) ii;
 
  if (periodic_conditionWE_YN==0) //AB If periodic conditions are not present, the  script is exactly the same as before (I just changed the order of some lines)
  {
  j = (int) jj;
  t = jj - ((double) j);
  
  jl=j;
  jr=j+1;    
  }
  else //AB If periodic conditions are present:
  {
  
  if ((jj>=0) && (cord_xy_vec[0] < v_field_lon[nlon-1])) //AB If the particle longitude is comprised between the first and the last lon delimiting the velocity field (excluding the exact value of the right limit), I do as before
  {j = (int) jj;
  t = jj - ((double) j);
  jl=j;
  jr=j+1;}
  else //AB If the particle longitude is not comprised between the first and the last lon delimiting the velocity field, I apply the periodic conditions
  {jl=nlon-1;
  jr=0;
  if (jj<0) //AB If the particle longitude is lower than the left limit of the vel field (eg -180.1)
  {t = (360 + cord_xy_vec[0] - v_field_lon[nlon-1]) * degree_resolution_inv;}
  else  //AB If the particle longitude is higher or equal to the right limit of the velocity field (eg 179.875 for the Glob_dt product)
  {t = (cord_xy_vec[0] - v_field_lon[nlon-1]) * degree_resolution_inv;}}
  
  }
  //  cout << "i " << i << endl;
  //  cout << "j " << j << endl;
  
  u = ii - ((double) i);  
  
  cord_t2 = cord_t/vel_field_timestep;//AB
  k = (int) cord_t2;//AB
  dfrac = cord_t2 - ((double) k);//AB
 
  
  v_xy[0] = ((1-u) * (1-t)) * ((1-dfrac) * v_field_U[k][i][jl] + dfrac * v_field_U[k+1][i][jl])
    + ((1-u) * t) * ((1-dfrac)*(v_field_U[k][i][jr]) + (dfrac)*(v_field_U[k+1][i][jr]))
    + (u * t) * ((1-dfrac)*(v_field_U[k][i+1][jr]) + (dfrac)*(v_field_U[k+1][i+1][jr]))
    + (u * (1-t)) * ((1-dfrac)*(v_field_U[k][i+1][jl]) + (dfrac)*(v_field_U[k+1][i+1][jl]));
  
  v_xy[1] = ((1-u) * (1-t)) * ((1-dfrac)*(v_field_V[k][i][jl]) + (dfrac)*(v_field_V[k+1][i][jl]))
    + ((1-u) * t) * ((1-dfrac)*(v_field_V[k][i][jr]) + (dfrac)*(v_field_V[k+1][i][jr]))
    + (u * t) * ((1-dfrac)*(v_field_V[k][i+1][jr]) + (dfrac)*(v_field_V[k+1][i+1][jr]))
    + (u * (1-t)) * ((1-dfrac)*(v_field_V[k][i+1][jl]) + (dfrac)*(v_field_V[k+1][i+1][jl]));

  
  // if (dfrac == 0){
  //   cout << "*******************************************************  " << t << " " << u <<  "\n";
  //   cout << "velocity(x,y) -> ( " << v_xy[0] << " , "<< v_xy[1] << " ) \n";
  //   cout << "field(x,y) -> ( " << vel_field_U[k][i][j] << " , "<<  vel_field_V[k][i][j] << " ) \n";
  //   cout << "field(x,y) -> ( " << vel_field_U[k][i+1][j] << " , "<<  vel_field_V[k][i+1][j] << " ) \n";
  //   cout << "field(x,y) -> ( " << vel_field_U[k][i][jp1] << " , "<<  vel_field_V[k][i][jp1] << " ) \n";
  //   cout << "field(x,y) -> ( " << vel_field_U[k][i+1][jp1] << " , "<<  vel_field_V[k][i+1][jp1] << " ) \n";
  // }
  
  // Units of the interpolated velocities are in degrees per day (inputs of the RK4)
  v_xy[0] = ((double) time_dir)*(v_xy[0] / (Radius * cos(cord_xy_vec[1] * deg2rad))) * (1. / deg2rad) * day2sec;
  v_xy[1] = ((double) time_dir)*(v_xy[1] / Radius) * (1. / deg2rad) * day2sec;
  
  return v_xy;
}

// FUNCTION RK4
vector<double> rungeKutta4(double start_cord_t,
			   vector<double> start_cord_xy_vec,
			   vector<double> velocity_xy_interp (double cord_t,
							      vector<double> cord_xy_vec))
{	
  //	total number of elements in the vector
  int n = start_cord_xy_vec.size();
	
  //	first step
  vector<double> k1;
  k1 = velocity_xy_interp(start_cord_t, start_cord_xy_vec);

  // cout << "start cord t " << start_cord_t << endl;
  // cout << "start cord vec " << start_cord_xy_vec[0] << " " << start_cord_xy_vec[1] << endl;
  // cout << "k1 " << k1[1] << endl;


  for (int i=0; i<n; ++i) {
    k1.at(i) *= int_step;        
  }

  //	second step
  vector<double> k2(start_cord_xy_vec);
  for (int i=0; i<n; ++i) {
    k2.at(i) += k1.at(i) / 2.;
  }


  k2 = velocity_xy_interp(start_cord_t + int_step / 2., k2);

  for (int i=0; i<n; ++i) {
    k2.at(i) *= int_step;
  }

  //	third step
  vector<double> k3(start_cord_xy_vec);
  for (int i=0; i<n; ++i) {
    k3.at(i) += k2.at(i) / 2.;
  }
  k3=velocity_xy_interp(start_cord_t + int_step / 2., k3);
  for (int i=0; i<n; ++i) {
    k3.at(i) *= int_step;
  }

  //	fourth step
  vector<double> k4(start_cord_xy_vec);
  for (int i=0; i<n; ++i) {
    k4.at(i) += k3.at(i);
  }
  k4=velocity_xy_interp(start_cord_t + int_step, k4);
  for (int i=0; i<n; ++i) {
    k4.at(i) *= int_step;
  }

  //	sum the weighted steps into yf and return the final y values
  vector<double> yf(start_cord_xy_vec);
  for (int i=0; i<n; ++i) {
    yf.at(i) = yf.at(i)
      + k1.at(i) / 6.
      + k2.at(i) / 3.
      + k3.at(i) / 3.
      + k4.at(i) / 6.;
  }	
  return yf;
}

// eigenvalues computation function
double eigenvalue(vector<double> Delta)									
	{

	double largest;
	vector<double> eigen(2);
	double trace;
	double deter;

	trace = Delta[0] + Delta[3];
	deter = Delta[0] * Delta[3] - Delta[1] * Delta[2];

	eigen[0] = (0.5 * (trace - sqrt(trace * trace - 4 * deter)));
	eigen[1] = (0.5 * (trace + sqrt(trace * trace - 4 * deter)));

	if(sqrt(eigen[0] * eigen[0]) > sqrt(eigen[1] * eigen[1]))
		largest = eigen[0];
	else
		largest = eigen[1];


	return largest;

}


// ftle computation function
double ftle_func(vector<double> x0,vector<double> xf,double t_ftle,double eigenvalue(vector<double> Delta))									
	{
     vector<double> jacob_m(4),transpose_jacob_m(4),Delta2(4);
     double ftle;
	if (t_ftle==0)
	{
	ftle=0;return ftle; 
	}
    jacob_m[0]=(xf[2]-xf[6])/(x0[2]-x0[6]);
    jacob_m[1]=(xf[4]-xf[8])/(x0[5]-x0[9]);
    jacob_m[2]=(xf[3]-xf[7])/(x0[2]-x0[6]);
    jacob_m[3]=(xf[5]-xf[9])/(x0[5]-x0[9]);
	
	transpose_jacob_m[0]=jacob_m[0];
    transpose_jacob_m[1]=jacob_m[2];
    transpose_jacob_m[2]=jacob_m[1];
    transpose_jacob_m[3]=jacob_m[3];

    Delta2[0]=(transpose_jacob_m[0]) * (jacob_m[0]) + (transpose_jacob_m[1]) * (jacob_m[2]);
    Delta2[1]=(transpose_jacob_m[0]) * (jacob_m[1]) + (transpose_jacob_m[1]) * (jacob_m[3]);
    Delta2[2]=(transpose_jacob_m[2]) * (jacob_m[0]) + (transpose_jacob_m[3]) * (jacob_m[2]);
    Delta2[3]=(transpose_jacob_m[2]) * (jacob_m[1]) + (transpose_jacob_m[3]) * (jacob_m[3]);
    
    ftle = (1/t_ftle)*log(sqrt(eigenvalue(Delta2)));
return ftle;
}
