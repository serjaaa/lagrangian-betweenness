//To update this file and the others in the same folder on niko:
//rsync -rv -e 'ssh -p 8022' ~/cpp/code_Betweenness/*.* abaudena@localhost:~/cpp/code_Betweenness/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>  
#include <math.h>

using namespace std;
#include "parameters.h"
#define PI 3.14159265

enum { /* 0 */ LAND_SEA_NORM,//AB
       /* 1 */ VEL_PRODUCT,
       /* 2 */ VERBOSE,
       /* 3 */ NODE_SIZE,
       /* 4 */ PARTICLE_LATSPACING,
       /* 5 */ NETWORK_LL_LON,
       /* 6 */ NETWORK_TR_LON,
       /* 7 */ NETWORK_LL_LAT,
       /* 8 */ NETWORK_TR_LAT,
       /* 9 */ LAYER_INDEX,
       /* 10 */ START_DATE,
       /* 11 */ TAU,
       /* 12 */ DELTA0_FTLE,
       /* 13 */ DELTA_FSLE_FIN,//AB_s
       /* 14 */ TRJ_TIMESTEP,
       /* 15 */ RAND_IDX,
       /* 16 */ OPZ_CLIM,//AB_e
       NPARAMETERS
};

//ABint domain_network;
int land_sea_norm;
int vel_product;
int verbose;
double node_size;
double particle_latspacing;
double network_ll_lat;//AB_s
double network_ll_lon;
double network_tr_lat;
double network_tr_lon;
double delta0_ftle;
double delta_fsle_fin;
int trj_timestep;
int rand_idx;
int opz_clim;//AB_e
int layer_index;
struct tm start_date= {0};
/* The structure tm contains nine members of type int (in any order), which are:
 * tm_sec	int	seconds after the minute 0-61
 * tm_min	int	minutes after the hour	 0-59
 * tm_hour	int	hours since midnight	 0-23
 * tm_mday	int	day of the month	 1-31
 * tm_mon	int	months since January	 0-11
 * tm_year	int	years since 1900	
 * tm_wday	int	days since Sunday	 0-6
 * tm_yday	int	days since January 1	 0-365
 * tm_isdst	int	Daylight Saving Time flag 
 */
double tau;

//AB_s
string ncdir,u_name_nc,v_name_nc;
string lon_name_nc,lat_name_nc;
string prod_name,prod_prefix;
int hour_resYN;
double vel_field_timestep;
int time_dir;
double int_step;
double long_corr;//longitudinal correction for particles position
//Total domain of the velocity field
double lon_vfield_min;
double lon_vfield_max;
double lat_vfield_min;
double lat_vfield_max;
int periodic_conditionWE_YN; //It indicates if the velocity field is periodic (1) or not (0) along the west-east direction (i.e., longitude). So basically it is 1 for global products, 0 for regional ones.
/* vfield domain used for the advection (can be a portion of the total one since it may not be necessary to use it all)*/
double degree_resolution,degree_resolution_inv,degree_resolution2_in_meters_inv;
double deg2m=111226.2636822; // converts from degrees to meters
double vfield_ll_lat;
double vfield_ll_lon;
double vfield_tr_lat;
double vfield_tr_lon;


//AB_e

/* DOMAIN network: this must be in grid_construction.cpp*/
/*AB_s double network_ll_lat;
double network_ll_lon;
double network_tr_lat;
double network_tr_lon;*///AB_e
double buffer_particles; //AB
double buffer_advection;//AB
double latmax;//AB

void listofparameters(string pname[], string ptype[])
{
  int parameter;

  cout << "In the file parameters there must be these parameters ";
  cout << "(it is not necessary the same order, but it is recommended):"<< endl;
  cout << "Number Name Type"<< endl;
  for(parameter=0; parameter<NPARAMETERS; parameter++)
    cout << parameter <<"->"<<pname[parameter]<<" "<<ptype[parameter]<<endl;
}

int readparams(string nfileparameters)
{
  string me="readparams()";
  ifstream fparameters(nfileparameters.c_str());

  int delimiter,end;
  string name, value;
  string line;


  int *pflag;

  int parameter;
  string pname[NPARAMETERS];
  string ptype[NPARAMETERS];
  
  if (!fparameters.is_open())
    {
      cout << me <<": Skipping unreadable file \"" << nfileparameters.c_str() << "\" "<<endl; 
      return 1;
    }

  pflag = (int*) calloc (NPARAMETERS,sizeof(int));

//AB  pname[DOMAIN_NETWORK]="domain_network";           ptype[DOMAIN_NETWORK]="int";
  pname[LAND_SEA_NORM]="land_sea_norm";             ptype[LAND_SEA_NORM]="int"; 
  pname[VEL_PRODUCT]="vel_product";           		ptype[VEL_PRODUCT]="int"; 
  pname[VERBOSE]="verbose";                         ptype[VERBOSE]="int"; 
  pname[NODE_SIZE]="node_size";                     ptype[NODE_SIZE]="double"; 
  pname[PARTICLE_LATSPACING]="particle_latspacing"; ptype[PARTICLE_LATSPACING]="double"; 
  pname[NETWORK_LL_LON]="network_ll_lon"; 			ptype[NETWORK_LL_LON]="double"; 
  pname[NETWORK_TR_LON]="network_tr_lon"; 			ptype[NETWORK_TR_LON]="double"; 
  pname[NETWORK_LL_LAT]="network_ll_lat"; 			ptype[NETWORK_LL_LON]="double"; 
  pname[NETWORK_TR_LAT]="network_tr_lat"; 			ptype[NETWORK_TR_LAT]="double"; 
  pname[LAYER_INDEX]="layer_index";                 ptype[LAYER_INDEX]="int"; 
  pname[START_DATE]="start_date";                   ptype[START_DATE]="double"; 
  pname[TAU]="tau";                                 ptype[TAU]="double"; 
  pname[DELTA0_FTLE]="delta0_ftle";           		ptype[DELTA0_FTLE]="double"; //AB_s
  pname[DELTA_FSLE_FIN]="delta_fsle_fin";           ptype[DELTA_FSLE_FIN]="double"; 
  pname[TRJ_TIMESTEP]="trj_timestep";               ptype[TRJ_TIMESTEP]="int"; 
  pname[RAND_IDX]="rand_idx";             		    ptype[RAND_IDX]="int"; 
  pname[OPZ_CLIM]="opz_clim";             		    ptype[OPZ_CLIM]="int"; //AB_e

  while(!fparameters.eof())
    {
      getline(fparameters,line);
      if (line[0] == '#') 
	continue;  /* ignore comment line which starts with # */

      delimiter = line.find("=");
      end = line.length();
      if (end == 0) 
	continue; /* ignore blank line */
      
      value = line.substr(delimiter+1, end);
      name = line.substr(0, delimiter);

      for(parameter=0; parameter<NPARAMETERS; parameter++)
	{
	  if (name.compare(pname[parameter]) == 0)
	    {
	      pflag[parameter]++;
	      if(pflag[parameter]>1)
		{
		  cout  << me << ": Parameter "<< name << " repeated"<<endl;
		  listofparameters(pname, ptype);
		  return 1;
		}
	       if ((delimiter+1) == end) 
		 {
		   cout << me << ": Parameter "<< name << " has no value"<<endl;
		   listofparameters(pname, ptype);
		   return 1;
		 }
	       switch (parameter) 
		{
/*AB_s		case DOMAIN_NETWORK:{
		  domain_network = atoi(value.c_str());}*///AB_e
		  break;
		case LAND_SEA_NORM:
		  land_sea_norm = atoi(value.c_str());
		  break;
		case VEL_PRODUCT:
                  vel_product = atoi(value.c_str());
		  break;
		case VERBOSE:
		  verbose = atoi(value.c_str());
		  break;
		case NODE_SIZE:
		  node_size = atof(value.c_str());
		  break;
		case PARTICLE_LATSPACING:
		  particle_latspacing = atof(value.c_str());
		  break;
		case NETWORK_LL_LON:
		  network_ll_lon = atof(value.c_str());
		  break;
		case NETWORK_TR_LON:
		  network_tr_lon = atof(value.c_str());
		  break;
		case NETWORK_LL_LAT:
		  network_ll_lat = atof(value.c_str());
		  break;
		case NETWORK_TR_LAT:
		  network_tr_lat = atof(value.c_str());
		  break;
		case LAYER_INDEX:
		  layer_index = atoi(value.c_str());
		  break;
		case DELTA0_FTLE://AB_s
		  delta0_ftle = atof(value.c_str());
		  break;
		case DELTA_FSLE_FIN:
		  delta_fsle_fin = atof(value.c_str());
		  break;
		case TRJ_TIMESTEP:
		  trj_timestep = atoi(value.c_str());
		break;
		case RAND_IDX:
		  rand_idx = atoi(value.c_str());
		break;
		case OPZ_CLIM:
		  opz_clim = atoi(value.c_str());
		break;//AB_e
		
		case START_DATE:
		  if(sscanf(value.c_str(),"%d-%d-%d;%d",&start_date.tm_mday,&start_date.tm_mon,&start_date.tm_year,&start_date.tm_hour)!=4)//AB
		    {
		      cout  << me << ": Date format in start_date is incorrect" << endl;
		      listofparameters(pname, ptype);
		      return 1;
		    }
		  start_date.tm_mon = start_date.tm_mon - 1; 
//AB	  start_date.tm_hour= 0;
		  start_date.tm_min = 0; 
		  start_date.tm_sec = 0;
		  break;
		case TAU:
                  tau = atof(value.c_str());
                  if (tau>=0)
                  {time_dir=1;}
                  else
                  {time_dir=-1;tau=-tau;}
		  break;
		default:
		  cout  << me << ": Unknown parameter "<< name <<endl ;
		  listofparameters(pname, ptype);
		  return 1;
		}
	       break;
	    }
	}
      if(parameter==NPARAMETERS)
	{
	  cout  << me << ": Unknown parameter "<< name <<endl ;
	  listofparameters(pname, ptype);
	  return 1;
	}
    }

  for(parameter=0; parameter<NPARAMETERS; parameter++)
	{
	  if(pflag[parameter]==0)
	    {
	      cout  << me << ": parameter "<< pname[parameter] << " is not defined"<<endl;
	      listofparameters(pname, ptype);
	      return 1;
	    }
	}

  fparameters.close();

if (vel_product==1)
    {
      prod_name="Med_combined";
      prod_prefix="Med_combined";
      degree_resolution = 0.04166666667; // 1/24 degrees
      ncdir = "/home/abaudena/matlab/Velocity_field/data_vel/Med_combined/";
      lon_name_nc="lon";lat_name_nc="lat";
      u_name_nc="uc";v_name_nc="vc";//how the vel field is called in the netcdf file
      int_step=1./48.; //time step of the integration (RK4) in days
      vel_field_timestep=1./24.; // vel field time resolution
      lon_vfield_min=-6.00; //total domain of the velocity field in the netcdf file
      lon_vfield_max=36.2917;
      lat_vfield_min=30.1875;
      lat_vfield_max=45.9792;
      periodic_conditionWE_YN=0; //It indicates if the velocity field is periodic (1) or not (0) along the west-east direction (i.e., longitude). So basically it is 1 for global products, 0 for regional ones.
      
    }
else if (vel_product==2)
    {
      prod_name="Glob_dt";
      prod_prefix="Glob_dt";
      degree_resolution = 0.25; // 1/4 degrees
      ncdir = "/home/abaudena/matlab/Velocity_field/data_vel/Glob_dt/";
      lon_name_nc="lon";lat_name_nc="lat";
      u_name_nc="ugos";v_name_nc="vgos";//how the vel field is called in the netcdf file
      int_step=1./8.; //time step of the integration (RK4) in days
      vel_field_timestep=1; // vel field time resolution
      lon_vfield_min=-179.8750; //total domain of the velocity field in the netcdf file
      lon_vfield_max=179.8750;
      lat_vfield_min=-89.8750;
      lat_vfield_max=89.8750;
      start_date.tm_hour= 0;
      periodic_conditionWE_YN=1; //It indicates if the velocity field is periodic (1) or not (0) along the west-east direction (i.e., longitude). So basically it is 1 for global products, 0 for regional ones.
    }
else if (vel_product==3)
    {
      prod_name="Glob_Ekman_dt";
      prod_prefix="ekman_dt";
      degree_resolution = 0.25; // 1/4 degrees
      ncdir = "/home/abaudena/cpp/code_Betweenness/Velocity_field/Glob_Ekman_dt/";
      lon_name_nc="lon";lat_name_nc="lat";
      u_name_nc="uo";v_name_nc="vo";//how the vel field is called in the netcdf file
      int_step=1./8.; //time step of the integration (RK4) in days
      vel_field_timestep=1; // vel field time resolution
      lon_vfield_min=-179.8750; //total domain of the velocity field in the netcdf file
      lon_vfield_max=179.8750;
      lat_vfield_min=-89.8750;
      lat_vfield_max=89.8750;
      start_date.tm_hour= 0;
      periodic_conditionWE_YN=1; //It indicates if the velocity field is periodic (1) or not (0) along the west-east direction (i.e., longitude). So basically it is 1 for global products, 0 for regional ones.
    }
else if (vel_product==4)
    {
      prod_name="GloryS12";
      prod_prefix="GloryS12";
      degree_resolution = 1./12.; // 1/12 degrees
      ncdir = "/home/abaudena/matlab/Velocity_field/data_vel/GloryS12/";
      lon_name_nc="lon";lat_name_nc="lat";
      u_name_nc="uo";v_name_nc="vo";//how the vel field is called in the netcdf file
      int_step=1./8.; //time step of the integration (RK4) in days
      vel_field_timestep=1; // vel field time resolution
      lon_vfield_min=-180; //total domain of the velocity field in the netcdf file
      lon_vfield_max=179.9167;
      lat_vfield_min=-80;
      lat_vfield_max=90;
      start_date.tm_hour= 12;
      periodic_conditionWE_YN=1; //It indicates if the velocity field is periodic (1) or not (0) along the west-east direction (i.e., longitude). So basically it is 1 for global products, 0 for regional ones.
    }
else
       {
	 cout << "Incorrect value of velocity product variable"<<endl ;
	 listofparameters(pname, ptype);
	 return 1;
       }

  degree_resolution_inv=1./degree_resolution;
  degree_resolution2_in_meters_inv=1./(2*degree_resolution*deg2m);
	// Max size of network grid should be smaller than velocity field (if 10 -> 0.6 deg margin / if 5 -> 0.3 deg margin)
	//AB Partilces can not be put at the boundary of the vel field domain, otherwise, it they exit during the RK4
	// the code crashes. So I set a buffer distance which is computed in this way. The max velocity is around 1.5m/s. 
	//The int_step is xx hours. So a particle can travel around 5.0 km/h maximum. Along latitude, buffer is 5.0 km/h*int_step*24/111
	//Along longitude, I have a latitudinal correction, thus I divide per the cosinus of the maximal latitude
	//The same approach applies below for the buffer for the advection
	buffer_particles=3.6*int_step*24/111;
	// I identify the max lat in the domain, to correct the buffer	
	if (lat_vfield_min*lat_vfield_min>lat_vfield_max*lat_vfield_max)
	{latmax=lat_vfield_min;}
	else
	{latmax=lat_vfield_max;}
	if (latmax*latmax>60*60) {latmax=60;}
	long_corr=1./cos ( latmax * PI / 180.0 );
//	cout <<"long_corr "<< long_corr<<endl;
	// I compute the buffer for the velocity field (which is based on advection)
	buffer_advection=6*24*tau/111;	
	// I compute the limits of the velocity field I will need to load
  	vfield_ll_lon = network_ll_lon-buffer_advection*long_corr; //AB I change v_field_lon[0] to lon_vfield_min+buffer_particles.. and so on, so that particles do not exit the domain (and crash the code)
 	vfield_tr_lon = network_tr_lon+buffer_advection*long_corr; 
 	vfield_ll_lat = network_ll_lat-buffer_advection; 
  	vfield_tr_lat = network_tr_lat+buffer_advection;  
  	
  	// Here I correct them in case they are outside the domain
  	if (vfield_ll_lat<lat_vfield_min||vfield_ll_lat>lat_vfield_max) {vfield_ll_lat=lat_vfield_min;}
	if (vfield_tr_lat>lat_vfield_max||vfield_tr_lat<lat_vfield_min) {vfield_tr_lat=lat_vfield_max;}
	if (periodic_conditionWE_YN==0) //AB If periodic conditions are not present, the script for the  longitude limits remains as before
	{
	if (vfield_ll_lon<lon_vfield_min||vfield_ll_lon>lon_vfield_max) {vfield_ll_lon=lon_vfield_min;}
	if (vfield_tr_lon>lon_vfield_max||vfield_tr_lon<lon_vfield_min) {vfield_tr_lon=lon_vfield_max;}
	}
  	else
  	{//AB If periodic conditions are present, and if the particle can possibly cross the 180E line, I load the whole velocity field to assure boundary periodic conditions
  	if (vfield_ll_lon<lon_vfield_min||vfield_tr_lon>lon_vfield_max) {vfield_ll_lon=lon_vfield_min;vfield_tr_lon=lon_vfield_max;}
  	}

  /* VERBOSE *///AB:I moved this if statement at the end
  if(verbose == 1)
    {
      cout << "PARAMETERS FROM FILE: "<< nfileparameters <<endl; 
      cout << " land_sea_norm = " <<land_sea_norm << endl;
      cout << " velocity product = "<<vel_product<< endl ;
      cout << " product name = "<<prod_name<< endl ;
      cout << " verbose = "<< verbose<< endl;
      cout << " node_size = "<< node_size<< endl;
      cout << " particle_latspacing = "<<particle_latspacing<<endl ;
      cout << " Domain low left lon = "<<network_ll_lon<<endl ;
      cout << " Domain top right lon = "<<network_tr_lon<<endl ;
      cout << " Domain low left lat = "<<network_ll_lat<<endl ;
      cout << " Domain top right lat = "<<network_tr_lat<<endl ;
      cout << " start_date = "<< start_date.tm_mday<<"-"<<start_date.tm_mon+1<<"-"<<start_date.tm_year<<endl ;
      cout << " tau = "<<tau<< endl ;
      cout << " int_step = "<<int_step<< endl ;
      cout << " delta0_ftle = "<<delta0_ftle<< endl ;
      cout << " delta_fsle_fin = "<<delta_fsle_fin<< endl ;
      cout << " trj_timestep = "<<trj_timestep<< endl ;
      cout << " rand_idx = "<<rand_idx<< endl ;
      cout << " opz_clim = "<<opz_clim<< endl ;
    }

  /* DOMAIN */

//AB  enum { /* 1 */ WHOLE_MEDSEA=1,
//AB	 /* 2 */ EASTERN_MEDSEA,
//AB	 /* 3 */ WESTERN_MEDSEA
//AB  };

/*AB_s    switch(domain_network) {
    case WHOLE_MEDSEA: //Whole MedSea
      {
//	cout <<"buffer"<<buffer_particles<<endl;
	network_ll_lat= 30.1875+buffer_particles;
	network_ll_lon= -6.0+buffer_particles*1.41421356;
	network_tr_lat= 45.9792-buffer_particles;
	network_tr_lon= 36.2917-buffer_particles*1.41421356;
    vfield_ll_lat=30.1875;
	vfield_ll_lon=-6.00;
	vfield_tr_lat=45.9792;
	vfield_tr_lon=36.2917;
      }
      break;
    case EASTERN_MEDSEA: //Eastern MedSea
      {
	network_ll_lat=31.0;
	network_ll_lon=10.0;
	network_tr_lat=44.0;
	network_tr_lon=36.0;
	vfield_ll_lat=30.1875;
	vfield_ll_lon=-6.0;
	vfield_tr_lat=45.9375;
	vfield_tr_lon=36.25;
      }
      break;
    case WESTERN_MEDSEA: //Western MedSea
      {
	network_ll_lat=31.0;
	network_ll_lon=-5.0;
	network_tr_lat=44.0;
	network_tr_lon=17.0;
	vfield_ll_lat=30.1875;
	vfield_ll_lon=-6.0;
	vfield_tr_lat=45.9375;
	vfield_tr_lon=36.25;
      }
      break;
     default: 
       {
	 cout << "Incorrect value of domain_network variable"<<endl ;
	 listofparameters(pname, ptype);
	 return 1;
       }
    }*///AB_s

    //deallocate memory:

    free (pflag);
    return 0;
}
