//To update this file and the others in the same folder on niko:
//rsync -rv -e 'ssh -p 8022' ~/cpp/code_Betweenness/*.* abaudena@localhost:~/cpp/code_Betweenness/

#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <ctime>  
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;
#include "parameters.h"

//Parameters 
static const int seconds_day = 86400;

// output file
char name_out[256];
ofstream  fout;

// Return this code to the OS in case of failure.
static const int NC_ERR = 2;

//  Velocity field
double* v_field_lat; 
double* v_field_lon;
double* v_field_depth;
double*** v_field_U;
double*** v_field_V;
int nlat;
int nlon;
int ndepth;
//AB_s
// Tracer maximum possible domain
  double start_tracer_long;
  double start_tracer_lat;
  double end_tracer_long; 
  double end_tracer_lat; 
//AB_e

int read_vfield(struct tm start_date, int n_velfield2load, int time_dir)//AB, int layer_index)
{
  struct tm *date = {0};
  time_t start_time, time; //time_t: value representing the number of seconds elapsed since 00:00 hours, Jan 1, 1970 UTC
  
//AB  char const ncdir[] = "/home/abaudena/matlab/Velocity_field/data_vel/Med_combined/";
  char ncfile[256];
  
  int i,j,t;// loop index

  NcError err(NcError::verbose_nonfatal);
  NcDim *ncdim_lat, *ncdim_lon, *ncdim_depth;
  NcVar *latvar, *lonvar, *depthvar;
  NcVar *Uvar, *Vvar;
  NcAtt *att;

//AB  double longitude_min, longitude_max;
  int ilonmin, ilonmax;
//AB  double latitude_min, latitude_max;
  int ilatmin, ilatmax;
 //AB double z_min, z_max;

  // Reads the grid lon, lat and depth from the first netcdf file

	if (vel_field_timestep<1)
	{
  	sprintf(ncfile, "%s%s_%04d%02d%02d_%02d.nc",ncdir.c_str(),prod_prefix.c_str(), start_date.tm_year,start_date.tm_mon+1,start_date.tm_mday,start_date.tm_hour);
	}
	else
	{
	sprintf(ncfile, "%s%04d/%s_%04d%02d%02d.nc",ncdir.c_str(),start_date.tm_year,prod_prefix.c_str(),start_date.tm_year, start_date.tm_mon+1, start_date.tm_mday);
	}
  cout << "reading first netcdf file " <<ncfile << endl;
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);
  
  // Check to see if the file was opened.
  if(!dataFile.is_valid())
    return NC_ERR;
  
  // Get dimensions of variables

  int dim_lat; 
  if (!(ncdim_lat = dataFile.get_dim("latitude")))
    return NC_ERR;
  dim_lat = ncdim_lat->size();


  int dim_lon; 
  if (!(ncdim_lon = dataFile.get_dim("longitude")))
    return NC_ERR;
  dim_lon = ncdim_lon->size();


  /*AB
  int dim_depth; 
  if (!(ncdim_depth = dataFile.get_dim("depth")))
    return NC_ERR;
  dim_depth = ncdim_depth->size();*/
 

  // Get pointers to the latitude and longitude variables.
  if (!(latvar = dataFile.get_var(lat_name_nc.c_str())))
    return NC_ERR;
  if (!(lonvar = dataFile.get_var(lon_name_nc.c_str())))
    return NC_ERR;
  /*AB
  if (!(depthvar = dataFile.get_var("depth")))
    return NC_ERR;*/


  // // Get max and min values of variables
  // if (!(att = latvar->get_att("valid_min")))
  //   return NC_ERR;
  // latitude_min = att->as_double(0);
  
  // delete att;  // Attributes and attribute values should be deleted by the caller
  //              // when no longer needed, to prevent memory 

  // if (!(att = latvar->get_att("valid_max")))
  //   return NC_ERR;
  // latitude_max = att->as_double(0);

  // delete att;  // Attributes and attribute values should be deleted by the caller
  //              // when no longer needed, to prevent memory leaks.

  // if (!(att = lonvar->get_att("valid_min")))
  //   return NC_ERR;
  // longitude_min = att->as_double(0);
  
  // delete att;  // Attributes and attribute values should be deleted by the caller
  //              // when no longer needed, to prevent memory leaks.

  // if (!(att = lonvar->get_att("valid_max")))
  //   return NC_ERR;
  // longitude_max = att->as_double(0);
  
  // delete att;   // Attributes and attribute values should be deleted by the caller
  //               // when no longer needed, to prevent memory leaks.

  // if (!(att = depthvar->get_att("valid_min")))
  //   return NC_ERR;
  // z_min = att->as_double(0);
  
  // delete att;   // Attributes and attribute values should be deleted by the caller
  //               // when no longer needed, to prevent memory leaks.

  // if (!(att = depthvar->get_att("valid_max")))
  //   return NC_ERR;
  // z_max = att->as_double(0);
  
  // delete att;   // Attributes and attribute values should be deleted by the caller
  //               // when no longer needed, to prevent memory leaks.
  

/*AB_s  longitude_min=-6.0;
  longitude_max=36.2917;

  latitude_min=30.1875;
  latitude_max=45.9792;*///AB_e

/*AB
  z_min=3.13;
  z_max=9.765777;*/


  ilatmin = ((vfield_ll_lat - lat_vfield_min)/degree_resolution); 
  ilatmax = ((vfield_tr_lat - lat_vfield_min)/degree_resolution);
  nlat = ilatmax-ilatmin+1;

  ilonmin = ((vfield_ll_lon - lon_vfield_min)/degree_resolution);
  ilonmax = ((vfield_tr_lon - lon_vfield_min)/degree_resolution);
  nlon = ilonmax-ilonmin+1;

  //AB ndepth = dim_depth;

  if(verbose==1)
    {
      cout << "PARAMETERS IN read_vfield" << endl;
      cout << " Low left vfield corner (lat,long)=(" << vfield_ll_lat << "," << vfield_ll_lon << ")"<<endl;
      cout << " Top right vfield corner (lat,long)=(" << vfield_tr_lat << "," << vfield_tr_lon << ")"<<endl;    
      cout << " ilatmax=" << ilatmax << endl;
      cout << " ilatmin=" << ilatmin << endl;
      cout << " nlat=" << nlat << endl;
      cout << " ilonmax=" << ilonmax << endl;
      cout << " ilonmin=" << ilonmin << endl;
      cout << " nlon=" << nlon << endl;
      cout << " longitude min:"<<lon_vfield_min<< endl; 
      cout << " longitude max:"<<lon_vfield_max<< endl; 
      cout << " latitude min:"<<lat_vfield_min<< endl; 
      cout << " latitude max:"<<lat_vfield_max<< endl; 
 /*AB
      cout << " z min:"<<z_min<< endl; 
      cout << " z max:"<<z_max<< endl; */
    }


  if (!latvar->set_cur(ilatmin))
    return NC_ERR;
  if (!lonvar->set_cur(ilonmin))
    return NC_ERR;

  v_field_lat = new double [nlat];
  v_field_lon = new double [nlon];
//AB   v_field_depth = new double [ndepth];

  // Get the lat/lon data from the file.
  if (!latvar->get(v_field_lat, nlat))
    return NC_ERR;
  if (!lonvar->get(v_field_lon, nlon))
    return NC_ERR;
/*AB
  if (!depthvar->get(v_field_depth, ndepth))
    return NC_ERR;*/
    
if (verbose==1)
{    
      cout << " Effective longitude min:"<<v_field_lon[0]<< endl; 
      cout << " Effective longitude max:"<<v_field_lon[nlon-1]<< endl; 
      cout << " Effective latitude min:"<<v_field_lat[0]<< endl; 
      cout << " Effective latitude max:"<<v_field_lat[nlat-1]<< endl; 
}
  dataFile.close();// close Netcdf file

/*AB
  if(verbose==1)
    {
      cout << " depth:"<<v_field_depth[layer_index]<< endl; 
    }
*/

  // reads velocities from files

  // Dynamically allocate velocity fields
  v_field_U = new double ** [n_velfield2load];
  v_field_V = new double ** [n_velfield2load];
  for (i = 0; i < n_velfield2load; i++)
    {
      v_field_U[i] = new double * [nlat];
      v_field_V[i] = new double * [nlat];
      for (j = 0; j < nlat; j++)
	{
	  v_field_U[i][j] = new double [nlon];
	  v_field_V[i][j] = new double [nlon];
	}
    }
 
  start_time = mktime(&start_date); // convert date to time in seconds
  for(t=0; t<n_velfield2load; t++)
    {
      time = start_time + time_dir*t*seconds_day*vel_field_timestep;
      date = gmtime(&time);

//AB_s
if (t==0)
{
/*       cout << " buffer_particles:"<<buffer_particles<< endl; 
      cout << " long_corr:"<<long_corr<< endl; */
	// Here I compute the limits for the boundary conditions, based on the limits of the velocity field I need to load
  	start_tracer_long = v_field_lon[0]+buffer_particles*long_corr; //AB I change v_field_lon[0] to lon_vfield_min+buffer_particles.. and so on, so that particles do not exit the domain (and crash the code)
 	end_tracer_long = v_field_lon[nlon-1]-buffer_particles*long_corr; 
 	start_tracer_lat = v_field_lat[0]+buffer_particles; 
  	end_tracer_lat = v_field_lat[nlat-1]-buffer_particles;
  	
  	//Here I correct the network limits, in case they are outside the domain, and I print a warning in case I have to correct them  
  	int count=0;
 if (periodic_conditionWE_YN==0) //AB If periodic conditions are not present, the  script is exactly the same as before 
  {
	if (network_ll_lon<start_tracer_long||network_ll_lon>end_tracer_long) {network_ll_lon=start_tracer_long;count+=1;}
	if (network_tr_lon>end_tracer_long||network_tr_lon<start_tracer_long) {network_tr_lon=end_tracer_long;count+=1;}
  }
	if (network_ll_lat<start_tracer_lat||network_ll_lat>end_tracer_lat) {network_ll_lat=start_tracer_lat;count+=1;}
	if (network_tr_lat>end_tracer_lat||network_tr_lat<start_tracer_lat) {network_tr_lat=end_tracer_lat;count+=1;}

  	if (count>0)
	{
	cout <<"Warning: network grid extension exceeds vel field domain"<<endl;
	cout <<"Therefore, network grid extension was reduced "<<endl;	
	}

/*	  cout << " start_tracer_long = "<<start_tracer_long<< endl ;
      cout << " end_tracer_long = "<<end_tracer_long<< endl ;
      cout << " start_tracer_lat = "<<start_tracer_lat<< endl ;
      cout << " end_tracer_lat = "<<end_tracer_lat<< endl ;*/

}	
//AB_e

      
		if (vel_field_timestep<1)
		{
		sprintf(ncfile, "%s%s_%04d%02d%02d_%02d.nc",ncdir.c_str(),prod_prefix.c_str(),date->tm_year, date->tm_mon+1, date->tm_mday, date->tm_hour);
		}
		else
		{
		sprintf(ncfile, "%s%04d/%s_%04d%02d%02d.nc",ncdir.c_str(),date->tm_year,prod_prefix.c_str(),date->tm_year, date->tm_mon+1, date->tm_mday);
		}
      if (verbose == 1)
	{
	  cout << " Reading " <<ncfile<<" -> ";
	}
      NcFile dataFile(ncfile, NcFile::ReadOnly);

      // Check to see if the file was opened.
      if(!dataFile.is_valid())
	return NC_ERR;

     if (!(Uvar = dataFile.get_var(u_name_nc.c_str())))
	return NC_ERR;
      if (!(Vvar  = dataFile.get_var(v_name_nc.c_str())))
	return NC_ERR;          

      for (i = 0; i < nlat; i++)
	{

	  if (!Uvar->set_cur(0, ilatmin+i, ilonmin))
	    return NC_ERR;
	  if (!Vvar->set_cur(0, ilatmin+i, ilonmin))
	    return NC_ERR;

	  if (!Uvar->get(&v_field_U[t][i][0], 1, 1, nlon))
	    return NC_ERR;
	  if (!Vvar->get(&v_field_V[t][i][0], 1, 1, nlon))
	    return NC_ERR;
	
	//Set nan values to 0
	 for (j = 0; j < nlon; j++)
	 {
	 if (isnan(v_field_U[t][i][j])==1)
	 	{
	 	v_field_U[t][i][j]=0;
	 	v_field_V[t][i][j]=0;
	 	}
	 }
	
	}

      dataFile.close();
      if (verbose == 1)
	{
	  cout << " OK" << endl;
	}
      

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // For debugging: prints the velocities

      // if(verbose==1)
      // 	{
      // 	  sprintf(name_out, "date_%04d-%02d-%02d.vfield", date->tm_year, date->tm_mon+1, date->tm_mday);
      // 	  fout.open(name_out);
      // 	  for(i=0; i<nlat; i++)
      // 	    {
      // 	      for(j=0; j<nlon; j++)
      // 		{
      // 		  if(v_field_U[t][i][j]<1e+20 && v_field_V[t][i][j]<1.e+20)
      // 		    {
      // 		      fout << v_field_lon[j] <<" "<<v_field_lat[i] <<" "<< v_field_U[t][i][j] <<" "<< v_field_V[t][i][j] << endl;
      // 		    }
      // 		}
      // 	    }
      // 	  fout.close();
      // 	}
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//cout << "u (1,1,"<<t<<")= "<<v_field_U[t][1][1]<<"  v (1,1,"<<t<<")= "<<v_field_V[t][1][1]<<endl;
    }
    
    ///////////////////////////////////////////////////////////////
    //TEST PART
    ///////////////////////////////////////////////////////////////
    if (opz_clim==1)
    {
    if (verbose==1){cout<<"Option climatology is active"<<endl;}
    double tmpU,tmpV,n_velfield2load_inv;
    n_velfield2load_inv=1./n_velfield2load;
    
    for (i = 0; i < nlat; i++)
    {
    for (j = 0; j < nlon; j++)
    {
    tmpU=0;tmpV=0;
    for(t=0; t<n_velfield2load; t++)
    {
    tmpU=tmpU+v_field_U[t][i][j]*n_velfield2load_inv;
    tmpV=tmpV+v_field_V[t][i][j]*n_velfield2load_inv;
	}//end t
	
	for(t=0; t<n_velfield2load; t++)
	{
	v_field_U[t][i][j]=tmpU;
	v_field_V[t][i][j]=tmpV;
	}//end t
	
	}//end j
	}//end i
	  
	}//end if
	///////////////////////////////////////////////////////////////
	//END TEST PART
	///////////////////////////////////////////////////////////////
  return 0;
}

void deallocate_vfield(int n_velfield2load)
{
  int i, j;
  //Deallocate

  delete[] v_field_lat;
  delete[] v_field_lon;
//AB  delete[] v_field_depth;

  for (i = 0; i < n_velfield2load; i++)
    {
      for (j = 0; j < nlat; j++)
	{
	  delete[] v_field_U[i][j];
	  delete[] v_field_V[i][j];
	}
      delete[] v_field_U[i];
      delete[] v_field_V[i];
    }
  delete[] v_field_U; 
  delete[] v_field_V;
}
