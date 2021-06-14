int readparams(string nfileparameters);

/* PARAMATERS */

extern  int domain_network;
extern  int land_sea_norm; 
extern  int vel_product;
extern  int verbose;       
extern  double node_size;  
extern  double particle_latspacing;
extern  double delta0_ftle;
extern  double delta_fsle_fin;
extern  int trj_timestep;
extern  int rand_idx;
extern  int opz_clim;
extern  int layer_index;
extern  struct tm start_date;
extern  double tau; 
extern  double int_step;

extern double degree_resolution,degree_resolution_inv,degree_resolution2_in_meters_inv;

extern double network_ll_lat;
extern double network_ll_lon;
extern double network_tr_lat;
extern double network_tr_lon;

extern double vfield_ll_lat;
extern double vfield_ll_lon;
extern double vfield_tr_lat;
extern double vfield_tr_lon;

//AB_s parameters added to use different products of vel field
extern string ncdir,u_name_nc,v_name_nc;
extern string lon_name_nc,lat_name_nc;
extern string prod_name,prod_prefix;
extern double buffer_particles;

extern double lon_vfield_min;
extern double lon_vfield_max;
extern double lat_vfield_min;
extern double lat_vfield_max;
extern int periodic_conditionWE_YN;
extern double long_corr;
extern double int_step;
extern int time_dir;
extern double vel_field_timestep;

//AB_e
