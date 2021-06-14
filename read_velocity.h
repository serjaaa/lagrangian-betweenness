/* FUNCTIONS */

int read_vfield(struct tm start_date, int n_velfield2load, int time_dir);//AB, int layer_index);// This function read the 2D velocity field since start_date to start_date+tau in a constant depth (layer_index) 

void deallocate_vfield(int n_velfield2load);// This function free the memory of V_field variables. It is necesary to call this function before the end of the program

/* VARIABLES */

extern double* v_field_lat; 
extern double* v_field_lon;
extern double* v_field_depth;
extern double*** v_field_U;
extern double*** v_field_V;
extern int nlat;
extern int nlon;
extern int ndepth;

//AB_s parameters added to use different products of vel field
extern double start_tracer_long;
extern double end_tracer_long;
extern double start_tracer_lat;
extern double end_tracer_lat;
