/* Compartments */
enum{Sc, Ic, Cc, Rc, Sy, Iy, Cy, Ry, Sa, Ia, Ca, Ra};

/* Global parameters */
enum{UPSILON_C, UPSILON_Y, UPSILON_A, 
     GAMMA_C, GAMMA_Y, GAMMA_A, 
     CHI_C, CHI_Y, CHI_A, 
     RHO_C, RHO_Y, RHO_A,
     PSI_C, PSI_Y, PSI_A,
     ALPHA_C, ALPHA_Y, ALPHA_A, SIGMA, 
     BETA_T1, BETA_T2, BETA_T3, BETA_T4};

/* Offsets in node local data (ldata) to parameters in the model */
enum {COUPLING, END_T1, END_T2, END_T3, END_T4, NEIGHBOR};

/* Environmental infectious pressure */
enum{PHI};

/* Decay of environmental infectious pressure with a forward Euler step.
 *
 * The time dependent beta is divided into four intervals of the year
 * where 0 <= day < 365
*/
double SimInf_forward_euler_linear_decay(
    double phi, int day,
    int end_t1, int end_t2, int end_t3, int end_t4,
    double beta_t1, double beta_t2, double beta_t3, double beta_t4);

/* Local spread of the environmental infectious pressure phi among
 * proximal nodes.
*/
double SimInf_local_spread(
    const double *neighbors,
    const double *phi,
    const int *u,
    const double N_i,
    const double phi_i,
    const int Nc,
    const double D);

/* Post-time Step function */
int ptsFun(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t);





