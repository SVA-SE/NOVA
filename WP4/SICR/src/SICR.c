#include "SimInf.h"
#include "SICR.h"

/* Decay of environmental infectious pressure with a forward Euler step.
 *
 * The time dependent beta is divided into four intervals of the year
 * where 0 <= day < 365
*/

double SimInf_forward_euler_linear_decay(
    double phi, int day,
    int end_t1, int end_t2, int end_t3, int end_t4,
    double beta_t1, double beta_t2, double beta_t3, double beta_t4)
{
  if (day < end_t2) {
    if (day < end_t1) {
      if (end_t1 < end_t4)
        return phi * (1.0 - beta_t1);
      if (day < end_t4) {
        if (end_t4 < end_t3)
          return phi * (1.0 - beta_t4);
        if (day < end_t3)
          return phi * (1.0 - beta_t3);
        return phi * (1.0 - beta_t4);
      }
      return phi * (1.0 - beta_t1);
    }
    
    return phi * (1.0 - beta_t2);
  }
  
  if (end_t3 < end_t1 || day < end_t3)
    return phi * (1.0 - beta_t3);
  
  if (end_t4 < end_t1 || day < end_t4)
    return phi * (1.0 - beta_t4);
  
  return phi * (1.0 - beta_t1);
}


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
    const double D)
{
  int j, k;
  double N_j, ls = 0.0;
  const double phi_i_N_i = phi_i * N_i;
  
  j = (int)*neighbors++;
  while (j >= 0) {
    /* Count number of individuals in node j */
    for (k = j * Nc, N_j = 0; k < (j + 1) * Nc; k++)
      N_j += u[k];
    
    if (N_j > 0.0)
      ls += ((phi[j] * N_j - phi_i_N_i) * D) / (N_i * (*neighbors));
    
    /* Move to next neighbor pair (index, distance) */
    neighbors++;
    j = (int)*neighbors++;
  }
  
  return ls;
}


/* post time step function */
int ptsFun(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t)
{
  
  const int day = (int)t % 365;
  const double I_c = u[Ic];
  const double I_y = u[Iy];
  const double I_a = u[Ia];
  const double C_c = u[Cc];
  const double C_y = u[Cy];
  const double C_a = u[Ca];
  const double S_i = u[Sc] + u[Sy] + u[Sa];
  const double I_i = I_c + I_y + I_a;
  const double C_i = C_c + C_y + C_a;
  const double R_i = u[Rc] + u[Ry] + u[Ra];
  const double N_i = S_i + I_i + C_i + R_i;
  const double phi = v[PHI];
  const int Nc = 12;
  
  /* Deterimine the pointer to the continuous state vector in the
  * first node. Use this to find phi at neighbours to the current
  * node. */
  const double *phi_0 = &v[-node];
  
  /* Deterimine the pointer to the compartment state vector in the
   * first node. Use this to find the number of individuals at
   * neighbours to the current node. */
  const int *u_0 = &u[-Nc*node];
  
  /* Time dependent beta in each of the four intervals of the
   * year. Forward Euler step. */
  v_new[PHI] = SimInf_forward_euler_linear_decay(
    phi, day,
    ldata[END_T1], ldata[END_T2], ldata[END_T3], ldata[END_T4],
    gdata[BETA_T1], gdata[BETA_T2], gdata[BETA_T3], gdata[BETA_T4]);
  
  /* Local spread among proximal nodes. */
  if (N_i > 0.0) {
    v_new[PHI] += ((gdata[ALPHA_C] * (I_c + gdata[SIGMA]*C_c)) +
                   (gdata[ALPHA_Y] * (I_y + gdata[SIGMA]*C_y)) +
                   (gdata[ALPHA_A] * (I_a + gdata[SIGMA]*C_a)) ) / N_i +
      SimInf_local_spread(&ldata[NEIGHBOR], phi_0, u_0,
                          N_i, phi, Nc, ldata[COUPLING]);
  }
  
  if (!isfinite(v_new[PHI]))
    return SIMINF_ERR_V_IS_NOT_FINITE;
  if (v_new[PHI] < 0.0)
    return SIMINF_ERR_V_IS_NEGATIVE;
  return phi != v_new[PHI]; /* 1 if needs update */
}
  
 
