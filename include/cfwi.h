#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __TOOL_H__
#define __TOOL_H__
#include "tool.h"
#endif

#ifndef __CFWI_H__
#define __CFWI_H__

/*MPI*/
int size,rank;

/*coefficient*/
int hori_grid,vert_grid,time_step;
int nx,nz,N,npd,mx,mz,nabc,nt,it;
int nsrc,isrc,shot_number,nrec,irec,rece_number;
int nfreq,ifreq,freq_step,iiter,niter,itera_number;
float freq_low,freq_int,freq_cen;
float dt,time_interval,dx,dz,hori_inter,vert_inter;
char src_path[string_length],rec_path[string_length];
char vp_path[string_length],sm_path[string_length];
char vs_path[string_length],rho_path[string_length];

/*model*/
struct matrix at,bt;
struct matrix vp,padvp,truvp,lapla;
struct matrix p_now,p_new,p_old;
struct matrix mem_p_x,mem_p_z;
struct matrix q_now,q_new,q_old;
struct matrix mem_q_x,mem_q_z;
struct matrix p_sig,p_sum,p_grd,p_dot;
struct matrix p_obs,p_ajs,p_syn;
struct matrix p_t,p_b,p_l,p_r,ck_old,ck_now;
struct matrix gk0,dk0;
struct location src_ori,src_sta,rec_ori,rec_sta;
struct matrix wavelet,cof,misfit_dat,misfit_mod;

/*MPI*/
void MPI_set(int argc,char *argv[]);

/*read coefficient from txt file*/
void coef_read(char file_path[string_length]);

/*allocate memory and set model value*/
void model_arrang();

/*release memory space*/
void model_delete();

/*Weighted coefficiency*/
float weight(int sub_i, int sub_N);

/*Ricker wavelet*/
void ricker_wavelet(int start_point,float frequency,float amplitude,struct matrix sub_wave);

/*update PML coefficient*/
void PML_update(float freq,struct matrix sub_at,struct matrix sub_bt);

/*Compute coefficient for PML boundary*/
void PML_coefficient(int power,float freq,struct matrix sub_at,struct matrix sub_bt);

/*solving wave equation for seismogram*/
void wave_eq_for_signal(int isrc,struct matrix sub_vp,struct matrix sub_wave,struct matrix seism);

/*solving wave equation for reconstruction*/
void wave_eq_for_recons(int isrc,
			struct matrix sub_vp,struct matrix sub_wave,struct matrix seism,
			struct matrix sub_ck_now,struct matrix sub_ck_old,
			struct matrix sub_p_t,struct matrix sub_p_b,struct matrix sub_p_l,struct matrix sub_p_r);

/*compute adjoint source*/
void adjoint_source(struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_p_ajs);

/*cross-correlation imaging condition*/
void banana_doughnut(int isrc,
		     struct matrix sub_vp,struct matrix sub_wave,
		     struct matrix sub_ck_now,struct matrix sub_ck_old,
		     struct matrix sub_p_t,struct matrix sub_p_b,struct matrix sub_p_l,struct matrix sub_p_r,
		     struct matrix sub_p_ajs,struct matrix sub_p_sig);

/*velocity kernel for cfwi*/
float cfwi_velo_kernel(int sub_ifreq,int sub_iiter,struct matrix sub_vp,struct matrix sub_wave,struct matrix sub_p);

/*velocity kernel for cfwi with adaptive matching filter*/
float amf_cfwi_velo_kernel(int sub_ifreq,int sub_iiter,struct matrix sub_vp,struct matrix sub_wave,struct matrix sub_p);

/*update velocity model*/
void update_model_for_cfwi(int sub_ifreq,int sub_iiter,float resi_value,
			   struct matrix sub_vp,struct matrix sub_wave,
			   struct matrix gradient,struct matrix sub_gk0,struct matrix sub_dk0);
/*update velocity model with amf*/
void amf_update_model_for_cfwi(int sub_ifreq,int sub_iiter,float resi_value,
			       struct matrix sub_vp,struct matrix sub_wave,
			       struct matrix gradient,struct matrix sub_gk0,struct matrix sub_dk0);

#endif


