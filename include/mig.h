#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __TOOL_H__
#define __TOOL_H__
#include "tool.h"
#endif

#ifndef __MIG_H__
#define __MIG_H__

/*MPI*/
int size,rank;

/*coefficient*/
int hori_grid,vert_grid,time_step;
int nx,nz,N,npd,mx,mz,nabc,nt,it;
int nsrc,isrc,shot_number,nrec,irec,rece_number;
int nfreq,ifreq,dfreq,freq_step;
int iiter,niter,itera_number;
float freq_low,freq_int,freq_cen;
float dt,time_interval,dx,dz,hori_inter,vert_inter;
char src_path[string_length],rec_path[string_length];
char vp_path[string_length],sm_path[string_length];
char vs_path[string_length],rho_path[string_length];

/*model*/
struct matrix at,bt;
struct matrix vp,padvp,truvp;
struct matrix pe,padpe;

struct matrix lapla_bac,lapla_per;

struct matrix p_src_now,p_src_new,p_src_old;
struct matrix mem_p_src_x,mem_p_src_z;

struct matrix p_rec_now,p_rec_new,p_rec_old;
struct matrix mem_p_rec_x,mem_p_rec_z;

struct matrix q_src_now,q_src_new,q_src_old;
struct matrix mem_q_src_x,mem_q_src_z;

struct matrix p_src_sig,p_rec_sig;
struct matrix p_src_dot,p_rec_dot;
struct matrix p_sum,p_grd;
struct matrix p_obs,p_ajs,p_syn;

struct matrix p_t,p_b,p_l,p_r,ck_p_old,ck_p_now;
struct matrix q_t,q_b,q_l,q_r,ck_q_old,ck_q_now;

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

/*solving born equation for seismogram*/
void born_eq_for_signal(int isrc,struct matrix sub_vp,struct matrix sub_pe,struct matrix sub_wave,struct matrix seism);

/*solving wave equation for reconstruction*/
void wave_eq_for_recons(int isrc,
			struct matrix sub_vp,struct matrix sub_wave,struct matrix seism,
			struct matrix sub_ck_p_now,struct matrix sub_ck_p_old,
			struct matrix sub_p_t,struct matrix sub_p_b,struct matrix sub_p_l,struct matrix sub_p_r);

/*solving born equation for reconstruction*/
void born_eq_for_recons(int isrc,
			struct matrix sub_vp,struct matrix sub_pe,struct matrix sub_wave,struct matrix seism,
			struct matrix sub_ck_p_now,struct matrix sub_ck_p_old,
			struct matrix sub_p_t,struct matrix sub_p_b,struct matrix sub_p_l,struct matrix sub_p_r,
			struct matrix sub_ck_q_now,struct matrix sub_ck_q_old,
			struct matrix sub_q_t,struct matrix sub_q_b,struct matrix sub_q_l,struct matrix sub_q_r);

/*compute adjoint source*/
void adjoint_source(struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_p_ajs);

/*data misfit with signal substraction*/
float misfit_sub(struct matrix sub_p_syn,struct matrix sub_p_obs);

/*zero-lag cross-correlation imaging condition*/
void cross_correlation(int isrc,
			struct matrix sub_vp,struct matrix sub_wave,
			struct matrix sub_ck_p_now,struct matrix sub_ck_p_old,
			struct matrix sub_p_t,struct matrix sub_p_b,struct matrix sub_p_l,struct matrix sub_p_r,
			struct matrix sub_p_ajs,struct matrix sub_p_src_sig);

/*reverse time migration*/
void reverse_time_migration(int sub_iiter,float p_1,float p_2,float p_3,float p_4,
			    struct matrix sub_vp,struct matrix sub_wave,struct matrix rtm_refl);

/*reflectivity kernel for lsrtm*/
float least_square_kernel(int sub_iiter,float p_1,float p_2,float p_3,float p_4,
			  struct matrix sub_vp,struct matrix sub_pe,struct matrix sub_wave,struct matrix lsm_refl);

/*least_square reverse time migration*/
void least_square_migration(float p1,float p2,float p3,float p4,
			    struct matrix sub_vp,struct matrix sub_wave,struct matrix lsm_refl);

/**update velocity mode for lsm*/
void update_image_for_lsm(int sub_iiter,float p1,float p2,float p3,float p4,float resi_value,
			   struct matrix sub_vp,struct matrix sub_pe,struct matrix sub_wave,
			   struct matrix gradient,struct matrix sub_gk0,struct matrix sub_dk0);

#endif


