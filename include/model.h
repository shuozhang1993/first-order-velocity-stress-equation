#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __TOOL_H__
#define __TOOL_H__
#include "tool.h"
#endif

#ifndef __MODEL_H__
#define __MODEL_H__

/*MPI*/
int size,rank;

/*coefficient*/
int hori_grid,vert_grid,time_step;
int nx,nz,N,npd,mx,mz,nabc,nt,it;
int nsrc,isrc,shot_number,nrec,irec,rece_number;
int nfreq,ifreq,freq_step,niter,itera_number;
float freq_low,freq_int,freq_cen;
float dt,time_interval,dx,dz,hori_inter,vert_inter;
char src_path[string_length],rec_path[string_length];
char vp_path[string_length],sm_path[string_length];
char vs_path[string_length],rho_path[string_length];

/*model*/
struct matrix wavelet,cof;
struct matrix vp,padvp;
struct matrix p_sig;
struct matrix p_now,p_new,p_old,lapla;
struct matrix at,bt,mem_x,mem_z;
struct location src_ori,src_sta,rec_ori,rec_sta;

/*MPI*/
void MPI_set(int argc,char *argv[]);

/*read coefficient from txt file*/
void coef_read(char file_path[string_length]);

/*allocate memory and set model value*/
void model_arrang();

/*release memory space*/
void model_delete();

/*Weighted coefficiency*/
float weight(int i, int N);

/*Ricker wavelet*/
void ricker_wavelet(int start_point,float frequency,float amplitude,struct matrix sub_wave);

/*update PML coefficient*/
void PML_update(float freq,struct matrix sub_at,struct matrix sub_bt);

/*Compute coefficient for PML boundary*/
void PML_coefficient(int power,float freq,struct matrix sub_at,struct matrix sub_bt);

/*solving wave equation for seismogram*/
void wave_eq_for_signal(int isrc,struct matrix sub_vp,struct matrix sub_wave,struct matrix seism);

#endif


