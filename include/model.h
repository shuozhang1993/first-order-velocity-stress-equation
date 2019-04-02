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
int N,npd,nabc;
int hori_grid,vert_grid;
int nx,nz,mx,mz,ix,iz;
float hori_inter,vert_inter,dx,dz;
int time_step,nt,it;
float time_interval,dt;
int shot_number,rece_number,nsrc,nrec,isrc,irec;
int freq_step,ifreq,nfreq;
float freq_low,freq_int,freq_cen;
int itera_number,iiter,niter;
char src_path[string_length],rec_path[string_length];
char vp_path[string_length],vpsm_path[string_length];
char rho_path[string_length],rhosm_path[string_length];
char vs_path[string_length],vssm_path[string_length];


/*model*/
struct matrix wavelet,cof;
struct matrix vp,padvp,rho,padrho,vs,padvs;
struct matrix p_signal;
struct matrix u,w,p,mem_u,mem_w,mem_p_x,mem_p_z;
struct matrix bu,lam;
struct matrix at,bt;
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
float weight(int i,int N);

/*Ricker wavelet*/
void ricker_wavelet(int strat_point,float frequency,float amplitude,struct matrix sub_wave);

/*update PML coefficient*/
void PML_update(float freq,struct matrix sub_at,struct matrix sub_bt);

/*Compute coefficient for PML boundary*/
void PML_coefficient(int power,float freq,struct matrix sub_at,struct matrix sub_bt);

/*Compute the physical property of the model*/
void field_transform(struct matrix sub_vp,struct matrix sub_rho);

/*solving first-order velocity-stress wave equation*/
void wave_eq_for_signal(char dir[],int isrc,struct matrix sub_vp,struct matrix sub_rho,struct matrix sub_wave,struct matrix seism);

#endif


