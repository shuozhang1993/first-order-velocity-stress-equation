#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __MODEL_H__
#define __MODEL_H__

/*MPI*/
int size,rank;

/*coefficient*/
int hori_grid,vert_grid,time_step;
int nx,nz,N,npd,x,z,boundary;
int shot_number,rece_number,freq_step,itera_number;
float freq_low,freq_int,freq_cen;
float lenth_phy,depth_phy,dt,time_interval,dh,hori_inter,vert_inter;
char src_path[string_length],rec_path[string_length];
char vp_path[string_length],vs_path[string_length],rho_path[string_length];

/*model*/
float vpmax;
struct matrix wavelet,cof;
struct matrix Vp,Rho,vp,rho,BU,LAM;
struct matrix u,w,p,phi_u,phi_w,phi_p_x,phi_p_z;
struct matrix ax,az,bx,bz;
struct location src_sta,rec_sta;

/*MPI*/
void MPI_set(int argc,char *argv[]);

/*read coefficient from txt file*/
void coef_read(char file_path[string_length]);

/*allocate memory and set model value*/
void model_arrang();

/*zero propagation parameter*/
void zero_propagation();

/*free memory space for propagation parameter*/
void free_propagation();

/*free memory for model value*/
void free_model();

/*PML setting*/
void PML_coef(float freq);

/*zero PML parameter*/
void PML_zero();

/*Source function*/
float ricker(float freq, float time);

/*Source wavelet*/
void sou_let(float freq);

/*Weighted coefficiency*/
float weight(int i, int N);

/*Stable adjudgement*/
int stability_condition();

/*Compute the physical property of the model*/
void field_transform();

/*Extend boundary area*/
void extend(struct matrix in,struct matrix out);

/*Compute coefficient for PML boundary*/
void PML_coefficient(int power,float freq);

/*Update stress field in acoustic case with PML*/
void PML_stress_acoustic(int flag);

/*Update displacement field in acoustic case*/
void PML_velo_acoustic(int flag);

#endif


