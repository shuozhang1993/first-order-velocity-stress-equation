#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __TOOL_H__
#define __TOOL_H__

/*Gaussian Filter*/
void bandpassfilter(float _1st,float _2nd,float _3rd,float _4th,struct filter lowpass);

/*laplacian operator remove low wavenumber component*/
void lapla_filter_z(struct matrix in,struct matrix out);

/*filter seismogram*/
void trace_filter(int sub_isrc,int hori_len,int len_1,int len_2,int len_3,int len_4,
                  struct matrix in,struct location src_loc,struct location rec_loc);

/*Remove the direct wave*/
void remove_direct(int n_shot,int len,int duration,float delta_t,float delta_h,float velo,
		   struct location src_loc,struct location rec_loc,struct matrix gather);

/*Data misfit*/
float abs_misfit_value(struct matrix mat1,struct matrix mat2);

/*timeshift signal*/
void timeshift(int shift,int h_dura,struct matrix trace,struct matrix volume);

#endif
