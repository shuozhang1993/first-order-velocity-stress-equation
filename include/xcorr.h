#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __TOOL_H__
#define __TOOL_H__
#include "tool.h"
#endif

#ifndef __XCORR_H__
#define __XCORR_H__

/*1D signal correlation*/
void xcorr_signal(float dt,struct matrix sig1,struct matrix sig2,struct matrix xcorr);

/*1D signal convolution*/
void xconv_signal(float dt,struct matrix sig1,struct matrix sig2,struct matrix xconv);

/*2D seismogram correlation by trace*/
void xcorr_seismogram(int hf_tau,float dt,struct matrix seis1,struct matrix seis2,struct matrix seis3);

/*2D seismogram convolution by trace*/
void xconv_seismogram(int hf_tau,float dt,struct matrix seis1,struct matrix seis2,struct matrix seis3);

#endif
