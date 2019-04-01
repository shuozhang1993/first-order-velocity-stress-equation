#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __FFT_H__
#define __FFT_H__

/*Get absolute value of complex number*/
void abs_complex(struct cpxmat in,struct matrix out);

/*Initialize memory space*/
void fft_zero(int l);

/*Fast forward fourier transform*/
void fft(int deminsion,struct matrix tnd, struct cpxmat snd);

/*Fast inverse fourier transform*/
void ifft(int deminsion,struct cpxmat snd,struct matrix tnd);

/*shift the frequency spectrum*/
void fftshift(struct cpxmat freq_be,struct cpxmat freq_af,char flag[]);
#endif
