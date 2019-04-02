#ifndef __FFTW_H__
#define __FFTW_H__
#include <fftw3.h>
#endif

#include "gene.h"
#include "fft.h"

fftw_complex *in,*out;
fftw_plan p;

/*Get absolute value of complex number*/
void abs_complex(struct cpxmat in,struct matrix out)
{
	int n; for(n=0;n<in.x*in.y*in.z;n++) out.value[n] = sqrt(in.real[n]*in.real[n]+in.imag[n]*in.imag[n]);
}

/*Initialize memory space*/
void fft_zero(int l)
{
	int i,j; for(i=0;i<l;i++) for(j=0;j<2;j++) {in[i][j] = 0.; out[i][j]= 0.;}
}

/*Fast forward fourier transform*/
void fft(int deminsion,struct matrix tnd, struct cpxmat snd)
{
	int i,j,n,len=tnd.x*tnd.y*tnd.z;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*len);
	out= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*len);

	fft_zero(len);

	if(deminsion-1==0) p = fftw_plan_dft_1d(len,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

	if(deminsion-2==0) p = fftw_plan_dft_2d(tnd.x,tnd.z,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

	for(n=0;n<len;n++) {in[n][0] = tnd.value[n];}

	fftw_execute(p);

	for(n=0;n<len;n++) {snd.real[n] = out[n][0]; snd.imag[n] = out[n][1];}

	fftw_destroy_plan(p); fftw_free(in); fftw_free(out);
}

/*Fast inverse fourier transform*/
void ifft(int deminsion,struct cpxmat snd,struct matrix tnd)
{
	int n,len=snd.x*snd.z*snd.y;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*len);
	out= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*len);

	fft_zero(len);

	if(deminsion-1==0) p = fftw_plan_dft_1d(len,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);

	if(deminsion-2==0) p = fftw_plan_dft_2d(snd.x,snd.z,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);

	for(n=0;n<len;n++) {in[n][0] = snd.real[n]; in[n][1] = snd.imag[n];}

	fftw_execute(p);

	for(n=0;n<len;n++) {tnd.value[n] = out[n][0]/len;}

	fftw_destroy_plan(p); fftw_free(in); fftw_free(out);
}

/*shift the negative frequency component*/
void fftshift(struct cpxmat freq_be,struct cpxmat freq_af,char flag[])
{
	int len=freq_be.x*freq_be.y*freq_be.z;
	int n,wdur,fdur;

	if(len%2==0) {wdur=len+1; fdur=len/2;}
	else {wdur=len; fdur=(len-1)/2;}
	
	for(n=0;n<fdur;n++) {freq_af.real[n]=freq_be.real[fdur+1+n];freq_af.imag[n]=freq_be.imag[fdur+1+n];}
	freq_af.real[fdur]=freq_be.real[fdur]; freq_af.imag[fdur]=freq_be.imag[fdur];
	for(n=0;n<fdur;n++) {freq_af.real[fdur+1+n]=freq_be.real[n];freq_af.imag[fdur+1+n]=freq_be.imag[n];}
}

