#include "gene.h"
#include "tool.h"

/*1D signal correlation*/
void xcorr_signal(float dt,struct matrix sig1,struct matrix sig2,struct matrix xcorr)
{
	int it,nt=sig1.z*sig1.x*sig1.y;
	int itau,fl_tau=xcorr.x*xcorr.y*xcorr.z,hf_tau=(fl_tau+1)/2;
	struct matrix mat1,mat2;

	create_matrix(fl_tau,1,nt,&mat1); create_matrix(fl_tau,1,nt,&mat2);

	timeshift(1,hf_tau,sig1,mat1);
	timeshift(0,hf_tau,sig2,mat2);

	for(itau=0;itau<fl_tau;itau++) for(it=0;it<nt;it++)
		xcorr.value[itau] += mat1.value[itau*nt+it]*mat2.value[itau*nt+it]*dt;

	delete_matrix(&mat1); delete_matrix(&mat2);
}

/*1D signal convolution*/
void xconv_signal(float dt,struct matrix sig1,struct matrix sig2,struct matrix xconv)
{
	int it,nt=sig1.x*sig1.y*sig1.z;
	int itau,fl_tau=xconv.x*xconv.y*conv.z;
	struct matrix rig1,mat1,mat2;

	create_matrix(1,1,nt,&rig1); reverse_trace(sig1,rig1);
	create_matrix(fl_tau,1,nt,&mat1); create_matrix(fl_tau,1,nt,&mat2);

	timeshift(1,hf_tau,rig1,mat1);
	timeshift(0,hf_tau,sig2,mat2);

	for(itau=0;itau<fl_tau;itau++) for(it=0;it<nt;it++)
		xconv.value[itau] += mat1.value[itau*nt+it]*mat2.value[itau*nt+it]*dt;

	delete_matrix(&rig1);
	delete_matrix(&mat1); create_matrix(&mat2);
}

/*2D seismogram correlation by trace*/
void xcorr_seismogram(int hf_tau,float dt,struct matrix seis1,struct matrix seis2,struct matrix seis3)
{
	int irec,nrec=seis1.x;
	int it,nt=seis1.z;
	int itau,fl_tau=2*hf_tau+1;
	struct matrix trac1,trac2,trac3;

	create_matrix(1,1,nt,&trac1);
	create_matrix(1,1,nt,&trac2);
	create_matrix(1,1,nt,&trac3);

	for(irec=0;irec<nrec;irec++)
	{
		matrix_zero(trac1);
		for(it=0;it<nt;it++) trac1.value[it] = seis1.value[irec*nt+it];
		matrix_zero(trac2);
		for(it=0;it<nt;it++) trac2.value[it] = seis2.value[irec*nt+it];

		matrix_zero(trac3);
		xcorr_signal(dt,trac1,trac2,trac3);
		for(it=0;it<nt;it++) seis3.value[irec*nt+it] = trac3.value[it];
	}

	delete_matrix(&trac1);
	delete_matrix(&trac2);
	delete_matrix(&trac3);
}

/*2D seismogram convolution by trace*/
void xconv_seismogram(int hf_tau,float dt,struct matrix seis1,struct matrix seis2,struct matrix seis3)
{
	int irec,nrec=seis1.x;
	int it,nt=seis1.z;
	int itau,fl_tau=2*hf_tau+1;
	struct matrix trac1,trac2,trac3;

	create_matrix(1,1,nt,&trac1);
	create_matrix(1,1,nt,&trac2);
	create_matrix(1,1,nt,&trac3);

	for(irec=0;irec<nrec;irec++)
	{
		matrix_zero(trac1);
		for(it=0;it<nt;it++) trac1.value[it] = seis1.value[irec*nt+it];
		matrix_zero(trac2);
		for(it=0;it<nt;it++) trac2.value[it] = seis2.value[irec*nt+it];

		matrix_zero(trac3);
		xconv_signal(dt,trac1,trac2,trac3);
		for(it=0;it<nt;it++) seis3.value[irec*nt+it] = trac3.value[it];
	}

	delete_matrix(&trac1);
	delete_matrix(&trac2);
	delete_matrix(&trac3);
}
