#include "gene.h"
#include "tool.h"

/*Gaussian Filter*/
void bandpassfilter(float _1st,float _2nd,float _3rd,float _4th,
					struct filter lowpass)
{
	if(_1st*_2nd*_3rd*_4th<0)
	{
		printf("\nBandPassFilter: The coefficients should be negative.\n"); return;
	}
	if(!(_4th>=_3rd && _3rd>=_2nd && _2nd>=_1st)) 
	{
		printf("\nBandPleaseFilter : Please set the coefficients in order.\n");
		printf("\n1st:%e, 2nd:%e, 3rd:%e, 4th:%e\n",_1st,_2nd,_3rd,_4th); return;
	}

    int p1st,p2nd,p3rd,p4th,dif_p1,dif_p2,dif_p3;

    p1st=(int)(_1st/lowpass.interval);p2nd=(int)(_2nd/lowpass.interval);
    p3rd=(int)(_3rd/lowpass.interval);p4th=(int)(_4th/lowpass.interval);

    dif_p1=p2nd-p1st; dif_p2=p3rd-p2nd; dif_p3=p4th-p3rd;
    for(n=0;n<lowpass.length;n++)
    {
    	if(n<p1st) lowpass.value[n] = 0.0;

    	if(n>=p1st && n<p2nd) lowpass.value[n] = (cos((PI/dif_p1)*(p2nd-n))+1)/2;

        if(n>=p2nd && n<p3rd) lowpass.value[n] = 1.0;

        if(n>=p3rd && n<p4th) lowpass.value[n] = (cos((PI/dif_p3)*(n-p3rd))+1)/2;

        if(n>=p4th) lowpass.value[n] = 0.0;
    }
}

/*laplacian operator remove low wavenumber component*/
void laplace_filter_z(struct matrix in,struct matrix out)
{
	struct matrix trace,lapla; float max_in,max_la;
	int ix,x=in.x,iz,z=in.z,iy,y=in.y;

	create_matrix(1,1,1,1,z,1,&trace);
	create_matrix(1,1,1,1,z,1,&lapla);

	for(ix=0;ix<x;ix++)
	{
		for(iz=0;iz<z;iz++) trace.value[iz] = in.value[ix*z+iz]; max_in = matrix_max(trace);
		for(iz=1;iz<z-1;iz++) lapla.value[iz] = -(trace.value[iz+1]+trace.value[iz-1]-2*trace.value[iz]);
		max_la = matrix_max(lapla)+max_in/1e10; for(iz=0;iz<z;iz++) lapla.value[iz] *= max_in/max_la;
		for(iz=0;iz<z;iz++) out.value[ix*z+iz] = lapla.value[iz];
	}

	delete_matrix(&lapla);
	delete_matrix(&trace);
}

/*filter seismogram*/
void trace_filter(int sub_isrc,int hori_len,int len_1,int len_2,int len_3,int len_4,
                  struct matrix in,struct location src_loc,struct location rec_loc)
{
        float src_mx=src_loc.value[sub_isrc*2];
        int irec,nrec=in.x,it,nt=in.z;hori_len++;
        struct filter taper_hf,taper_al;

        create_filter(hori_len,1.0,&taper_hf); create_filter(hori_len,1.0,&taper_al);
        bandpassfilter(len_1,len_2,len_3,len_4,taper_hf);

        for(i=0;i<src_mx;i++) taper_al.value[i] = taper_hf.value[abs(i-src_mx)];
        for(i=src_mx;i<hori_len;i++) taper_al.value[i] = taper_hf.value[abs(i-src_mx)];

        for(it=0;it<nt;it++) for(irec=0;irec<nrec;irec++)
                in.value[irec*nt+it] *= taper_al.value[rec_loc.value[irec*2]];

        delete_filter(&taper_hf); delete_filter(&taper_al);
}


/*Remove the direct wave*/
void remove_direct(int n_shot,int len,int duration,float delta_t,float delta_x,float velo,
				   struct location src_loc,struct location rec_loc,struct matrix gather)
{
	int i,j,first; struct filter removefilter; FILE *fp; create_filter(gather.z+len+duration+1,1,&removefilter);

	for(i=0;i<gather.x;i++)
	{
		first = (int)fabsf((src_loc.value[n_shot*src_loc.y]-rec_loc.value[i*rec_loc.y])*delta_x/(velo*delta_t));

		if(first>=gather.z) first = gather.z; 
		
		bandpassfilter(first+len,first+len+duration,removefilter.length,removefilter.length,removefilter);

		for(j=0;j<gather.z;j++) gather.value[i*gather.z+j] *= removefilter.value[j];
	}

	delete_filter(&removefilter);
}

/*Data misfit*/
float abs_misfit_value(struct matrix mat1,struct matrix mat2)
{
	int n; float value = 0.0; 
	for(n=0;n<mat1.x*mat1.z*mat1.y;n++) value += fabs(mat1.value[n]-mat2.value[n]); return value;
}

/*time shift signal*/
void timeshift(int shift,int h_dura,struct matrix trace,struct matrix volume)
{
        int n,it,nt=trace.z;
        int a_dura = h_dura*2+1,delay;

        /*left*/
        for(n=0;n<h_dura;n++)
        {
                delay = (h_dura-n)*shift;
                for(it=0;it<(nt-delay);it++) volume.value[n*nt+it] = trace.value[it+delay];
        }

        /*middel*/for(it=0;it<nt;it++) volume.value[h_dura*nt+it] = trace.value[it];

        /*right*/
        for(n=h_dura+1;n<a_dura;n++)
        {
                delay = (n-h_dura)*shift;
                for(it=delay;it<nt;it++) volume.value[n*nt+it] = trace.value[it-delay];
	}
}
