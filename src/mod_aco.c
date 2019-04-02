#ifndef __MPI_H__
#define __MPI_H__
#include <mpi.h>
#endif

#include "gene.h"
#include "tool.h"
#include "model.h"

int mute_direct=0;
int centra_freq=0;

int main(int argc,char *argv[])
{
	MPI_set(argc,argv); float freq;

	coef_read("./model/mod/2_coef.txt"); model_arrang(); MPI_Barrier(MPI_COMM_WORLD);

	struct matrix co_vp,co_rho,p_direct; create_matrix(nrec,1,1,1,nt,dt,&p_direct);
	create_matrix(mx,dx,1,1,mz,dz,&co_vp); create_matrix(mx,dx,1,1,mz,dz,&co_rho);

	for(ix=0;ix<mx;ix++) for(iz=0;iz<mz;iz++) co_vp.value[ix*mz+iz] = vp.value[ix*mz+1];
	for(ix=0;ix<mx;ix++) for(iz=0;iz<mz;iz++) co_rho.value[ix*mz+iz] = rho.value[ix*mz+1];

	if(centra_freq==1) nfreq=1;

	for(ifreq=0;ifreq<nfreq;ifreq++)/*FREQUENCY LOOP*/
        {   
                if(centra_freq==1) freq = freq_cen;
                else freq = freq_low + ifreq * freq_int;

                PML_update(freq,at,bt); ricker_wavelet(0,freq,200,wavelet);

                if(rank==0) printf("\nComponent: %03d, frequency: %02.3fHz\n",ifreq+1,freq);

                for(isrc=rank;isrc<nsrc;isrc+=size)/*SHOT LOOP*/
                {   
                        wave_eq_for_signal("forward",isrc,vp,rho,wavelet,p_signal);

                        if(mute_direct==1) wave_eq_for_signal("forward",isrc,co_vp,co_rho,wavelet,p_direct);

                        if(centra_freq==1) {sprintf(filename,"./model/sig/p_%02d.obs",isrc+1); matrix_write(0,p_signal,filename);}
                        else   {sprintf(filename,"./model/sig/p%02d_%02d.obs",ifreq+1,isrc+1); matrix_write(0,p_signal,filename);}

                        if(mute_direct==1)
                        {   
                                for(n=0;n<nrec*nt;n++) p_signal.value[n] -= p_direct.value[n];
                                if(centra_freq==0)
                                {   
                                        sprintf(filename,"./model/sig/p%02d_%02d.ref",ifreq+1,isrc+1);
                                                matrix_write(0,p_signal,filename);
                                }   
                                else
                                {   
                                        sprintf(filename,"./model/sig/p_%02d.ref",isrc+1);
                                                matrix_write(0,p_signal,filename);
                                }   
                        }   
                }/*SHOT LOOP END*/
        }/*FREQUENCY LOOP END*/
        if(rank==0) {printf("\n-----------------------Fintie difference finishes----------------------------\n");}

        delete_matrix(&co_vp); delete_matrix(&co_rho); delete_matrix(&p_direct);

        model_delete(); MPI_Finalize(); return(0);
}
