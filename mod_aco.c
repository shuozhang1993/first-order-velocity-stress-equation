#ifndef __MPI_H__
#define __MPI_H__
#include <mpi.h>
#endif

#include "gene.h"
#include "model.h"

int main(int argc,char *argv[])
{
	MPI_set(argc,argv); int i,j,n,k,nf,ns,rece;char name[string_length];float freq; FILE *fp;

	coef_read("./model/mod/1_coef.txt"); model_arrang(); MPI_Barrier(MPI_COMM_WORLD);

	struct matrix p_gather; create_matrix(rece_number,1,time_step,&p_gather); 

/*-------------------------------Forward process start-------------------------------------*/

	for(nf=0;nf<freq_step;nf++)/*FREQUENCY LOOP*/
	{
		freq = freq_low + nf * freq_int; sou_let(freq); PML_coef(freq); if(rank==0) printf("\nComponent: %d, frequency: %.3fHz\n",nf+1,freq);

		for(ns=rank;ns<shot_number;ns+=size)/*SHOT LOOP*/
		{
			zero_propagation();/*zero propagation parameter*/

			for(n=0;n<time_step;n++)/*TIME LOOP*/
			{
				/*Add source wavelet*/ p.value[src_sta.value[ns*src_sta.y+0]*p.z+src_sta.value[ns*src_sta.y+1]] += wavelet.value[n];

				/*First order acoustic wave-equation*/PML_velo_acoustic(0); PML_stress_acoustic(0);

				for(rece=0;rece<rece_number;rece++)
					if(rec_sta.value[rece*rec_sta.y+0] != src_sta.value[ns*src_sta.y+0] || rec_sta.value[rece*rec_sta.y+1] != src_sta.value[ns*src_sta.y+1])
						p_gather.value[rece*p_gather.z+n] = p.value[rec_sta.value[rece*rec_sta.y+0]*p.z+rec_sta.value[rece*rec_sta.y+1]];
			}/*TIME LOOP END*/

			/*Seismic record*/ sprintf(name,"./model/sig/p%02d_%02d.obs",nf+1,ns+1);matrix_write(0,p_gather,name);

		}/*SHOT LOOP END*/PML_zero();
	}/*FREQUENCY LOOP END*/
	if(rank==0) {printf("\n-----------------------Fintie difference finishes----------------------------\n");}

	free_propagation(); free_model(); MPI_Finalize(); return(0);
}
