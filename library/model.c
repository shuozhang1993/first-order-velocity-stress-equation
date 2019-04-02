#include "gene.h"
#include "tool.h"
#include "model.h"

#ifndef ___MPI_H__
#define ___MPI_H__
#include <mpi.h>
#endif

const int diff_order=2;/*order of differential*/
const int abs_width=200;/*width of absorption boundary*/

/*#####################memory arrangement####################*/

/*MPI*/
void MPI_set(int argc,char *argv[])
{
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
}

/*read coefficient from txt file*/
void coef_read(char file_path[])
{
	int count; FILE *fp=fopen(file_path,"rt");

	if(rank==0 && fp==NULL) {printf("\nfail to read the coefficient text file.\n"); exit(1);}

	count=fscanf(fp,"%d\n",&hori_grid);   count=fscanf(fp,"%d\n",&vert_grid);
	count=fscanf(fp,"%f\n",&hori_inter);  count=fscanf(fp,"%f\n",&vert_inter);

	count=fscanf(fp,"%d\n",&time_step);   count=fscanf(fp,"%f\n",&time_interval);

	count=fscanf(fp,"%d\n",&shot_number); count=fscanf(fp,"%d\n",&rece_number);

	count=fscanf(fp,"%d\n",&freq_step);   count=fscanf(fp,"%f\n",&freq_int);
	count=fscanf(fp,"%f\n",&freq_low);    count=fscanf(fp,"%f\n",&freq_cen);

	count=fscanf(fp,"%d\n",&itera_number);

	count=fscanf(fp,"%s\n",src_path);     count=fscanf(fp,"%s\n",rec_path);
	count=fscanf(fp,"%s\n",vp_path);      count=fscanf(fp,"%s\n",vpsm_path);

	if(rank==0 && fp!=NULL) printf("\nSuccessfully read the coefficient from text file.\n");

	fclose(fp); fp = NULL;

	N = diff_order/2; npd = abs_width; nabc = abs_width/2;
	mx = hori_grid; nx = mx+2*nabc; dx = hori_inter;
	mz = vert_grid; nz = mz+2*nabc; dz = vert_inter;
	nt = time_step; dt = time_interval;
	nfreq = freq_step;
	nsrc = shot_number;
	nrec = rece_number;
}

/*allocate memory and set model value*/
void model_arrang()
{
	create_matrix(mx,dx,1,1,mz,dz,&vp);  create_matrix(nx,dx,1,1,nz,dz,&padvp);
	create_matrix(mx,dx,1,1,mz,dz,&rho); create_matrix(nx,dx,1,1,nz,dz,&padrho);

	create_matrix(nx,dx,1,1,nz,dz,&bu);  create_matrix(nx,dx,1,1,nz,dz,&lam);

	create_matrix(nx,dx,1,1,nz,dz,&u);   create_matrix(nx,dx,1,1,nz,dz,&w); create_matrix(nx,dx,1,1,nz,dz,&p);

	create_matrix(nx,1,1,1,nz,1,&at); create_matrix(nx,1,1,1,nz,1,&bt);

	create_matrix(nx,dx,1,1,nz,dz,&mem_u); create_matrix(nx,dx,1,1,nz,dz,&mem_w);
	create_matrix(nx,dx,1,1,nz,dz,&mem_p_x); create_matrix(nx,dx,1,1,nz,dz,&mem_p_z);

	create_location(nsrc,1,2,1,&src_ori); location_read(src_ori,src_path);
	create_location(nsrc,1,2,1,&src_sta); for(isrc=0;isrc<2*nsrc;isrc++) src_sta.value[isrc] = src_ori.value[isrc] + nabc;
	create_location(nrec,1,2,1,&rec_ori); location_read(rec_ori,rec_path);
	create_location(nrec,1,2,1,&rec_sta); for(irec=0;irec<2*nrec;irec++) rec_sta.value[irec] = rec_ori.value[irec] + nabc;

	matrix_read(0,vp,vp_path); for(n=0;n<mx*mz;n++) vp.value[n] /= 1000.0;
	matrix_ones(rho);

	create_matrix(nrec,1,1,1,nt,dt,&p_signal);

	create_matrix(nt,dt,1,1,1,1,&wavelet); create_matrix(N+1,1,1,1,1,1,&cof); for(i=0;i<N;i++) cof.value[i]=weight(i,2*N);
}

/*release memory space*/
void model_delete()
{
	delete_matrix(&vp); delete_matrix(&padvp); delete_matrix(&rho); delete_matrix(&padrho);

	delete_matrix(&bu); delete_matrix(&lam);

	delete_matrix(&u); delete_matrix(&mem_u);
	delete_matrix(&w); delete_matrix(&mem_w);
	delete_matrix(&p); delete_matrix(&mem_p_x); delete_matrix(&mem_p_z);

	delete_matrix(&at); delete_matrix(&bt);

	delete_location(&src_ori); delete_location(&src_sta); delete_location(&rec_ori); delete_location(&rec_sta);

	delete_matrix(&wavelet); delete_matrix(&cof);

	delete_matrix(&p_signal);
}

/*Weighted coefficiency*/
float weight(int i, int N)
{
	struct matrix a; create_matrix(6,1,1,1,1,1,&a); float value;

	switch(N)
	{
		case 2:  a.value[0]=1;
			break;
		case 4:  a.value[0]=1.12500;a.value[1]=-0.0416660;
			break;
		case 6:  a.value[0]=1.17187;a.value[1]=-0.0651042;a.value[2]=0.0046875;
			break;
		case 8:  a.value[0]=1.19629;a.value[1]=-0.0797526;a.value[2]=0.0095703;a.value[3]=-0.00069754;
			break;
		case 10: a.value[0]=1.21124;a.value[1]=-0.0897217;a.value[2]=0.0138428;a.value[3]=-0.00176566;a.value[4]=0.000118680;
			break;
		case 12: a.value[0]=1.22134;a.value[1]=-0.0969315;a.value[2]=0.0174477;a.value[3]=-0.00296729;a.value[4]=0.000359005;a.value[5]=-0.0000218478;
			break;
	} value = 0.0; value=a.value[i]; return(value);
}

/*Ricker wavelet*/
void ricker_wavelet(int start_point,float frequency,float amplitude,struct matrix sub_wave)
{
        float temp,max; struct matrix rick; create_matrix(nt,1,1,1,1,1,&rick);

        for(it=0;it<nt;it++)
        {
                temp = (PI*frequency*(it*dt-(1/frequency)));
                rick.value[it] = (1-2*(temp*temp))*exp(-temp*temp);
        }
        max = matrix_max(rick);
        for(it=0;it<nt;it++) rick.value[it] *= amplitude/max;

        for(it=start_point;it<nt;it++)
                sub_wave.value[it] = rick.value[it-start_point];

        delete_matrix(&rick);
}

/*update PML coefficient*/
void PML_update(float freq,struct matrix sub_at,struct matrix sub_bt)
{
        matrix_zero(sub_at); matrix_zero(sub_bt); PML_coefficient(2,freq,sub_at,sub_bt);
}

/*Compute coefficient for PML boundary*/
void PML_coefficient(int power,float freq,struct matrix sub_at,struct matrix sub_bt)
{
	float d0,Rcoef=.001,abscissa_in_PML,abscissa_normalized;
	float thick,thi_l,thi_r,thi_t,thi_b,A_max,xval,zval,coef;
	struct matrix Dx,Dz,Ax,Az,ax,bx,az,bz;

	create_matrix(nx,1,1,1,1,1,&Dx); create_matrix(nx,1,1,1,1,1,&Ax);
	create_matrix(nx,1,1,1,1,1,&ax); create_matrix(nx,1,1,1,1,1,&bx);
	create_matrix(nz,1,1,1,1,1,&Dz); create_matrix(nz,1,1,1,1,1,&Az);
	create_matrix(nz,1,1,1,1,1,&az); create_matrix(nz,1,1,1,1,1,&bz);

	thick = dx * npd/2; thi_l = thick; thi_r = (nx+npd/2)*dx; thi_t = thick; thi_b = (nz+npd/2)*dz;

	d0 = (-1.) * (float)(power+1) * matrix_max(vp) * (float)log((double)Rcoef) / (2. * thick);
	A_max = 3.1415*freq;

	for(i=N/2;i<nx-N/2;i++)
	{
		xval = dx * (i-N/2);
		abscissa_in_PML = thi_l - xval;
		if(abscissa_in_PML>=0)
		{
			abscissa_normalized = abscissa_in_PML / thick;
			Dx.value[i] = d0*pow(abscissa_normalized,power);
			Ax.value[i] = A_max * (1.-abscissa_normalized);
		}
		abscissa_in_PML = xval - thi_r;
		if(abscissa_in_PML>=0)
		{
			abscissa_normalized = abscissa_in_PML / thick;
			Dx.value[i] = d0*pow(abscissa_normalized,power);
			Ax.value[i] = A_max * (1.-abscissa_normalized);
		}

		if(Ax.value[i]<0) Ax.value[i] = 0; bx.value[i] = exp(- (Dx.value[i] + Ax.value[i]) * dt);
		if(fabs(Dx.value[i])>1e-6) ax.value[i] = (Dx.value[i]*(bx.value[i]-1.)) / (Dx.value[i]+Ax.value[i]);
	}

	matrix_zero(sub_at); matrix_ones(sub_bt);

	for(k=nabc-1;k>=0;k--)
        {
                coef=ax.value[k];
                for(i=0;i<nx;i++) {sub_at.value[(i)*nz+k] = coef; sub_at.value[(i)*nz+nz-1-k] = coef;}
                for(j=0;j<nz;j++) {sub_at.value[(k)*nz+j] = coef; sub_at.value[(nx-1-k)*nz+j] = coef;}
                coef=bx.value[k];
                for(i=0;i<nx;i++) {sub_bt.value[(i)*nz+k] = coef; sub_bt.value[(i)*nz+nz-1-k] = coef;}
                for(j=0;j<nz;j++) {sub_bt.value[(k)*nz+j] = coef; sub_bt.value[(nx-1-k)*nz+j] = coef;}
        }

	delete_matrix(&Dx); delete_matrix(&Ax); delete_matrix(&ax); delete_matrix(&bx);
	delete_matrix(&Dz); delete_matrix(&Az); delete_matrix(&az); delete_matrix(&bz);
}

/*Compute the physical property of the model*/
void field_transform(struct matrix sub_vp,struct matrix sub_rho)
{
	matrix_zero(padvp); padding(nabc,vp,padvp);
	matrix_zero(padrho);padding(nabc,rho,padrho);

	for(ix=0;ix<nx;ix++) for(iz=0;iz<nz;iz++)
	{
		bu.value[ix*nz+iz] = 1 / padrho.value[ix*nz+iz];
		lam.value[ix*nz+iz] = pow(padvp.value[ix*nz+iz],2) * padrho.value[ix*nz+iz];
	}
}

/*solving first-order velocity-stress wave equation*/
void wave_eq_for_signal(char dir[],int isrc,struct matrix sub_vp,struct matrix sub_rho,struct matrix sub_wave,struct matrix seism)
{
	float u_x,w_z,p_x,p_z;int flag;

	matrix_zero(seism);
	matrix_zero(u); matrix_zero(mem_u);
	matrix_zero(w); matrix_zero(mem_w);
	matrix_zero(p); matrix_zero(mem_p_x); matrix_zero(mem_p_z);

	if(strcmp(dir,"forward")==0) flag=0;
	if(strcmp(dir,"backward")==0) flag=1;

	field_transform(sub_vp,sub_rho);

	for(it=0;it<nt;it++)
	{
		p.value[(src_sta.value[isrc*2+0])*nz+(src_sta.value[isrc*2+1])] += sub_wave.value[it]*dt;

		for(ix=N;ix<nx-N;ix++) for(iz=N;iz<nz-N;iz++)
		{
			p_x=0.0; p_z=0.0;
			for(k=0;k<N;k++)
			{
				p_x += cof.value[k]*(p.value[(ix+k+1)*nz+iz]-p.value[(ix-k)*nz+iz])/dx;
				p_z += cof.value[k]*(p.value[(ix)*nz+iz+1+k]-p.value[(ix)*nz+iz-k])/dz;
			}

			mem_p_x.value[ix*nz+iz] = bt.value[ix*nz+iz]*mem_p_x.value[ix*nz+iz] + at.value[ix*nz+iz]*p_x;
			mem_p_z.value[ix*nz+iz] = bt.value[ix*nz+iz]*mem_p_z.value[ix*nz+iz] + at.value[ix*nz+iz]*p_z;

			u.value[ix*nz+iz] += (flag*(-2)+1)*dt*bu.value[ix*nz+iz]*(p_x+mem_p_x.value[ix*nz+iz]);
			w.value[ix*nz+iz] += (flag*(-2)+1)*dt*bu.value[ix*nz+iz]*(p_z+mem_p_z.value[ix*nz+iz]);
		}

		for(ix=N;ix<nx-N;ix++) for(iz=N;iz<nz-N;iz++)
		{
			u_x=0.0; w_z=0.0;
			for(k=0;k<N;k++)
			{
				u_x += cof.value[k]*(u.value[(ix+k)*nz+iz]-u.value[(ix-1-k)*nz+iz])/dx;
				w_z += cof.value[k]*(w.value[(ix)*nz+iz+k]-w.value[(ix)*nz+iz-1-k])/dz;
			}
			mem_u.value[ix*nz+iz] = bt.value[ix*nz+iz]*mem_u.value[ix*nz+iz] + at.value[ix*nz+iz]*u_x;
			mem_w.value[ix*nz+iz] = bt.value[ix*nz+iz]*mem_w.value[ix*nz+iz] + at.value[ix*nz+iz]*w_z;

			p.value[ix*nz+iz] += (flag*(-2)+1)*dt*lam.value[ix*nz+iz]*(u_x+mem_u.value[ix*nz+iz]+w_z+mem_w.value[ix*nz+iz]);
		}

                for(irec=0;irec<nrec;irec++)
                        if(!(src_sta.value[isrc*2+0] == rec_sta.value[irec*2+0] && src_sta.value[isrc*2+1] == rec_sta.value[irec*2+1]))
                                {seism.value[irec*nt+it] = p.value[(rec_sta.value[irec*2+0])*nz+(rec_sta.value[irec*2+1])];}
	}
}


