#include "gene.h"
#include "model.h"

#ifndef ___MPI_H__
#define ___MPI_H__
#include <mpi.h>
#endif

const int field_type=1;/*0 for elastic, 1 for acoustic*/
const int taper_type=1;/*0 for taper, 1 for no taper*/
const int image_type=1;/*0 for regular, 1 for separate*/
const int bound_type=1;/*0 for sponge, 1 for PML*/
const int separation=1;/*0 for no P/S separation, 1 for P/S separation*/
const int diff_order=8;/*order of differential*/
const int abs_width=100;/*width of absorption boundary*/

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
	int i,count; FILE *fp=fopen(file_path,"rt"); if(rank==0 && fp==NULL) {printf("\nfail to read the coefficient text file.\n");exit(1);}

	count=fscanf(fp,"%d\n",&hori_grid);   count=fscanf(fp,"%d\n",&vert_grid);
	count=fscanf(fp,"%f\n",&hori_inter);  count=fscanf(fp,"%f\n",&vert_inter);   count=fscanf(fp,"%f\n",&lenth_phy);  count=fscanf(fp,"%f\n",&depth_phy);
	count=fscanf(fp,"%d\n",&time_step);   count=fscanf(fp,"%f\n",&time_interval);count=fscanf(fp,"%d\n",&shot_number);count=fscanf(fp,"%d\n",&rece_number);
	count=fscanf(fp,"%d\n",&freq_step);   count=fscanf(fp,"%f\n",&freq_int);     count=fscanf(fp,"%f\n",&freq_low);   count=fscanf(fp,"%f\n",&freq_cen);
	count=fscanf(fp,"%d\n",&itera_number);

	count=fscanf(fp,"%s\n",src_path);     count=fscanf(fp,"%s\n",rec_path);
	count=fscanf(fp,"%s\n",vp_path);      count=fscanf(fp,"%s\n",vs_path);       count=fscanf(fp,"\n%s\n",rho_path);

	if(rank==0 && fp!=NULL) printf("\nSuccessful read the coefficient from text file.\n");

	fclose(fp); fp = NULL;

	nx= hori_grid; nz = vert_grid; N = diff_order; npd = abs_width;
	x = hori_grid+diff_order+abs_width; z = vert_grid+diff_order+abs_width; boundary = N/2+npd/2;
	dt= time_interval; dh = vert_inter;
}

/*allocate memory and set model value*/
void model_arrang()
{
	int n;

	create_matrix(hori_grid,1,vert_grid,&Vp); create_matrix(hori_grid,1,vert_grid,&Rho); create_matrix(x,1,z,&vp); create_matrix(x,1,z,&rho);

	create_matrix(x,1,z,&BU); create_matrix(x,1,z,&LAM); create_matrix(x,1,z,&u); create_matrix(x,1,z,&w); create_matrix(x,1,z,&p);

	create_matrix(x,1,z,&phi_u); create_matrix(x,1,z,&phi_w); create_matrix(x,1,z,&phi_p_x); create_matrix(x,1,z,&phi_p_z);

	create_location(shot_number,2,&src_sta); location_read(src_sta,src_path); for(n=0;n<src_sta.x*src_sta.y;n++) src_sta.value[n] += boundary;
	create_location(rece_number,2,&rec_sta); location_read(rec_sta,rec_path); for(n=0;n<rec_sta.x*rec_sta.y;n++) rec_sta.value[n] += boundary;

	matrix_read(0,Vp,vp_path); /*matrix_read(0,Rho,rho_path);*/ matrix_ones(Rho); extend(Vp,vp); extend(Rho,rho); field_transform();

	vpmax = matrix_max(Vp); if(rank==0 && stability_condition()==1) exit(1);

	create_matrix(diff_order/2,1,1,&cof); for(n=0;n<diff_order/2;n++) cof.value[n] = weight(n,N); create_matrix(time_step,1,1,&wavelet);
}

/*zero propagation parameter*/
void zero_propagation()
{
	matrix_zero(u); matrix_zero(w); matrix_zero(p); matrix_zero(phi_u); matrix_zero(phi_w); matrix_zero(phi_p_x); matrix_zero(phi_p_z);
}

/*free memory space for propagation parameter*/
void free_propagation()
{
	delete_matrix(&u); delete_matrix(&w); delete_matrix(&p); delete_matrix(&phi_u); delete_matrix(&phi_w); delete_matrix(&phi_p_x); delete_matrix(&phi_p_z);
}

/*free memory for model value*/
void free_model()
{
	delete_matrix(&Vp); delete_matrix(&Rho); delete_matrix(&vp); delete_matrix(&rho); delete_matrix(&BU); delete_matrix(&LAM);
	delete_matrix(&cof); delete_matrix(&wavelet);
	delete_location(&src_sta); delete_location(&rec_sta);
}

/*PML setting*/
void PML_coef(float freq)
{

	create_matrix(x,1,1,&ax); create_matrix(x,1,1,&bx); create_matrix(z,1,1,&az); create_matrix(z,1,1,&bz); PML_coefficient(2,freq);

	MPI_Bcast(ax.value,x,MPI_FLOAT,0,MPI_COMM_WORLD); MPI_Bcast(az.value,z,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(bx.value,x,MPI_FLOAT,0,MPI_COMM_WORLD); MPI_Bcast(bz.value,z,MPI_FLOAT,0,MPI_COMM_WORLD);
}

/*zero PML parameter*/
void PML_zero()
{
	delete_matrix(&ax); delete_matrix(&az); delete_matrix(&bx); delete_matrix(&bz);
}

/*####################source wavelet####################*/

/*Source function*/
float ricker(float freq, float time)
{
	float s,coef,index;	
	coef=1.0-2.0*(3.14*3.14)*(freq*freq)*((time-1.0/freq)*(time-1.0/freq));
	index=(3.14*3.14)*(freq*freq)*((time-1.0/freq)*(time-1.0/freq));
	s=coef*exp(-1.0*index);
	return(s);
}

/*Source wavelet*/
void sou_let(float freq)
{
	int i; for(i=0;i<time_step;i++) wavelet.value[i] = ricker(freq,i*dt);
}

/*Weighted coefficiency*/
float weight(int i, int N)
{
	float a[6],out;
	int j;
	for(j=0;j<6;j++)
		a[j]=0.0;
	switch(N)
	{
		case 2:  a[0]=1;
			break;
		case 4:  a[0]=1.125;a[1]=-0.041666;
			break;
		case 6:  a[0]=1.17187;a[1]=-0.0651042;a[2]=0.0046875;
			break;
		case 8:  a[0]=1.19629;a[1]=-0.0797526;a[2]=0.0095703;a[3]=-0.00069754;
			break;
		case 10: a[0]=1.21124;a[1]=-0.0897217;a[2]=0.0138428;a[3]=-0.00176566;a[4]=0.000118680;
			break;
		case 12: a[0]=1.22134;a[1]=-0.0969315;a[2]=0.0174477;a[3]=-0.00296729;a[4]=0.000359005;a[5]=-0.0000218478;
			break;
	}
	out=a[i];
	return(out);
}

/*Stable adjudgement*/
int stability_condition()
{	
	float res;
	
	res = (dh*1)/(sqrt(diff_order)*vpmax);

	if(dt >= res)
	{
		printf("\nStable condition is %.8fs.\n",res);
		printf("\nThis dt-dh relationship can't stablize computation.\n");
		printf("\nFinite Difference will generate dispersion!\n");
		printf("\nCorrect the time and space interval is necessary.\n");
		printf("\n----------------------Reverse Time Migration finish--------------------------\n");
		return(1);
	}
	else 
	{
		printf("\nStable condition is %.8fs.\n",res);
		printf("\nThe numerical model passes test and can be continue.\n");
		return(0);
	}
}

/*Compute the physical property of the model*/
void field_transform()
{
	int i,j; float a=dt;

	for(i=0;i<x;i++)/*original one*/
	for(j=0;j<z;j++)
	{
		BU.value[i*BU.z+j] = a / rho.value[i*rho.z+j];
		LAM.value[i*LAM.z+j] = a * vp.value[i*vp.z+j] * vp.value[i*vp.z+j] * rho.value[i*rho.z+j];
	}
}

/*Extend boundary area*/
void extend(struct matrix in,struct matrix out)
{
	int i,j;

	/*fulfill central area*/
	for(i=0;i<hori_grid;i++) for(j=0;j<vert_grid;j++) out.value[(boundary+i)*out.z+boundary+j] = in.value[i*in.z+j];

	for(j=0;j<boundary;j++) for(i=boundary;i<hori_grid+boundary;i++)
	{
		out.value[i*out.z+j] = out.value[i*out.z+boundary];/*Upper*/
		out.value[(i)*out.z+vert_grid+boundary+j] = out.value[i*out.z+(vert_grid+boundary-1)];/*Bottom*/
	}

	for(j=0;j<vert_grid+2*boundary;j++) for(i=0;i<boundary;i++)
	{
		out.value[i*out.z+j] = out.value[boundary*out.z+j];/*Left*/
		out.value[(hori_grid+boundary+i)*out.z+j] = out.value[(hori_grid+boundary-1)*out.z+j];/*Right*/
	}
}

/*Compute coefficient for PML boundary*/
void PML_coefficient(int power,float freq)
{
	int i,j;
	float d0,Rcoef=.001,abscissa_in_PML,abscissa_normalized;
	float thick,thi_l,thi_r,thi_t,thi_b,A_max,xval,zval;
	struct matrix Dx,Dz,Ax,Az;

	create_matrix(nx+2*boundary,1,1,&Dx); create_matrix(nz+2*boundary,1,1,&Dz);
	create_matrix(nx+2*boundary,1,1,&Ax); create_matrix(nz+2*boundary,1,1,&Az);

	thick = dh * npd/2; thi_l = thick; thi_r = (nx+npd/2)*dh; thi_t = thick; thi_b = (nz+npd/2)*dh;

	d0 = (-1.) * (float)(power+1) * vpmax * (float)log((double)Rcoef) / (2. * thick);
	A_max = 3.1415*freq;

	for(i=N/2;i<nx+2*boundary-N/2;i++)/*x direction*/
	{
		xval = dh * (i-N/2);
		/*left boundary*/
		abscissa_in_PML = thi_l - xval;
		if(abscissa_in_PML>=0)
		{
			abscissa_normalized = abscissa_in_PML / thick;
			Dx.value[i] = d0*pow(abscissa_normalized,power);
			Ax.value[i] = A_max * (1.-abscissa_normalized);
		}
		/*right boundary*/
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
	for(j=N/2;j<nz+2*boundary-N/2;j++)/*z direction*/
	{
		zval = dh * (j-N/2);
		/*top boundary*/
		abscissa_in_PML = thi_t - zval;
		if(abscissa_in_PML>=0)
		{
			abscissa_normalized = abscissa_in_PML / thick;
			Dz.value[j] = d0*pow(abscissa_normalized,power);
			Az.value[j] = A_max * (1.-abscissa_normalized); 
		}
		/*bottom boundary*/
		abscissa_in_PML = zval - thi_b;
		if(abscissa_in_PML>=0)
		{
			abscissa_normalized = abscissa_in_PML / thick;
			Dz.value[j] = d0*pow(abscissa_normalized,power);
			Az.value[j] = A_max * (1.-abscissa_normalized); 
		}

		if(Az.value[j]<0) Az.value[j] = 0; bz.value[j] = exp(-(Dz.value[j]+Az.value[j])*dt);
		if(fabs(Dz.value[j])>1e-6) az.value[j] = (Dz.value[j]*(bz.value[j]-1.)) / (Dz.value[j]+Az.value[j]);
	}

//	for(i=0;i<(int)x/2;i++) {ax.value[x-1-i]=ax.value[i]; bx.value[x-1-i]=bx.value[i];}
//	for(j=0;j<(int)z/2;j++) {az.value[z-1-j]=az.value[j]; bz.value[z-1-j]=bz.value[j];}

	delete_matrix(&Dx); delete_matrix(&Dz); delete_matrix(&Ax); delete_matrix(&Az);
}

/*Update stress field in acoustic case with PML*/
void PML_stress_acoustic(int flag)
{
	int i,j,k;
	float u_x,w_z;

	/*Stress from velocity*/	
	for(i=N/2;i<x-N/2;i++)
	for(j=N/2;j<z-N/2;j++)
	{
		/*Nth order FD*/
		u_x=0.;w_z=0.;
		for(k=0;k<N/2;k++)
		{
			u_x = u_x + cof.value[k]*(u.value[(i+k)*u.z+j]-u.value[(i-1-k)*u.z+j])/dh;
			w_z = w_z + cof.value[k]*(w.value[(i)*w.z+j+k]-w.value[(i)*w.z+j-1-k])/dh;
		}
		phi_u.value[i*phi_u.z+j] = bx.value[i]*phi_u.value[i*phi_u.z+j] + ax.value[i]*u_x;
		phi_w.value[i*phi_w.z+j] = bz.value[j]*phi_w.value[i*phi_w.z+j] + az.value[j]*w_z;
		p.value[i*p.z+j] += (flag*(-2)+1)*LAM.value[i*LAM.z+j]*(u_x+phi_u.value[i*phi_u.z+j]+w_z+phi_w.value[i*phi_w.z+j]);
	}	
}

/*Update displacement field in acoustic case*/
void PML_velo_acoustic(int flag)
{
	int i,j,k;
	float p_x,p_z;

	for(i=N/2;i<x-N/2;i++)
	for(j=N/2;j<z-N/2;j++)
	{
		p_x=0.0;p_z=0.0;
		for(k=0;k<N/2;k++)
		{
			p_x = p_x + cof.value[k]*(p.value[(i+1+k)*p.z+j]-p.value[(i-k)*p.z+j])/dh;
			p_z = p_z + cof.value[k]*(p.value[(i)*p.z+j+1+k]-p.value[(i)*p.z+j-k])/dh;
//			p_x = p_x + cof.value[k]*(p.value[(i+k)*p.z+j]-p.value[(i-1-k)*p.z+j])/dh;
//			p_z = p_z + cof.value[k]*(p.value[(i)*p.z+j+k]-p.value[(i)*p.z+j-1-k])/dh;
		}
		phi_p_x.value[i*phi_p_x.z+j] = bx.value[i]*phi_p_x.value[i*phi_p_x.z+j] + ax.value[i]*p_x;
		phi_p_z.value[i*phi_p_z.z+j] = bz.value[j]*phi_p_z.value[i*phi_p_z.z+j] + az.value[j]*p_z;
		u.value[i*u.z+j] += (flag*(-2)+1)*(BU.value[i*BU.z+j])*(p_x+phi_p_x.value[i*phi_p_x.z+j]);
		w.value[i*w.z+j] += (flag*(-2)+1)*(BU.value[i*BU.z+j])*(p_z+phi_p_z.value[i*phi_p_z.z+j]);				
	}
}

