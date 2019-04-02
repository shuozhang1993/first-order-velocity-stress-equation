#include "gene.h"

/*Transform float matrix to absolute value*/
void matrix_fabs(struct matrix mat_in,struct matrix mat_out)
{
	for(n=0;n<mat_in.x*mat_in.y*mat_in.z;n++) mat_out.value[n]=mat_in.value[n];
}

/*Get norm of the matrix*/
float matrix_norm(struct matrix mat)
{
	float norm=0.0; for(i=0;i<mat.x*mat.y*mat.z;i++) norm += fabs(mat.value[i]); return(norm);
}

/*Pick up maximum value*/
float matrix_max(struct matrix mat)
{
	float maxi=0; maxi = mat.value[0]; for(n=0;n<mat.x*mat.y*mat.z;n++) if(maxi<mat.value[n]) maxi=mat.value[n]; return(maxi);	
}

/*Pick up maximum absolute value*/
float matrix_abs_max(struct matrix mat)
{
	float maxi=fabs(mat.value[0]); for(n=0;n<mat.x*mat.y*mat.z;n++) if(maxi<fabs(mat.value[n])) maxi=fabs(mat.value[n]); return(maxi);	
}

/*Pick up minimum value*/
float matrix_min(struct matrix mat)
{
	float mini=0; mini = mat.value[0]; for(n=0;n<mat.x*mat.y*mat.z;n++) if(mini>mat.value[n]) mini=mat.value[n]; return(mini);	
}

/*Get difference of two matrices*/
void matrix_diff(struct matrix mat1,struct matrix mat2,struct matrix mat3)
{
	for(n=0;n<mat1.x*mat1.y*mat1.z;n++) mat3.value[n] = mat1.value[n] - mat2.value[n];
}

/*Get summation of two 1D martices*/
void matrix_sum(struct matrix mat1,struct matrix mat2,struct matrix mat3)
{
	for(n=0;n<mat1.x*mat1.y*mat1.z;n++) mat3.value[n] = mat1.value[n] + mat2.value[n];
}

/*Reset matrix*/
void matrix_zero(struct matrix mat)
{
	for(n=0;n<mat.x*mat.y*mat.z;n++) mat.value[n]=0.0;
}

/*Identity matrix*/
void matrix_ones(struct matrix mat)
{
	for(n=0;n<mat.x*mat.y*mat.z;n++) mat.value[n] = 1.0;
}

/*----------Data structure I/O----------*/

/*Read binary file*/
void matrix_read(int delay,struct matrix mat,char filename[])
{
	int count; fp = fopen(filename,"rb");

	for(i=delay;i<mat.x-delay;i++)
		for(j=delay;j<mat.y-delay;j++)
			for(k=delay;k<mat.z-delay;k++)
				count = fread(&mat.value[j*mat.x*mat.z+i*mat.z+k],sizeof(float),1,fp);
	fclose(fp); fp = NULL;
}


/*Write binary file*/
void matrix_write(int delay,struct matrix mat,char filename[])
{
	int count; fp = fopen(filename,"wb");

	for(i=delay;i<mat.x-delay;i++)
		for(j=delay;j<mat.y-delay;j++)
			for(k=delay;k<mat.z-delay;k++) 
				count = fwrite(&mat.value[j*mat.x*mat.z+i*mat.z+k],sizeof(float),1,fp);

	fclose(fp); fp = NULL;
}

/*Read 2D data from text file*/
void location_read(struct location station,char filename[])
{
	int count; fp = fopen(filename,"rt");

	for(i=0;i<station.x;i++)
	{
		for(j=0;j<station.y;j++)
			count = fscanf(fp,"%d",&station.value[i*station.y+j]);
		count = fscanf(fp,"\n");
	}

	fclose(fp); fp = NULL;
}

/*----------Data structure creation and deletion----------*/

/*Create matrix*/
void create_matrix(int a,float da,int b,float db,int c,float dc,struct matrix *mat)
{
	mat->x=a; mat->dx=da; mat->y=b; mat->dy=db; mat->z=c; mat->dz=dc;
	mat->value=(float*)malloc(sizeof(float)*a*b*c);  for(n=0;n<a*b*c;n++) mat->value[n] = 0.0;
}

/*Delete matrix*/
void delete_matrix(struct matrix *mat)
{
	free((void*)mat->value); mat->x=0; mat->dx=0.0; mat->y=0; mat->dy=0.0;  mat->z=0; mat->dz=0.0;
}

/*Create complex matrix*/
void create_cpxmat(int a,float da,int b,float db,int c,float dc,struct cpxmat *cpx)
{
	cpx->x=a; cpx->dx=da; cpx->y=b; cpx->dy=db; cpx->z=c; cpx->dz=dc;
	cpx->real=(float*)malloc(sizeof(float)*a*b*c); cpx->imag=(float*)malloc(sizeof(float)*a*b*c);
	for(n=0;n<a*b*c;n++) {cpx->real[n] = 0.0; cpx->imag[n] = 0.0;}
}

/*Delete complex matrix*/
void delete_cpxmat(struct cpxmat *cpx)
{
	free((void*)cpx->real); free((void*)cpx->imag);
	cpx->x=0; cpx->dx=0.0; cpx->y=0; cpx->dy=0.0; cpx->z=0; cpx->dz=0.0;
}

/*Create location*/
void create_location(int a,float da,int b,float db,struct location *sta)
{
	sta->x=a; sta->dx=da; sta->y=b; sta->dy=db;
	sta->value=(int*)malloc(sizeof(int)*a*b); for(n=0;n<a*b;n++) sta->value[n] = 0.0;
}

/*Delete location*/
void delete_location(struct location *sta)
{
	free((void*)sta->value); sta->x=0; sta->dx=0.0; sta->y=0; sta->dy=0.0;
}

/*Create filter*/
void create_filter(int a,float b,struct filter *fil)
{
	fil->length=a; fil->interval=b; fil->value=(float*)malloc(sizeof(float)*a); for(n=0;n<a;n++) fil->value[n]=0.0;
}

/*Delete filter*/
void delete_filter(struct filter *fil)
{
	free((void*)fil->value); fil->length = 0; fil->interval = 0.0;
}

/*----------Trace opeartion----------*/

/*Reverse order of 1D trace*/
void reverse_trace(struct matrix ip,struct matrix op)
{
	int il,len=ip.x*ip.z*ip.y;

	for(il=0;il<len;il++) op.value[il] = ip.value[len-1-il];
}

/*----------Matrix operation----------*/

/*Transpose*/
void transpose(struct matrix mat1,struct matrix mat2)
{
	for(i=0;i<mat2.x;i++) for(j=0;j<mat2.z;j++) mat2.value[i*mat2.z+j] = mat1.value[j*mat1.z+i];	
}

/*Transite value*/
void transition(struct matrix mat1,struct matrix mat2)
{
	for(n=0;n<mat1.x*mat1.y*mat1.z;n++) mat2.value[n] = mat1.value[n];
}

/*Extend boundary area*/
void padding(int nabc,struct matrix vel,struct matrix padvel)
{
	int mx=vel.x,mz=vel.z,nx=padvel.x,nz=padvel.z;
	for(i=0;i<mx;i++) for(j=0;j<mz;j++) padvel.value[(i+nabc)*nz+(j+nabc)] = vel.value[(i)*mz+(j)];

	for(i=0;i<nx;i++) for(j=0;j<nabc;j++) padvel.value[(i)*nz+(j)] = padvel.value[(i)*nz+(nabc)];
	for(i=0;i<nx;i++) for(j=0;j<nabc;j++) padvel.value[(i)*nz+(j+mz+nabc)] = padvel.value[(i)*nz+(mz+nabc-1)];
	for(i=0;i<nabc;i++) for(j=0;j<nz;j++) padvel.value[(i)*nz+(j)] = padvel.value[(nabc)*nz+(j)];
	for(i=0;i<nabc;i++) for(j=0;j<nz;j++) padvel.value[(i+mx+nabc)*nz+(j)] = padvel.value[(mx+nabc-1)*nz+(j)];
}

/*Cross product between two matrices*/
void crossprod(struct matrix mat1,struct matrix mat2,struct matrix mat3)
{
	float sum=0.0;; 
	for(i=0;i<mat1.x;i++) for(j=0;j<mat2.z;j++)
	{
		sum = 0.0;
		for(k=0;k<mat2.x;k++) 
			sum = mat1.value[i*mat1.z+k]*mat2.value[k*mat2.z+j];
		mat3.value[i*mat3.z+j] = sum;
	}
}

/*Inner product between two matrices*/
void innerprod(struct matrix mat1,struct matrix mat2,struct matrix mat3)
{
	for(n=0;n<mat1.x*mat1.y*mat1.z;n++) mat3.value[n] = mat1.value[n]*mat2.value[n];
}

/*Product between scalar and matrix*/
void constprod(float constant,struct matrix mat)
{
	for(n=0;n<mat.x*mat.y*mat.z;n++) mat.value[n] *= constant;
}

/*Dot product between two matrices*/
float dotprod(struct matrix mat1,struct matrix mat2)
{
	int n; float sum=0.0; for(n=0;n<mat1.x*mat1.y*mat1.z;n++) sum += mat1.value[n] * mat2.value[n]; return(sum);
}

/*Normalize matrix 1 with matrix 2*/
void normal_mat(struct matrix mat1,struct matrix mat2)
{
	float max1,max2;	max1 = matrix_abs_max(mat1); max2 = matrix_abs_max(mat2); 
	for(n=0;n<mat1.x*mat1.y*mat1.z;n++) mat1.value[n] = (max2/max1) * mat1.value[n];
}
