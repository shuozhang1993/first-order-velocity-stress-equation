#ifndef __GENE_H__
#define __GENE_H__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#ifndef __FFTW_H__
#define __FFTW_H__
//#include <fftw3.h>
#endif

struct matrix
{
	int x; float dx;
	int y; float dy;
	int z; float dz;
	float *value;
};

struct cpxmat
{
	int x; float dx;
	int y; float dy;
	int z; float dz;
	float *real;
	float *imag;
};

struct location
{
	int x; float dx;
	int y; float dy;
	int *value;
};

struct filter
{
	int length;
	float interval;
	float *value;
};

#define string_length 100
#define PI 3.1415926

FILE *fp; char filename[string_length];
int i,j,k,n;

# if 0
typedef  unsigned char      boolean;     /* Boolean value type. */
typedef  unsigned long int  uint32;      /* Unsigned 32 bit value */
typedef  unsigned short     uint16;      /* Unsigned 16 bit value */
typedef  unsigned char      uint8;       /* Unsigned 8  bit value */
typedef  signed long int    int32;       /* Signed 32 bit value */
typedef  signed short       int16;       /* Signed 16 bit value */
typedef  signed char        int8;        /* Signed 8  bit value */
# endif

/*Compare and pick up maximal value*/
#define  MAX(x,y) (((x) > (y)) ? (x) : (y))
/*Compare and pick up minimal value*/
#define  MIN(x,y) (((x) < (y)) ? (x) : (y))
/*Return element number of array*/
#define ARR_SIZE(a)  (sizeof((a))/sizeof((a[0])))
/*Return element number of matrix*/
#define MAT_SIZE(a)  (sizeof((a))/sizeof((a[0][0])))

/*Transform float matrix to absolute value*/
void matrix_fabs(struct matrix mat_int,struct matrix mat_out);

/*Get norm of the matrix*/
float matrix_norm(struct matrix mat);

/*Pick up maximum value*/
float matrix_max(struct matrix mat);

/*Pick up maximum absolute value*/
float matrix_abs_max(struct matrix mat);

/*Pick up minimum value*/
float matrix_min(struct matrix mat);

/*Get difference of two matrices*/
void matrix_diff(struct matrix mat1,struct matrix mat2,struct matrix mat3);

/*Get summation of two 1D martices*/
void matrix_sum(struct matrix mat1,struct matrix mat2,struct matrix mat3);

/*Reset matrix*/
void matrix_zero(struct matrix mat);

/*Identity matrix*/
void matrix_ones(struct matrix mat);

/*Data structure I/O*/

/*Read binary file*/
void matrix_read(int delay,struct matrix mat,char filename[]);

/*Write binary file*/
void matrix_write(int delay,struct matrix mat,char filename[]);

/*Read 2D data from text file*/
void location_read(struct location mat,char filename[]);

/*Data structure creation and deletion*/

/*Create matrix*/
void create_matrix(int a,float da,int b,float db,int c,float dc,struct matrix *mat);

/*delete matrix*/
void delete_matrix(struct matrix *mat);

/*Create complex matrix*/
void create_cpxmat(int a,float da,int b,float db,int c,float dc,struct cpxmat *cpx);

/*Delete complex matrix*/
void delete_cpxmat(struct cpxmat *cpx);

/*Create location*/
void create_location(int a,float da,int b,float db,struct location *sta);

/*Delete location*/
void delete_location(struct location *sta);

/*Create filter*/
void create_filter(int a,float b,struct filter *fil);

/*Delete filter*/
void delete_filter(struct filter *fil);

/*Trace operation*/

/*Reverse order of trace*/
void reverse_trace(struct matrix ip,struct matrix op);

/*Matrix operation*/

/*Transpose*/
void transpose(struct matrix mat1,struct matrix mat2);

/*Transite value*/
void transition(struct matrix mat1,struct matrix mat2);

/*Extend boundary area*/
void padding(int bond,struct matrix in,struct matrix out);

/*Cross product between two matrices*/
void crossprod(struct matrix mat1,struct matrix mat2,struct matrix mat3);

/*Inner product between two matrices*/
void innerprod(struct matrix mat1,struct matrix mat2,struct matrix mat3);

/*Product between scalar and matrix*/
void constprod(float constant,struct matrix mat);

/*Dot product between two matrices*/
float dotprod(struct matrix mat1,struct matrix mat2);

/*Normalize matrix 1 with matrix 2*/
void normal_mat(struct matrix mat1,struct matrix mat2);

# endif
