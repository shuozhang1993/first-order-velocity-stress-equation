#ifndef __GENE_H__
#define __GENE_H__
#include "gene.h"
#endif

#ifndef __TOOL_H__
#define __TOOL_H__
#include "tool.h"
#endif

#ifndef __AMF_H__
#define __AMF_H__

/*compute adjoint source from adaptive matching filter*/
void amf_adj_source(struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_p_ajs);

/*data misfit by adaptive matching filter*/
float misfit_amf(struct matrix sub_p_syn,struct matrix sub_p_obs);

/*compute the traveltime in vector based on amf*/
/*traveltime difference has same dimensionality as seismic recording*/
void amf_traveltime_vector(float dt,struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_d_trav_time);

/*compute the traveltime in scalor based on amf*/
/*traveltime difference is a scalor for each trace*/
void amf_traveltime_scalar(float dt,struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_d_trav_time);

/*compute adjoint source for conventional traveltime measurement*/
void amf_scalar_traveltime_adj_source(float dt,struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_p_ajs);

/*compute adjoint source for traveltime from amf with panlty function*/
void amf_vector_traveltime_adj_source(float dt,struct matrix sub_p_syn,struct matrix sub_p_obs,struct matrix sub_p_ajs);

#endif
