#ifndef __ACESCC_TO_ACES_H_INCLUDED__
#define __ACESCC_TO_ACES_H_INCLUDED__

//
// ACES Color Space Conversion - ACEScc to ACES
//
// converts ACEScc (AP1 w/ ACESlog encoding) to 
//          ACES2065-1 (AP0 w/ linear encoding)
//

//import "ACESlib.Transform_Common";


__DEVICE__ inline float ACEScc_to_lin( float in)
{
  if (in < -0.3013698630f) // (9.72-15)/17.52
    return (_powf( 2.0f, in * 17.52f - 9.72f) - _powf( 2.0f, -16.0f)) * 2.0f;
  else if ( in < (_log2f(HALF_MAX) + 9.72f) / 17.52f )
    return _powf( 2.0f, in * 17.52f - 9.72f);
  else // (in >= (_log2f(HALF_MAX)+9.72f)/17.52f)
    return HALF_MAX;
}

__DEVICE__ inline float3 ACEScc_to_ACES( float3 ACEScc)
{
    float3 lin_AP1;
    lin_AP1.x = ACEScc_to_lin( ACEScc.x);
    lin_AP1.y = ACEScc_to_lin( ACEScc.y);
    lin_AP1.z = ACEScc_to_lin( ACEScc.z);

    float3 ACES = mult_f3_f44( lin_AP1, AP1_2_AP0_MAT);
	return ACES;  
}

#endif