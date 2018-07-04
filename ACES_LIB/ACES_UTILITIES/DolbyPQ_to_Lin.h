#ifndef __DOLBYPQ_TO_LIN_H_INCLUDED__
#define __DOLBYPQ_TO_LIN_H_INCLUDED__

// 
// Generic transform from SMPTE ST2084 to linear
// 

//import "ACESlib.Utilities_Color";


__DEVICE__ inline float3 DolbyPQ_to_Lin( float3 In)
{
  float3 out;
  out.x = ST2084_2_Y( In.x );
  out.y = ST2084_2_Y( In.y );
  out.z = ST2084_2_Y( In.z );
  return out;
}

#endif