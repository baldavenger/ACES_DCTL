#ifndef __ACES_PROXY12_TO_ACES_H_INCLUDED__
#define __ACES_PROXY12_TO_ACES_H_INCLUDED__

//
// ACES Color Space Conversion - ACESproxy (12-bit) to ACES
//
// converts ACESproxy (AP1 w/ ACESproxy encoding) to 
//          ACES2065-1 (AP0 w/ linear encoding)
//

//import "ACESlib.Transform_Common";


__DEVICE__ inline float ACESproxy12_to_lin( int in)
{
float StepsPerStop = 200.0f;
float MidCVoffset = 1700.0f;
int CVmin = 256;
int CVmax = 3760;

return _powf( 2.0f, ( in - MidCVoffset)/StepsPerStop - 2.5f);
}


__DEVICE__ inline float ACESproxy12_to_ACES( int In[3])
{
    
int ACESproxy[3];
ACESproxy[0] = In.x * 4095.0f;
ACESproxy[1] = In.y * 4095.0f;
ACESproxy[2] = In.z * 4095.0f;

float3 lin_AP1;
lin_AP1.x = ACESproxy12_to_lin( ACESproxy[0]);
lin_AP1.y = ACESproxy12_to_lin( ACESproxy[1]);
lin_AP1.z = ACESproxy12_to_lin( ACESproxy[2]);

float3 ACES = mult_f3_f44( lin_AP1, AP1_2_AP0_MAT);
return ACES;

}

#endif