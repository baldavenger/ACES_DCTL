#ifndef __IDT_SONY_SLOG1_SGAMUT_10_H_INCLUDED__
#define __IDT_SONY_SLOG1_SGAMUT_10_H_INCLUDED__

//
// IDT for Sony Cameras - 10 bits
// Provided by Sony Electronics Corp.
//


__DEVICE__ inline float3 IDT_Sony_SLog1_SGamut_10( float3 In)
{

mat3 SGAMUT_TO_ACES_MTX = { { 0.754338638f,  0.021198141f, -0.009756991f },
							{ 0.133697046f,  1.005410934f,  0.004508563f },
							{ 0.111968437f, -0.026610548f,  1.005253201f } };
							
const float B = 64.0f;
const float AB = 90.0f;
const float W = 940.0f;

// Prepare input values based on application bit depth handling
float3 SLog;
SLog.x = In.x * 1023.0f;
SLog.y = In.y * 1023.0f;
SLog.z = In.z * 1023.0f;

// 10-bit Sony S-log to linear S-gamut
float3 lin;
lin.x = SLog1_to_lin( SLog.x, B, AB, W);
lin.y = SLog1_to_lin( SLog.y, B, AB, W);
lin.z = SLog1_to_lin( SLog.z, B, AB, W);

// S-Gamut to ACES matrix
float3 aces = mult_f3_f33( lin, SGAMUT_TO_ACES_MTX);

return aces;

}

#endif