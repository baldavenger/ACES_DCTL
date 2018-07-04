#ifndef __IDT_SONY_SLOG2_SGAMUT_DAYLIGHT_12_H_INCLUDED__
#define __IDT_SONY_SLOG2_SGAMUT_DAYLIGHT_12_H_INCLUDED__

//
// IDT for Sony Cameras - 12 bits - Daylight (5500K)
// Provided by Sony Electronics Corp.
//


__DEVICE__ inline float3 IDT_Sony_SLog2_SGamut_Daylight_12( float3 In)
{

mat3 SGAMUT_DAYLIGHT_TO_ACES_MTX = { { 0.8764457030f,  0.0774075345f,  0.0573564351f},
  									{ 0.0145411681f,  0.9529571767f, -0.1151066335f},
  									{ 0.1090131290f, -0.0303647111f,  1.0577501984f} };
							
const float B = 256.0f;
const float AB = 360.0f;
const float W = 3760.0f;

// Prepare input values based on application bit depth handling
float3 SLog;
SLog.x = In.x * 4095.0f;
SLog.y = In.y * 4095.0f;
SLog.z = In.z * 4095.0f;

// 12-bit Sony S-log to linear S-gamut
float3 lin;
lin.x = SLog2_to_lin( SLog.x, B, AB, W);
lin.y = SLog2_to_lin( SLog.y, B, AB, W);
lin.z = SLog2_to_lin( SLog.z, B, AB, W);

// S-Gamut to ACES matrix
float3 aces = mult_f3_f33( lin, SGAMUT_DAYLIGHT_TO_ACES_MTX);

return aces;

}

#endif