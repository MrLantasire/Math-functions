#ifndef __KINEMATICS_H
#define __KINEMATICS_H

#include "vector.h"

// Определение пространственных углов между векторами вокруг заданных осей
extern unsigned char Define_spatial_angles_of_vectors(  vector32_3D_t initial_vec, vector32_3D_t final_vec, vector32_3D_t first_axis, 
                                                        vector32_3D_t second_axis, float out_first_angle[2], float out_second_angle[2]);

#endif //__KINEMATICS_H
