#ifndef __PLANE_CURVE_APPROXIMATION_H
#define __PLANE_CURVE_APPROXIMATION_H

// Аппроксимация набора точек окружностью
float Circle_approximation(float *points, unsigned short points_num, float tolerance, float *x_centr, float *y_centr);

// Аппроксимация набора точек прямой (параметрический вид)
extern float Line_approximation(float *points, unsigned short points_num, float tolerance, float *x_point, float *y_point, float *x_vector, float *y_vector);

#endif //__PLANE_CURVE_APPROXIMATION_H
