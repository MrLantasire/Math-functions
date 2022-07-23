#ifndef __EQUATION_H
#define __EQUATION_H

#include "complex.h"

// Решение квадратного уравнения
extern float Solve_quadratic_equation(float a, float b, float c, complex32_t out[2]);

// Решение кубического уравнения
extern float Solve_cubic_equation(float a, float b, float c, float d, complex32_t out[3]);

#endif //__EQUATION_H
