#ifndef __POLYGON_H
#define __POLYGON_H

#include <stdbool.h>

// Функция определения нахождения точки внутри простого многоугольника
extern bool Is_point_in_polygon(float coordinate_x, float coordinate_y, float *polygon, unsigned char amount_of_vertices, float tolerance);

// Функция вычисления площади простого многоугольника
extern float Calculate_polygon_area(float *polygon, unsigned char amount_of_vertices);

#endif //__POLYGON_H
