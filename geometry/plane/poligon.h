#ifndef __POLIGON_H
#define __POLIGON_H

#include <stdbool.h>

// Функция определения нахождения точки внутри простого многоугольника
extern bool Is_point_in_poligon(float coordinate_x, float coordinate_y, float *poligon, unsigned char amount_of_vertices, float tolerance);

#endif //__POLIGON_H
