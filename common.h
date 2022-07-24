#ifndef __COMMON_H
#define __COMMON_H
#define _USE_MATH_DEFINES

#include <math.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))           // Модуль
#define GRAD_TO_RAD(x) ( (x) * M_PI / 180.)     // Перевод градусы в радианы
#define RAD_TO_GRAD(x) ( (x) * 180. / M_PI)     // Перевод радианы в градусы

#endif //__COMMON_H
