#include <malloc.h>
#include "plane_curve_approximation.h"
#include "common.h"
#include "algebra\linear\matrix.h"

/***********************************************************
* Функция нахождения аппроксимирующей окружности, описываемой уравнением вида
* (x - x_center)^2 + (y - y_centr) = r^2 по методу наименьших квадратов
* Входные данные:
* points - Указатель на массив точек с координатами
* хранящимися в двумерном массиве points с размерами points_num * 2 (строка {x,y}).
* points_num - Количество точек в массиве
* Выходные данные:
* x_centr - Координата X центра окружности
* y_centr - Координата Y центра окружности
* r - Радиус окружности
************************************************************/
float Circle_approximation(float *points, unsigned short points_num, float tolerance, float *x_centr, float *y_centr)
{
    float r = 0.0;                              // Радиус окружности
    float *matrix = calloc(12, sizeof(float));  // Матрица системы уравннений

    // Обнуление выходных параметров 
    *x_centr = 0.0;
    *y_centr = 0.0;

    // Обнуление множителей
    for (unsigned short i = 0; i < 12U; i++)
    {
        matrix[i] = 0.0;
    }

    // Вычисление множителей для системы уравнений
    float xx = 0.0;
    float yy = 0.0;
    for (unsigned short i = 0; i < points_num; i++)
    {
        xx = points[i * 2] * points[i * 2];
        yy = points[i * 2 + 1] * points[i * 2 + 1];

        matrix[0] += xx;                                    // x^2
        matrix[1] += points[i * 2] * points[i * 2 + 1];     // xy
        matrix[2] += points[i * 2];                         // x
        matrix[3] -= points[i * 2] * (xx + yy);             // -x(x^2 + y^2)
        matrix[5] += yy;                                    // y^2
        matrix[6] += points[i * 2 + 1];                     // y
        matrix[7] -= points[i * 2 + 1] * (xx + yy);         // -y(x^2 + y^2)
        matrix[11] -= xx + yy;                              // -(x^2 + y^2)
    }
    matrix[4] = matrix[1];                                  // xy
    matrix[8] = matrix[2];                                  // x
    matrix[9] = matrix[6];                                  // y
    matrix[10] = (float) points_num;                        // n

    if (Gaussian_elimination(matrix, 3, 4, tolerance, matrix) == 3U)
    {
        // Решение системы уравнений методом Гаусса
        r = matrix[11] / matrix[10];
        *y_centr = (matrix[7] - matrix[6] * r) / matrix[5];
        *x_centr = (matrix[3] - matrix[2] * r - matrix[1] * (*y_centr)) / matrix[0];

        // Вычисление результирующих значений
        r = sqrt((*x_centr) * (*x_centr) + (*y_centr) * (*y_centr) - 4.0 * r);
        *x_centr /= -2.0;
        *y_centr /= -2.0;
    }

    free(matrix);
    return r;
}

/***********************************************************
* Функция нахождения аппроксимирующей прямой в параметрическом виде, описываемой уравнением вида
* x = x_point + x_vector * t и y = y_point + y_vector * t, где t - параметр
* Входные данные:
* points - Указатель на массив точек с координатами
* хранящимися в двумерном массиве points с размерами points_num * 2 (строка {x,y}).
* points_num - Количество точек в массиве
* Выходные данные:
* x_point - Координата X начальной точки прямой (центр масс всех точек, при условии, что массы точек равны)
* y_point - Координата Y начальной точки прямой (центр масс всех точек, при условии, что массы точек равны)
* x_vector - Координата X направляющего вектора прямой
* y_vector - Координата Y направляющего вектора прямой
* r - Длина направляющего вектора
************************************************************/
float Line_approximation(float *points, unsigned short points_num, float tolerance, float *x_point, float *y_point, float *x_vector, float *y_vector)
{
    float r = 0.0;  // Длина направляющего вектора


    return r;
}
