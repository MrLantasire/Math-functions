#include <malloc.h>
#include "plane_curve_approximation.h"
#include "common.h"
#include "algebra\linear\matrix.h"
#include "algebra\equation.h"

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
float Circle_approximation(float *points, unsigned short points_num, float *x_centr, float *y_centr)
{
    float r = 0.0;                              // Радиус окружности
    float *matrix = calloc(12, sizeof(float));  // Матрица системы уравнений

    // Обнуление выходных параметров 
    *x_centr = 0.0;
    *y_centr = 0.0;

    // Обнуление множителей
    for (unsigned short i = 0; i < 12U; i++)
    {
        matrix[i] = 0.0;
    }

    // Вычисление множителей для системы уравнений
    for (unsigned short i = 0; i < points_num; i++)
    {
        matrix[4] = points[i * 2] * points[i * 2];                  // Использование для промежуточных вычислений (xi^2)
        matrix[8] = points[i * 2 + 1] * points[i * 2 + 1];          // Использование для промежуточных вычислений (yi^2)

        matrix[0] += matrix[4];                                     // x^2
        matrix[1] += points[i * 2] * points[i * 2 + 1];             // xy
        matrix[2] += points[i * 2];                                 // x
        matrix[3] -= points[i * 2] * (matrix[4] + matrix[8]);       // -x(x^2 + y^2)
        matrix[5] += matrix[8];                                     // y^2
        matrix[6] += points[i * 2 + 1];                             // y
        matrix[7] -= points[i * 2 + 1] * (matrix[4] + matrix[8]);   // -y(x^2 + y^2)
        matrix[11] -= matrix[4] + matrix[8];                        // -(x^2 + y^2)
    }
    matrix[4] = matrix[1];                                          // xy
    matrix[8] = matrix[2];                                          // x
    matrix[9] = matrix[6];                                          // y
    matrix[10] = (float) points_num;                                // n

    if (Gaussian_elimination(matrix, 3, 4, 1.E-37, matrix) == 3U)
    {
        // Решение системы уравнений методом Гаусса
        r = matrix[11] / matrix[10];
        *y_centr = (matrix[7] - matrix[6] * r) / matrix[5];
        *x_centr = (matrix[3] - matrix[2] * r - matrix[1] * (*y_centr)) / matrix[0];

        // Вычисление результирующих значений
        r = sqrt((*x_centr) * (*x_centr) + (*y_centr) * (*y_centr) - 4.0 * r) / 2.0;
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
* return - Длина направляющего вектора
************************************************************/
float Line_approximation(float *points, unsigned short points_num, float *x_point, float *y_point, float *x_vector, float *y_vector)
{
    float *matrix = calloc(6, sizeof(float));   // Матрица системы уравнений
    complex32_t *solutions = calloc(2, sizeof(complex32_t));

    // Обнуление выходных параметров 
    *x_point = 0.0;
    *y_point = 0.0;
    *x_vector = 0.0;
    *y_vector = 0.0;

    // Обнуление матрицы
    for(unsigned short i = 0; i < 6; i++)
    {
        matrix[i] = 0.0;
    }

    // Определение координат начальной точки прямой как среднее арифметическое всех координат
    for (unsigned short i = 0; i < points_num; i++)
    {
        *x_point += points[i * 2];
        *y_point += points[i * 2 + 1];
    }
    *x_point /= points_num;
    *y_point /= points_num;

    // Вычисление тензора моментов инерции
    for (unsigned short i = 0; i < points_num; i++)
    {
        // Вспомогательные вычисления
        matrix[4] = points[i * 2] - (*x_point);
        matrix[5] = points[i * 2 + 1] - (*y_point);
        for (unsigned short j = 0; j < 2; j++)
        {
            matrix[j * 2 + j] += matrix[4] * matrix[4] + matrix[5] * matrix[5] - matrix[4 + j] * matrix[4 + j];
            matrix[1 + j] -= matrix[4] * matrix[5];
        }
    }

    if(!(Solve_quadratic_equation(-1.0, (matrix[0] + matrix[3]), (matrix[1] * matrix[2] - matrix[0] * matrix[3]), solutions ) < 0.0))
    {
        // Решениями системы будут максимальный и минимальный моменты инерции (выбирается минимальный)
        if (ABS(solutions[0].real) > ABS(solutions[1].real))
        {
            matrix[0] -= solutions[1].real;
            matrix[3] -= solutions[1].real;
        }
        else
        {
            matrix[0] -= solutions[0].real;
            matrix[3] -= solutions[0].real;
        }

        // Вычисление направляющего вектора прямой
        if (Gaussian_elimination(matrix, 2U, 2U, 1.E-37,matrix))
        {
            if(!(ABS(matrix[1]) > 0.0))
            {
                *x_vector = 0.0;
                *y_vector = 1.0;
            }
            else
            {
                if(!(ABS(matrix[0]) > 0.0))
                {
                    *x_vector = 1.0;
                    *y_vector = 0.0;
                }
                else
                {
                    *y_vector = 1.0;
                    *x_vector = -1.0 * (*y_vector) * matrix[1] / matrix[0];
                }
            }
        }
        // Если ранг матрицы равен нулю, то бесконечное множество решений
        else
        {
            // Как вариант одного из бесконечного множества решений
            *x_vector = 1.0;
            *y_vector = 0.0;
        }
    }

    free(matrix);
    free(solutions);
    return sqrt((*x_vector) * (*x_vector) + (*y_vector) * (*y_vector));
}
