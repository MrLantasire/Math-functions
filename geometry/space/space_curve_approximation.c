#include <malloc.h>
#include "space_curve_approximation.h"
#include "common.h"
#include "algebra\linear\matrix.h"
#include "algebra\equation.h"

/***********************************************************
* Функция нахождения аппроксимирующей прямой в параметрическом виде, описываемой уравнением вида
* x = x_point + x_vector * t и y = y_point + y_vector * t и z = z_point + z_vector * t,  где t - параметр
* Входные данные:
* points - Указатель на массив точек с координатами
* хранящимися в двумерном массиве points с размерами points_num * 3 (строка {x,y,z}).
* points_num - Количество точек в массиве
* Выходные данные:
* lines_param - Указатель на массив параметров аппроксимирующей прямой размером 2 * 3, где
* первая строка - координаты начальной точки {x,y,z}; вторая строка - координаты направляющего вектора {x,y,z}
* return - Длина направляющего вектора
************************************************************/
float Line_approximation_3D(float *points, unsigned short points_num, float *lines_param)
{
    float *matrix = calloc(12, sizeof(float));   // Матрица системы уравнений
    complex32_t *solutions = calloc(3, sizeof(complex32_t));

    // Обнуление матрицы и выходных параметров
    for(unsigned short i = 0; i < 6; i++)
    {
        matrix[2 * i] = 0.0;
        matrix[2 * i + 1] = 0.0;
        lines_param[i] = 0.0;
    }

    // Определение координат начальной точки прямой как среднее арифметическое всех координат
    for (unsigned short i = 0; i < points_num; i++)
    {
        lines_param[0] += points[i * 3];
        lines_param[1] += points[i * 3 + 1];
        lines_param[2] += points[i * 3 + 2];
    }
    lines_param[0] /= points_num;
    lines_param[1] /= points_num;
    lines_param[2] /= points_num;

    // Вычисление тензора моментов инерции
    for (unsigned short i = 0; i < points_num; i++)
    {
        // Вспомогательные вычисления
        matrix[9] = points[i * 3] - lines_param[0];
        matrix[10] = points[i * 3 + 1] - lines_param[1];
        matrix[11] = points[i * 3 + 2] - lines_param[2];
        matrix[2] = matrix[9] * matrix[9] + matrix[10] * matrix[10] + matrix[11] * matrix[11];
        for (unsigned short j = 0; j < 3; j++)
        {
            matrix[j * 3 + j] += matrix[2] - matrix[9 + j] * matrix[9 + j];
            matrix[j * 3 + (j + 1) % 3] -= matrix[9 + j] * matrix[9 + (j + 1) % 3];
        }
    }
    matrix[2] = matrix[6];
    matrix[3] = matrix[1];
    matrix[7] = matrix[5];

    // Вычисление инвариантов
    matrix[9] = 0.0;
    matrix[10] = 0.0;
    matrix[11] = 0.0;
    for (unsigned short i = 0; i < 3; i++)
    {
        // Первый инвариант
        matrix[9] += matrix[i * 3 + i];
        // Второй инвариант
        matrix[10] += matrix[i * 3 + i] * matrix[((i + 1) % 3) * 4];
        matrix[10] -= matrix[i * 3 + (i + 1) % 3] * matrix[i * 3 + (i + 1) % 3];
        // Третий инвариант
        matrix[11] += matrix[i * 3] * matrix[((i + 1) % 3) * 3 + 1] * matrix[((i + 2) % 3) * 3 + 2];
        matrix[11] -= matrix[i * 3] * matrix[((i + 1) % 3) * 3 + 2] * matrix[((i + 2) % 3) * 3 + 1];
    }

    if(!(Solve_cubic_equation(-1.0, matrix[9], (-1.0 * matrix[10]), matrix[11], solutions) < 0.0))
    {
        // Решениями системы будут 3 момента инерции (выбирается минимальный)
        matrix[9] = ABS(solutions[0].real) > ABS(solutions[1].real) ? solutions[1].real : solutions[0].real;
        if(ABS(matrix[9]) > ABS(solutions[2].real)) 
            matrix[9] = solutions[2].real;

        // Изменение исходной матрицы (преобразование диагональных элементов)
        matrix[0] -= matrix[9]; 
        matrix[4] -= matrix[9]; 
        matrix[8] -= matrix[9]; 

        // Вычисление направляющего вектора прямой
        if (Gaussian_elimination(matrix, 3U, 3U, 1.E-37,matrix))
        {
            if(!(ABS(matrix[0]) > 0.0))
            {
                lines_param[3] = 1.0;
                lines_param[4] = 0.0;
                lines_param[5] = 0.0;
            }
            else
            {
                if (!(ABS(matrix[4]) > 0.0))
                {
                    lines_param[4] = 1.0;
                    lines_param[5] = 0.0;
                }
                else
                {
                    lines_param[5] = 1.0;
                    lines_param[4] = -1.0 * lines_param[5] * matrix[5] / matrix[4];
                }

                lines_param[3] = -1.0 * (lines_param[4] * matrix[1] + lines_param[5] * matrix[2]) / matrix[0];
            }
        }
        // Если ранг матрицы равен нулю, то бесконечное множество решений
        else
        {
            // Как вариант одного из бесконечного множества решений
            lines_param[3] = 1.0;
            lines_param[4] = 0.0;
            lines_param[5] = 0.0;
        }
    }

    free(matrix);
    free(solutions);
    return sqrt(lines_param[3] * lines_param[3] + lines_param[4] * lines_param[4] + lines_param[5] * lines_param[5]);
}
