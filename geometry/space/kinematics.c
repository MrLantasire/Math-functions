#include <malloc.h>
#include "kinematics.h"
#include "common.h"
#include "matrix.h"
#include "equation.h"

// Определение пространственных углов между векторами вокруг заданных осей
unsigned char Define_spatial_angles_of_vectors( vector32_3D_t initial_vec, vector32_3D_t final_vec, vector32_3D_t first_axis, 
                                                vector32_3D_t second_axis, float out_first_angle[2], float out_second_angle[2])
{
    unsigned char out = 0;                      // 1 - есть решение, 0 - нет решения

    vector32_3D_t cross_vec = {0};              // Вектор-перпендикуляр осям вращения
    vector32_3D_t inter_vec = {0};              // Промежуточный вектор
    vector32_3D_t projected_vec = {0};          // Спроецированный вектор
    vector32_3D_t cross_point = {0};            // Координаты точки пересечения плоскостей
    float *matrix = calloc(12, sizeof(float));  // Матрица системы уравнений

    // Обнуление выходных значений
    out_first_angle[0] = 0.0;
    out_first_angle[1] = 0.0;
    out_second_angle[0] = 0.0;
    out_second_angle[1] = 0.0;

    // Нормализация входных векторов
    initial_vec = Vector_normalize(initial_vec);
    final_vec = Vector_normalize(final_vec);
    first_axis = Vector_normalize(first_axis);
    second_axis = Vector_normalize(second_axis);

    // Вычисление вектора-перпендикуляра осям вращения
    cross_vec = Vector_cross_mul(first_axis, second_axis);
    cross_vec = Vector_normalize(cross_vec);

    // Заполнение матрицы
    // Первая плоскость
    matrix[0] = first_axis.i;
    matrix[1] = first_axis.j;
    matrix[2] = first_axis.k;
    matrix[3] = first_axis.i * initial_vec.i + first_axis.j * initial_vec.j + first_axis.k * initial_vec.k;
    // Вторая плоскость
    matrix[4] = second_axis.i;
    matrix[5] = second_axis.j;
    matrix[6] = second_axis.k;
    matrix[7] = second_axis.i * final_vec.i + second_axis.j * final_vec.j + second_axis.k * final_vec.k;
    // Третья плоскость
    matrix[8] = cross_vec.i;
    matrix[9] = cross_vec.j;
    matrix[10] = cross_vec.k;
    matrix[11] = 0.0;

    // Если ранг матрицы меньше 3-х, то нет точки пересечения плоскостей (плоскости параллельны, или совпадают, или вырождены)
    // Возможны частные решения для заданного набора входных данных (количество решений либо 0, либо бесконечность)
    if (Gaussian_elimination(matrix, 3U, 4U,  1.E-37, matrix) == 3U)
    {
        // Вычисление координат точки пересечения плоскостей
        cross_point.k = matrix[11] / matrix[10];
        cross_point.j = (matrix[7] - matrix[6] * cross_point.k) / matrix[5];
        cross_point.i = (matrix[3] - matrix[2] * cross_point.k - matrix[1] * cross_point.j) / matrix[0]; 

        // Заполнение коэффициентов квадратного уравнения
        matrix[0] = cross_vec.i * cross_vec.i + cross_vec.j * cross_vec.j + cross_vec.k * cross_vec.k;
        matrix[1] = 2.0 * (cross_vec.i * cross_point.i + cross_vec.j * cross_point.j + cross_vec.k * cross_point.k);
        matrix[2] = -1.0 + cross_point.i * cross_point.i + cross_point.j * cross_point.j + cross_point.k * cross_point.k;

        complex32_t *solutions = calloc(2, sizeof(complex32_t));    // Матрица для решений квадратного уравнения

        // Если только комплексные корни, значит нет решений
        if (!(Solve_quadratic_equation(matrix[0], matrix[1], matrix[2], solutions) < 0.0))
        {
            // Дублирование значений начального и конечного векторов
            matrix[0] = initial_vec.i;
            matrix[1] = initial_vec.j;
            matrix[2] = initial_vec.k;
            matrix[3] = final_vec.i;
            matrix[4] = final_vec.j;
            matrix[5] = final_vec.k;

            for (unsigned char i = 0; i < 2; i++)
            {
                // Координаты вектора между первым и вторым поворотами
                inter_vec = Vector_mul(cross_vec, solutions[i].real);
                inter_vec = Vector_add(cross_point, inter_vec);

                // Проецирование вектора
                projected_vec = Project_vector_onto_normal(inter_vec, first_axis);
                initial_vec = Project_vector_onto_normal(initial_vec, first_axis);
                out_first_angle[i] = Define_angle_of_vectors(initial_vec, projected_vec);

                if (Vector_mixed_mul(initial_vec, projected_vec, first_axis) < 0.0)
                    out_first_angle[i] *= -1.0;

                projected_vec = Project_vector_onto_normal(inter_vec, second_axis);
                final_vec = Project_vector_onto_normal(final_vec, second_axis);
                out_second_angle[i] = Define_angle_of_vectors(final_vec, projected_vec);

                if (Vector_mixed_mul(projected_vec, final_vec, second_axis) < 0.0)
                    out_second_angle[i] *= -1.0;
                
                // Возвращение значений в векторы
                initial_vec.i = matrix[0];
                initial_vec.j = matrix[1];
                initial_vec.k = matrix[2];
                final_vec.i = matrix[3];
                final_vec.j = matrix[4];
                final_vec.k = matrix[5];
            }

            // Есть решение
            out = 1U;
        }

        free(solutions);
    }

    free(matrix);
    return out;
}
