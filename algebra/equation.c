#include "equation.h"
#include "common.h"

/***********************************************************
* Функция решения крадратного уравнения вида ax^2 + bx + c = 0.
* Входные данные:
* a - Первый множитель
* b - Второй множитель
* c - Свободный член
* Выходные данные:
* out[] - Массив корней уравнения в комплексном виде
* discriminant - Дискриминант уравнения
************************************************************/
float Solve_quadratic_equation(float a, float b, float c, complex32_t out[2])
{
    float discriminant = b * b - 4.0 * a * c;                                           // Дискриминант уравнения
    complex32_t sqrt_of_discr = Complex_root((complex32_t) {discriminant, 0.0}, 2, 0);  // Корень из дискриминанта в комплексном виде
    complex32_t complex_factors_division = {(b / (-2.0 * a)), 0.0};                     // -b/2a в комплексном виде

    // Деление корня из дискриминанта на 2a
    sqrt_of_discr = Complex_div(sqrt_of_discr, (complex32_t) {(2.0 * a), 0.0});

    // Вычисление корней уравнения в комплексном виде
    out[0] = Complex_add(complex_factors_division, sqrt_of_discr);
    out[1] = Complex_sub(complex_factors_division, sqrt_of_discr);

    return discriminant;
}

/***********************************************************
* Функция решения кубического уравнения вида ax^3 + bx^2 + cx + d =0
* с помощью формулы Кардано.
* Входные данные:
* a - Первый множитель
* b - Второй множитель
* c - Третий множитель
* d - Свободный член
* Выходные данные:
* out[] - Массив корней уравнения в комплексном виде
* discriminant - Дискриминант уравнения
************************************************************/
float Solve_cubic_equation(float a, float b, float c, float d, complex32_t out[3])
{
    float discriminant = 0.0;

    return discriminant;
}
