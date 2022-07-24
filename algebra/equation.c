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
    complex32_t complex_Q = {0.0,0.0};
    complex32_t alpha = {0.0,0.0};
    complex32_t betta = {0.0,0.0};
    float p = (3.0 * a * c - b * b) / (3.0 * a * a);            // Коэффициент p = (3ac - b^2)/(3a^2)
                                                                // Коэффициент q = (2b^3 - 9abc + 27a^2d)/(27a^3)
    float q = (2.0 * b * b * b - 9.0 * a * b * c +
                27.0 * a * a * d) / (27.0 * a * a * a);
    complex32_t complex_var = {0.0,0.0};                        // Комплексная переменная для промежуточных вычислений

    p /= 3.0;
    q /= 2.0;
    complex_Q.real = (p * p * p + q * q);
    discriminant = -108.0 * complex_Q.real;
    complex_Q = Complex_root(complex_Q, 2, 0);
    complex_var.real = -q;

    // Вычисление корней
    alpha = Complex_root(Complex_add(complex_var, complex_Q), 3, 0);

    float min_dif = 0.0;
    for (unsigned char i = 0; i < 3; i++)
    {
        complex32_t var_betta = Complex_root(Complex_sub(complex_var, complex_Q), 3, i);
        float dif = ABS(Complex_mul(alpha, var_betta).real + p);
    
        if ( (dif < min_dif) || (i == 0) )
        {
            betta = var_betta;
            min_dif = dif;
        }
    }

    out[0] = Complex_add(alpha, betta);

    out[1] = Complex_div(out[0], (complex32_t) {-2.0,0.0});
    out[2] = out[1];

    complex_var.real = 0.0;
    complex_var.imag = sqrt(3.0) / 2.0;
    complex_var = Complex_mul(Complex_sub(alpha, betta), complex_var);

    out[1] = Complex_add(out[1], complex_var);
    out[2] = Complex_sub(out[2], complex_var);

    complex_var.real = b / (3.0 * a);
    complex_var.imag = 0.0;
    for (unsigned char i = 0; i < 3; i++)
    {
        out[i] = Complex_sub(out[i], complex_var);
    }

    return discriminant;
}
