#include "complex.h"
#include "common.h"

/***********************************************************
* Функция сложения двух комплексных чисел.
* Входные данные:
* a - Первое комплексное слагаемое
* b - Второе комплексное слагаемое
* Выходные данные:
* Сумма комплексных чисел
************************************************************/
complex32_t Complex_add(complex32_t a, complex32_t b)
{
    return (complex32_t) {(a.real + b.real), (a.imag + b.imag)};
}

/***********************************************************
* Функция вычитания двух комплексных чисел.
* Входные данные:
* a - Комплексное уменьшаемое
* b - Комплексное вычитаемое
* Выходные данные:
* Разность комплексных чисел
************************************************************/
complex32_t Complex_sub(complex32_t a, complex32_t b)
{
    return (complex32_t) {(a.real - b.real), (a.imag - b.imag)};
}

/***********************************************************
* Функция умножения двух комплексных чисел.
* Входные данные:
* a - Первый комплексный множитель
* b - Всторой комплексный множитель
* Выходные данные:
* Произведение комплексных чисел
************************************************************/
complex32_t Complex_mul(complex32_t a, complex32_t b)
{
    return (complex32_t) {(a.real * b.real - a.imag * b.imag), (a.real * b.imag + a.imag * b.real)};
}

/***********************************************************
* Функция деления двух комплексных чисел.
* Входные данные:
* a - Комплексное делимое
* b - Комплексный делитель
* Выходные данные:
* Частное комплексных чисел
************************************************************/
complex32_t Complex_div(complex32_t a, complex32_t b)
{
    float denominator = (b.real * b.real + b.imag * b.imag);
    if (denominator > 0.0)
        return (complex32_t) {((a.real * b.real + a.imag * b.imag) / denominator), ((a.imag * b.real - a.real * b.imag) / denominator)};
    else
        return (complex32_t) {INFINITY, INFINITY};
}

/***********************************************************
* Функция возведения комплексного числа в степень.
* Входные данные:
* z - Комплексное число
* degree - Целочисленная степень
* Выходные данные:
* Комплексное число в степени
************************************************************/
complex32_t Complex_pow(complex32_t z, signed char degree)
{
    float abs = Complex_abs(z);
    float arg = Complex_arg(z);
    abs = powf(abs, degree);
    arg *= degree;
    return (complex32_t) {(abs * cos(arg)),(abs * sin(arg))};
}

/***********************************************************
* Функция извлечения корня степени degree из комплексного числа.
* Входные данные:
* z - Комплексное число
* degree - Целочисленная степень не равная 0
* root_pos - Номер корня комплексного числа, при равным 0 - главное значение корня
* Выходные данные:
* Корень из комплексного числа
************************************************************/
complex32_t Complex_root(complex32_t z, signed char degree, unsigned char root_pos)
{
    if (degree == 0)
        return (complex32_t ) {NAN, NAN};
        
    float abs = Complex_abs(z);
    float arg = Complex_arg(z);
    abs = powf(abs, (1.0 / degree));
    arg = (arg + 2.0 * M_PI * ((float) (root_pos % ABS(degree)))) / ((float) degree);
    return (complex32_t) {(abs * cos(arg)),(abs * sin(arg))};
}

/***********************************************************
* Функция вычисления аргумента (fi) комплексного числа.
* Входные данные:
* z - Комплексное число
* Выходные данные:
* Аргумент комплексного числа
************************************************************/
float Complex_arg(complex32_t z)
{
    return atan2(z.imag, z.real);
}

/***********************************************************
* Функция вычисления модуля (абсолютная величина) комплексного числа.
* Входные данные:
* z - Комплексное число
* Выходные данные:
* Модуль комплексного числа
************************************************************/
float Complex_abs(complex32_t z)
{
    return (sqrt((z.real * z.real + z.imag * z.imag)));
}

/***********************************************************
* Функция вычисления комплексно-сопряжённого числа.
* Входные данные:
* z - Комплексное число
* Выходные данные:
* Комплексно-сопряжённое число
************************************************************/
complex32_t Complex_conjugate(complex32_t z)
{
    return (complex32_t) {z.real, (-1. * z.imag)};
}

/***********************************************************
* Функция определения равенства комплексных чисел.
* Входные данные:
* a - Первое комплексное число
* b - Второе комплексное число
* tolerance - Допуск сравнения
* Выходные данные:
* True - если комплексные числа равны, False - комплексные числа не равны
************************************************************/
bool Is_complex_equal(complex32_t a, complex32_t b, float tolerance)
{
    return ( (ABS(a.real - b.real) <= tolerance) && (ABS(a.imag - b.imag) <= tolerance) );
}

/***********************************************************
* Функция определения комплексно ли число (имеет мнимую часть отличную от 0).
* Входные данные:
* a - Число в комплексном виде
* Выходные данные:
* True - число комплексное, False - число не комплексное
************************************************************/
extern bool Is_number_complex(complex32_t a)
{
    return (ABS(a.imag) > 0.0);
}
