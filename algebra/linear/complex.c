#include "complex.h"
#include "common.h"

// Сложение комплексных чисел
complex32_t Complex_add(complex32_t a, complex32_t b)
{
    return (complex32_t) {(a.real + b.real), (a.imag + b.imag)};
}

// Вычитание комплексных чисел
complex32_t Complex_sub(complex32_t a, complex32_t b)
{
    return (complex32_t) {(a.real - b.real), (a.imag - b.imag)};
}

// Умножение комплексных чиел
complex32_t Complex_mul(complex32_t a, complex32_t b)
{
    return (complex32_t) {(a.real * b.real - a.imag * b.imag), (a.real * b.imag + a.imag * b.real)};
}

// Деление комплексных чиел
complex32_t Complex_div(complex32_t a, complex32_t b)
{
    float denominator = (b.real * b.real + b.imag * b.imag);
    return (complex32_t) {((a.real * b.real + a.imag * b.imag) / denominator), ((a.imag * b.real - a.real * b.imag) / denominator)};
}

// Возведение комплексного числа в степень
complex32_t Complex_pow(complex32_t z, signed char degree)
{
    float abs = Complex_abs(z);
    float arg = Complex_arg(z);
    abs = powf(abs, degree);
    arg *= degree;
    return (complex32_t) {(abs * cos(arg)),(abs * sin(arg))};
}

// Извлечение корня из комплексного числа
complex32_t Complex_root(complex32_t z, signed char degree, unsigned char root_pos)
{
    float abs = Complex_abs(z);
    float arg = Complex_arg(z);
    abs = powf(abs, (1.0 / degree));
    arg = (arg + 2.0 * M_PI * ((float) (root_pos % ABS(degree)))) / ((float) degree);
    return (complex32_t) {(abs * cos(arg)),(abs * sin(arg))};
}

// Вычисление аргумента комплексного числа
float Complex_arg(complex32_t z)
{
    return atan2(z.imag, z.real);
}

// Вычисление модуля комплексного числа
float Complex_abs(complex32_t z)
{
    return (sqrt((z.real * z.real + z.imag * z.imag)));
}

// Вычисление комплексно-сопряжённого числа
complex32_t Complex_conjugate(complex32_t z)
{
    return (complex32_t) {z.real, (-1. * z.imag)};
}
