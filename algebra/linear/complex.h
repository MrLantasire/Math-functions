#ifndef __COMPLEX_H
#define __COMPLEX_H

typedef struct 
{
    float real;
    float imag;
} complex32_t;

// Сложение комплексных чисел
extern complex32_t Complex_add(complex32_t a, complex32_t b);

// Вычитание комплексных чисел
extern complex32_t Complex_sub(complex32_t a, complex32_t b);

// Умножение комплексных чиел
extern complex32_t Complex_mul(complex32_t a, complex32_t b);

// Деление комплексных чиел
extern complex32_t Complex_div(complex32_t a, complex32_t b);

// Возведение комплексного числа в степень
extern complex32_t Complex_pow(complex32_t z, signed char degree);

// Извлечение корня из комплексного числа
extern complex32_t Complex_root(complex32_t z, signed char degree, unsigned char root_pos);

// Вычисление аргумента комплексного числа
extern float Complex_arg(complex32_t z);

// Вычисление модуля комплексного числа
extern float Complex_abs(complex32_t z);

// Вычисление комплексно-сопряжённого числа
extern complex32_t Complex_conjugate(complex32_t z);

#endif //__COMPLEX_H
