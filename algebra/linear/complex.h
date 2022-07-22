#ifndef __COMPLEX_H
#define __COMPLEX_H

typedef struct 
{
    float real;
    float imag;
} complex32_t;

// Сложение комплексных чисел
extern complex32_t Complex_sum(complex32_t a, complex32_t b);

// Умножение комплексных чиел
extern complex32_t Complex_mul(complex32_t a, complex32_t b);

#endif //__COMPLEX_H
