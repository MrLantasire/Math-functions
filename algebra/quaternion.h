#ifndef __QUATERNION_H
#define __QUATERNION_H

#include <stdbool.h>
#include "vector.h"

// Структура кватерниона
typedef struct 
{
    float scalar;
    vector32_3D_t vector;
} quaternion32_t;

// Сложение кватернионов
extern quaternion32_t Quaternion_add(quaternion32_t a, quaternion32_t b);

// Вычитание кватернионов
extern quaternion32_t Quaternion_sub(quaternion32_t a, quaternion32_t b);

// Умножение кватернионов
extern quaternion32_t Quaternion_mul(quaternion32_t a, quaternion32_t b);

// Деление кватернионов
extern quaternion32_t Quaternion_div(quaternion32_t a, quaternion32_t b);

// Определение равенства кватернионов
extern bool Is_quaternion_equal(quaternion32_t a, quaternion32_t b, float tolerance);

// Вычисление сопряженного кватерниона
extern quaternion32_t Quaternion_conjugate(quaternion32_t q);

// Вычисление обратного кватерниона
extern quaternion32_t Quaternion_reciprocal(quaternion32_t q);

// Вычисление нормализованного кватерниона
extern quaternion32_t Quaternion_normalize(quaternion32_t q);

// Определение является ли кватернион чисто векторным
extern bool Is_quaternion_vector(quaternion32_t q);

// Определение является ли кватернион чисто скалярным
extern bool Is_quaternion_scalar(quaternion32_t q);

// Вычисление модуля (нормы) кватерниона
extern float Quaternion_abs(quaternion32_t q);

#endif //__QUATERNION_H
