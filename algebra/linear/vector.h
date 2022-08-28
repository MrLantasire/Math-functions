#ifndef __VECTOR_H
#define __VECTOR_H

#include <stdbool.h>

// Структура вектора в 3-х мерном пространстве
typedef struct
{
    float i;
    float j;
    float k;
} vector32_3D_t;

// Сложение векторов
extern vector32_3D_t Vector_add(vector32_3D_t a, vector32_3D_t b);

// Вычитание векторов
extern vector32_3D_t Vector_sub(vector32_3D_t a, vector32_3D_t b);

// Умножение вектора на число
extern vector32_3D_t Vector_mul(vector32_3D_t a, float b);

// Векторное произведение векторов
extern vector32_3D_t Vector_cross_mul(vector32_3D_t a, vector32_3D_t b);

// Скалярное произведение векторов
extern float Vector_scalar_mul(vector32_3D_t a, vector32_3D_t b);

// Смешанное произведение векторов
extern float Vector_mixed_mul(vector32_3D_t a, vector32_3D_t b, vector32_3D_t c);

// Вычисление длины вектора
extern float Vector_abs(vector32_3D_t a);

// Вычисление нормализованного вектора
extern vector32_3D_t Vector_normalize(vector32_3D_t a);

// Вычисление угла между векторами
extern float Define_angle_of_vectors(vector32_3D_t a, vector32_3D_t b);

// Вычисление матрицы поворота вокруг заданного вектора
extern void Get_rotation_matrix(vector32_3D_t n, float angle, float *out_matrix);

// Определение равенства векторов
extern bool Is_vector_equal(vector32_3D_t a, vector32_3D_t b, float tolerance);

#endif //__VECTOR_H
