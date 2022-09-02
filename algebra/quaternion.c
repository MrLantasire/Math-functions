#include "quaternion.h"
#include "common.h"

/***********************************************************
* Функция сложения двух кватернионов.
* Входные данные:
* a - Первый кватернион
* b - Второй кватернион
* Выходные данные:
* Сумма кватернионов
************************************************************/
quaternion32_t Quaternion_add(quaternion32_t a, quaternion32_t b)
{
    return (quaternion32_t) {(a.scalar + b.scalar), Vector_add(a.vector, b.vector)};
}

/***********************************************************
* Функция вычитания кватернионов.
* Входные данные:
* a - Уменьшаемое
* b - Вычитаемое
* Выходные данные:
* Разность кватернионов
************************************************************/
quaternion32_t Quaternion_sub(quaternion32_t a, quaternion32_t b)
{
    return (quaternion32_t) {(a.scalar - b.scalar), Vector_sub(a.vector, b.vector)};
}

/***********************************************************
* Функция умножения кватернионов.
* Входные данные:
* a - Первый кватернион
* b - Второй кватернион
* Выходные данные:
* Произведение кватернионов
************************************************************/
quaternion32_t Quaternion_mul(quaternion32_t a, quaternion32_t b)
{
    quaternion32_t c = {0};
    c.scalar = a.scalar * b.scalar - Vector_scalar_mul(a.vector, b.vector);
    c.vector = Vector_cross_mul(a.vector, b.vector);
    a.vector = Vector_mul(a.vector, b.scalar);
    b.vector = Vector_mul(b.vector, a.scalar);
    c.vector = Vector_add(c.vector, a.vector);
    c.vector = Vector_add(c.vector, b.vector);
    return c;
}

/***********************************************************
* Функция деления кватернионов.
* Входные данные:
* a - Делимое
* b - Делитель
* Выходные данные:
* Частное кватернионов
************************************************************/
quaternion32_t Quaternion_div(quaternion32_t a, quaternion32_t b)
{
    float denominator = Quaternion_abs(b);
    denominator *= denominator;

    if (denominator > 0.0)
    {
        b = Quaternion_conjugate(b);
        a = Quaternion_mul(a, b);
        a.scalar /= denominator;
        a.vector.i /= denominator;
        a.vector.j /= denominator;
        a.vector.k /= denominator;
        return a;
    }
    else
    {
        return (quaternion32_t) {INFINITY, (vector32_3D_t) {INFINITY, INFINITY, INFINITY}};
    }
}

/***********************************************************
* Функция сравнения кватернионов.
* Входные данные:
* a - Первый кватернион
* b - Второй кватернион
* tolerance - Допуск сравнения
* Выходные данные:
* True - если кватернионы равны, False - кватернионы не равны
************************************************************/
bool Is_quaternion_equal(quaternion32_t a, quaternion32_t b, float tolerance)
{
    return ( (ABS(a.scalar - b.scalar) <= tolerance) && Is_vector_equal(a.vector, b.vector, tolerance) );
}

/***********************************************************
* Функция определения сопряженного кватерниона.
* Входные данные:
* q - Кватернион
* Выходные данные:
* Сопряженный кватернион
************************************************************/
quaternion32_t Quaternion_conjugate(quaternion32_t q)
{
    return (quaternion32_t) {q.scalar, Vector_mul(q.vector, -1.0)};
}

/***********************************************************
* Функция определения обратного кватерниона.
* Входные данные:
* q - Кватернион
* Выходные данные:
* Обратный кватернион (1 / q)
************************************************************/
quaternion32_t Quaternion_reciprocal(quaternion32_t q)
{
    float denominator = q.scalar * q.scalar + q.vector.i * q.vector.i + q.vector.j * q.vector.j + q.vector.k * q.vector.k;
    
    if (denominator > 0.0)
    {
        q = Quaternion_conjugate(q);
        q.scalar /= denominator;
        q.vector.i /= denominator;
        q.vector.j /= denominator;
        q.vector.k /= denominator;
    }
    
    return q;
}

/***********************************************************
* Функция вычисления нормализованного кватерниона.
* Входные данные:
* q - Кватернион
* Выходные данные:
* Нормализованный кватернион |q| = 1
************************************************************/
extern quaternion32_t Quaternion_normalize(quaternion32_t q)
{
    float denominator = Quaternion_abs(q);
    
    if (denominator > 0.0)
    {
        q.scalar /= denominator;
        q.vector.i /= denominator;
        q.vector.j /= denominator;
        q.vector.k /= denominator;
    }
    
    return q;
}

/***********************************************************
* Функция определения является ли кватернион чисто векторным.
* Входные данные:
* q - Кватернион
* Выходные данные:
* True - если скалярная часть = 0, а векторная != 0
************************************************************/
bool Is_quaternion_vector(quaternion32_t q)
{
    return (!(ABS(q.scalar) > 0.0) && (Vector_abs(q.vector) > 0.0));
}

/***********************************************************
* Функция определения является ли кватернион чисто скалярным.
* Входные данные:
* q - Кватернион
* Выходные данные:
* True - если скалярная часть != 0, а векторная = 0
************************************************************/
bool Is_quaternion_scalar(quaternion32_t q)
{
    return ((ABS(q.scalar) > 0.0) && !(Vector_abs(q.vector) > 0.0));
}

/***********************************************************
* Функция вычисления нормы кватерниона.
* Входные данные:
* q - Кватернион
* Выходные данные:
* Норма кватерниона |q|
************************************************************/
float Quaternion_abs(quaternion32_t q)
{
    return (sqrt(q.scalar * q.scalar + q.vector.i * q.vector.i + q.vector.j * q.vector.j + q.vector.k * q.vector.k));
}

/***********************************************************
* Функция кватерниона поворота ({cos(a), v*sin(a)}).
* Входные данные:
* axis - Вектор оси поворота
* angle - Угол поворота [rad]
* Выходные данные:
* Кватернион поворота
************************************************************/
quaternion32_t Quaternion_rotor(vector32_3D_t axis, float angle)
{
    axis = Vector_normalize(axis);
    axis = Vector_mul(axis, sin(angle));
    return (quaternion32_t) {cos(angle), axis};
}

/***********************************************************
* Функция поворота кватерниона на угол.
* Входные данные:
* q - Исходный кватернион
* axis - Вектор оси поворота
* angle - Угол поворота [rad]
* Выходные данные:
* Повернутый кватернион
************************************************************/
quaternion32_t Rotate_quaternion(quaternion32_t q, vector32_3D_t axis, float angle)
{
    quaternion32_t r = Quaternion_rotor(axis, (angle / 2.0));
    q = Quaternion_mul(r,q);
    r = Quaternion_reciprocal(r);
    q = Quaternion_mul(q,r);
    return q;
}
