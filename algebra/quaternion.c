#include "quaternion.h"
#include "common.h"

// Сложение кватернионов
quaternion32_t Quaternion_add(quaternion32_t a, quaternion32_t b)
{
    return (quaternion32_t) {(a.scalar + b.scalar), Vector_add(a.vector, b.vector)};
}

// Вычитание кватернионов
quaternion32_t Quaternion_sub(quaternion32_t a, quaternion32_t b)
{
    return (quaternion32_t) {(a.scalar - b.scalar), Vector_sub(a.vector, b.vector)};
}

// Умножение кватернионов
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

// Деление кватернионов
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

// Определение равенства кватернионов
bool Is_quaternion_equal(quaternion32_t a, quaternion32_t b, float tolerance)
{
    return ( (ABS(a.scalar - b.scalar) <= tolerance) && Is_vector_equal(a.vector, b.vector, tolerance) );
}

// Вычисление сопряженного кватерниона
quaternion32_t Quaternion_conjugate(quaternion32_t q)
{
    return (quaternion32_t) {q.scalar, Vector_mul(q.vector, -1.0)};
}

// Вычисление обратного кватерниона
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

// Вычисление нормализованного кватерниона
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

// Определение является ли кватернион чисто векторным
bool Is_quaternion_vector(quaternion32_t q)
{
    return (!(ABS(q.scalar) > 0.0) && (Vector_abs(q.vector) > 0.0));
}

// Определение является ли кватернион чисто скалярным
bool Is_quaternion_scalar(quaternion32_t q)
{
    return ((ABS(q.scalar) > 0.0) && !(Vector_abs(q.vector) > 0.0));
}

// Вычисление модуля (нормы) кватерниона
float Quaternion_abs(quaternion32_t q)
{
    return (sqrt(q.scalar * q.scalar + q.vector.i * q.vector.i + q.vector.j * q.vector.j + q.vector.k * q.vector.k));
}
