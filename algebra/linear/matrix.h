#ifndef __MATRIX_H
#define __MATRIX_H

// Функция преобразования матрицы к треугольному виду
extern unsigned char Gaussian_elimination(float *input_matrix, unsigned char num_rows, unsigned char num_colums, float tolerance, float *output_matrix);

// Функция вычисления определителя матрицы n*n
extern float Define_determinant(float *matrix, unsigned char size, float tolerance);

// Функция умножения матриц
extern void Matrices_mul(float *first_matrix, float *second_matrix, unsigned char num_out_rows, unsigned char num_common, unsigned char num_out_colums, float *out_matrix);

// Функция умножения матрицы на число
extern void Matrix_mul_by_num(float *matrix, unsigned char num_rows, unsigned char num_colums, float factor);

#endif //__MATRIX_H
