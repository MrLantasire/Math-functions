#ifndef __MATRIX_H
#define __MATRIX_H

// Функция преобразования матрицы к треугольному виду
extern unsigned char Gaussian_elimination(float *input_matrix, unsigned char num_rows, unsigned char num_colums, float tolerance, float *output_matrix);

// Функция вычисления определителя матрицы n*n
extern float Define_determinant(float *matrix, unsigned char size, float tolerance);

#endif //__MATRIX_H
