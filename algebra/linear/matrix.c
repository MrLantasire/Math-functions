#include <malloc.h>
#include "matrix.h"
#include "common.h"

static void Copy_row(void *input, void *output, size_t size)
{
    for (size_t i = 0; i < size; i++)
    {
        *((char*) output + i) = *((char*) input + i);
    }
}

static void Swap_row(void *input, void *output, size_t size)
{
    char inter = 0;
    for (size_t i = 0; i < size; i++)
    {
       inter = *((char*) output + i);
       *((char*) output + i) = *((char*) input + i);
       *((char*) input + i) = inter;
    }
}

/***********************************************************
* Функция преобразования матрицы к треугольному виду.
* Входные данные:
* input_matrix - Указатель на входящую матрицу
* num_rows - Количество строк матрицы
* num_colums - Количество столбцов матрицы
* tolerance - Допуск
* Выходные данные:
* output_matrix - Указатель на выходящую матрицу
* rank - Ранг матрицы
************************************************************/
unsigned char Gaussian_elimination(float *input_matrix, unsigned char num_rows, unsigned char num_colums, float tolerance, float *output_matrix)
{   
    unsigned char rank = 0U;    // Переменная ранга матрицы
    unsigned char offset = 0U;  // Переменная смещения по столбцу матрицы (смещение)
    unsigned char shift = 0U;   // Переменная смещения по строке матрицы (сдвиг)
    float denominator = 0.0;    // Знаменатель для вычитания строк
    float multiplier = 0.0;     // Множитель для вычитания строк
    char exchange_flag = 1;     // Флаг перестановки строк

    // Копирование матрицы
    if (input_matrix != output_matrix)
    {
        for(unsigned char i = 0U; i < num_rows; i++)
        {
            Copy_row((input_matrix + i * num_colums), (output_matrix + i * num_colums), (sizeof(float) * num_colums));
        }
    }

    for(unsigned char i = 0U; i < num_rows; i++)
    {
        denominator = 0.0;
        offset = 0U;
        exchange_flag = 1;

        // Поиск строки с первым значением не равным 0
        while ((ABS(denominator) < tolerance) && ((shift + i) < num_colums))
        {
            denominator = output_matrix[i + shift + (offset + i) * num_colums];

            if((ABS(denominator) < tolerance))
                offset++;
            else
                rank++;

            if (offset + i >= num_rows)
            {
                shift++;
                offset = 0U;
            }
        }

        // Обмен строк местами
        if((offset != 0U) && ((offset + i) < num_rows))
        {
            Swap_row((output_matrix + i * num_colums), (output_matrix + (offset + i) * num_colums), sizeof(float) * num_colums);
            exchange_flag *= -1;
        }

        // Вычитание текущей строки из последующих строк
        for(unsigned char j = i; j < num_rows; j++)
        {
            if (j == i)
            {
                // Обнуление всех значений меньших tolerance в тукущей строке
                for (unsigned char k = 0U; k < num_colums; k++)
                {
                    // Обнуление значений меньших tolerance
                    if (ABS(output_matrix[k + j * num_colums]) < tolerance)
                        output_matrix[k + j * num_colums] = 0.0;
                    // Домножение на exchange_flag строки
                    else
                        output_matrix[k + j * num_colums] *= (float) exchange_flag;
                }
            }
            else
            {
                multiplier = (float) exchange_flag * output_matrix[i + shift + j * num_colums];

                if(!(ABS(multiplier) < tolerance))
                {
                    for (unsigned char k = 0U; k < num_colums; k++)
                    {
                        // Вычисление следующей строчки матрицы
                        output_matrix[k + j * num_colums] -= output_matrix[k + i * num_colums] * multiplier / denominator;
                        // Обнуление значений меньших tolerance
                        if (ABS(output_matrix[k + j * num_colums]) < tolerance)
                            output_matrix[k + j * num_colums] = 0.0;
                    }
                }
            }
        }
    }

    return rank;
}

/***********************************************************
* Функция вычисления определителя матрицы n*n.
* Входные данные:
* matrix - Указатель на матрицу
* size - Размер матрицы
* Выходные данные:
* out - Определитель матрицы
************************************************************/
float Define_determinant(float *matrix, unsigned char size, float tolerance)
{
    float out = 0.0;
    float *temp_matrix = calloc((size * size), sizeof(float));

    if(Gaussian_elimination(matrix, size, size,tolerance, temp_matrix) == size)
    {
        out = 1.0;
        for(unsigned char i = 0U; i < size; i++)
            out *= temp_matrix[i + i * size];
    }

    free(temp_matrix);
    return out;
}

/***********************************************************
* Функция умножения матриц.
* Входные данные:
* first_matrix - Указатель на первую матрицу размером num_out_rows * num_common
* second_matrix - Указатель на вторую матрицу размером num_common * num_out_colums
* num_out_rows - Количество строк первой матрицы
* num_common - Количество строк второй матрицы и количество столбцов первой матрицы
* num_out_colums - Количество столбцов второй матрицы
* Выходные данные:
* out_matrix - Указатель на выходную матрицу размером num_out_rows * num_out_colums
* значение out_matrix не должно равняться ни какому из значений first_matrix или second_matrix
************************************************************/
void Matrices_mul(float *first_matrix, float *second_matrix, unsigned char num_out_rows, unsigned char num_common, unsigned char num_out_colums, float *out_matrix)
{
    for (unsigned char i = 0; i < num_out_rows; i ++)
    {
        for (unsigned char j = 0; j < num_out_colums; j ++)
        {
            // Обнуление значения в выходной матрице
            out_matrix[i * num_out_colums + j] = 0.0;
            for (unsigned char k = 0; k < num_common; k ++)
            {
                out_matrix[i * num_out_colums + j] += first_matrix[i * num_common + k] * second_matrix[num_out_colums * k + j];
            }
        }
    }
}

/***********************************************************
* Функция умножает все значения матрицы на число.
* Входные данные:
* matrix - Указатель на матрицу
* num_rows - Количество строк матрицы
* num_colums - Количество столбцов матрицы
* factor - Множитель матрицы
************************************************************/
void Matrix_mul_by_num(float *matrix, unsigned char num_rows, unsigned char num_colums, float factor)
{
    for (unsigned short i = 0; i < (num_rows * num_colums); i++)
    {
        matrix[i] *=  factor;
    }
}
