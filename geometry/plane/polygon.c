#include "polygon.h"
#include "common.h"

/***********************************************************
* Функция определения нахождения точки внутри ПРОСТОГО многоугольника.
* Функция выдает значение "true", если точка с координатами (coordinate_x,coordinate_y)
* находится внутри многоугольника с вершинами в точках с координатами,
* хранящимися в двумерном массиве polygon с размерами amount_of_vertices * 2 (строка {x,y}).
* Направляющий вектор стороны многоугольника  - {x1,y1} --> {x0,y0}.
* Входные данные:
* coordinate_x и coordinate_y - Координаты точки
* polygon - Указатель на массив с координатами вершин многоугольника
* amount_of_vertices - Количество вершин многоугольника
* tolerance - точность вычисления
* Выходные данные:
* out - Нахождение точки (True - внутри многоугольника, False - снаружи многоугольника)
************************************************************/
bool Is_point_in_polygon(float coordinate_x, float coordinate_y, float *polygon, unsigned char amount_of_vertices, float tolerance)
{
    bool out = true;            // Выход функции
    unsigned char count = 0;    // Счетчик пересечений
    float a = 0.0;              // Направляющий вектор стороны многоугольника по оси Х
    float b = 0.0;              // Направляющий вектор стороны многоугольника по оси Y
    float t0 = 0.0;             // Параметр точки пересечения первой прямой
    float t1 = 0.0;             // Параметр точки пересечения второй прямой
    signed char sign = 0;       // Параметр для исключения сингулярности

    // Установка sign для первой точки
    b = polygon[((amount_of_vertices - 1) * 2) + 1] - polygon[1];
    if ( b > 0.0) 
        sign = 1;
    else 
        sign = -1;

    for (unsigned char i = 0; i < amount_of_vertices; i++ )
    {
        a = polygon[i * 2] - polygon[( ( (i + 1) % amount_of_vertices) ) * 2];
        b = polygon[(i * 2) + 1] - polygon[(( ( (i + 1) % amount_of_vertices) ) * 2) + 1];

        if ( (b * (float) sign) > 0.0 )
        {
            if ( (ABS(polygon[(i * 2) + 1] - coordinate_y) + 0.5 * tolerance) < tolerance )
            {
                if ( ((polygon[i * 2] - coordinate_x + tolerance * ((float) sign)) > 0.0) ) 
                    count += 1;
            }
        }

        if ( ABS(b) > tolerance)
        {
            if ( b > 0.0) 
                sign = 1;
            else 
                sign = -1;

            t1 = (coordinate_y - polygon[(( ( (i + 1) % amount_of_vertices) ) * 2) + 1])/b;
            if ( !((t1 > 1.0) || (t1 < 0.0)))
            {
                t0 = polygon[( ( (i + 1) % amount_of_vertices) ) * 2] - coordinate_x + a * t1;

                if ( !((t0 + tolerance * ((float) sign)) < 0.0) ) 
                    count += 1;
            }
        }
    }

    if ((count % 2) == 0) 
        out = false;

    return out;
}

/***********************************************************
* Функция вычисления площади ПРОСТОГО многоугольника
* Входные данные:
* polygon - Указатель на массив с координатами вершин многоугольника
* хранящимися в двумерном массиве polygon с размерами amount_of_vertices * 2 (строка {x,y}).
* amount_of_vertices - Количество вершин многоугольника
* Выходные данные:
* out - Площадь многоугольника 
* (out > 0 - расположение вершин в массиве против ЧС; out < 0 - расположение вершин по ЧС )
************************************************************/
float Calculate_polygon_area(float *polygon, unsigned char amount_of_vertices)
{
    float out = 0.0;

    for (unsigned char i = 0; i < amount_of_vertices; i++)
    {
        out += polygon[i * 2] * polygon[(( ( (i + 1) % amount_of_vertices) ) * 2) + 1];
        out -= polygon[(i * 2) + 1] * polygon[( ( (i + 1) % amount_of_vertices) ) * 2];
    }

    return (out / 2.0);
}
