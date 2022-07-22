#include "poligon.h"
#include "common.h"

/***********************************************************
* Функция определения нахождения точки внутри ПРОСТОГО многоугольника.
* Функция выдает значение "true", если точка с координатами (coordinate_x,coordinate_y)
* находится внутри многоугольника с вершинами в точках с координатами,
* хранящимися в двумерном массиве poligon с размерами amount_of_vertices * 2 (строка {x,y}).
* Направляющий вектор стороны многоугольника  - {x1,y1} --> {x0,y0}.
* Входные данные:
* coordinate_x и coordinate_y - Координаты точки
* poligon - Указатель на массив с координатами вершин многоугольника
* amount_of_vertices - Количество вершин многоугольника
* tolerance - точность вычисления
* Выходные данные:
* out - Нахождение точки (True - внутри многоугольника, False - снаружи многоугольника)
************************************************************/
bool Is_point_in_poligon(float coordinate_x, float coordinate_y, float *poligon, unsigned char amount_of_vertices, float tolerance)
{
    bool out = true;            // Выход функции
    unsigned char count = 0;    // Счетчик пересечений
    float a = 0.0;              // Направляющий вектор стороны многоугольника по оси Х
    float b = 0.0;              // Направляющий вектор стороны многоугольника по оси Y
    float t0 = 0.0;             // Параметр точки пересечения первой прямой
    float t1 = 0.0;             // Параметр точки пересечения второй прямой
    signed char sign = 0;       // Параметр для исключения сингулярности

    // Установка sign для первой точки
    b = poligon[((amount_of_vertices - 1) * 2) + 1] - poligon[1];
    if ( b > 0.0)
    {
        sign = 1;
    }
    else
    {
        sign = -1;
    }

    for (unsigned char i = 0; i < amount_of_vertices; i++ )
    {
        a = poligon[i * 2] - poligon[( ( (i + 1) % amount_of_vertices) ) * 2];
        b = poligon[(i * 2) + 1] - poligon[(( ( (i + 1) % amount_of_vertices) ) * 2) + 1];

        if ( (b * (float) sign) > 0.0 )
        {
            if ( (ABS(poligon[(i * 2) + 1] - coordinate_y) + 0.5 * tolerance) < tolerance )
            {
                if ( ((poligon[i * 2] - coordinate_x + tolerance * ((float) sign)) > 0.0) )
                {
                    count += 1;
                }
            }
        }

        if ( ABS(b) > tolerance)
        {
            if ( b > 0.0)
            {
                sign = 1;
            }
            else
            {
                sign = -1;
            }

            t1 = (coordinate_y - poligon[(( ( (i + 1) % amount_of_vertices) ) * 2) + 1])/b;
            if ( !((t1 > 1.0) || (t1 < 0.0)))
            {
                t0 = poligon[( ( (i + 1) % amount_of_vertices) ) * 2] - coordinate_x + a * t1;
                if ( !((t0 + tolerance * ((float) sign)) < 0.0) )
                {
                    count += 1;
                }
            }
        }
    }

    if ((count % 2) == 0)
    {
        out = false;
    }

    return out;
}
