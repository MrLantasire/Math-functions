#include <stdio.h>
#include "common.h"
#include "matrix.h"
#include "complex.h"
#include "equation.h"
#include "gis.h"
#include "polygon.h"
#include "plane_curve_approximation.h"
#include "space_curve_approximation.h"
#include "vector.h"
#include "quaternion.h"

static void Show_matrices(void);
static void Show_complex_numbers(void);
static void Show_equations(void);
static void Show_gis(void);
static void Show_polygon(void);
static void Show_plane_curve_approximation(void);
static void Show_space_curve_approximation(void);
static void Show_vector(void);
static void Show_quaternion(void);

int main()
{   
    Show_matrices();
    Show_complex_numbers();
    Show_equations();
    Show_gis();
    Show_polygon();
    Show_plane_curve_approximation();
    Show_space_curve_approximation();
    Show_vector();
    Show_quaternion();

    return 0;
}

// Примеры работы с функциями модуля "matrix"
static void Show_matrices(void)
{
    
    // Пример квадратной матрицы
    static float instance_sq_matrix[4][4] = {   {5.0,4.0,3.0,2.0},
                                                {1.0,2.0,2.0,5.0},
                                                {1.0,5.0,1.0,5.0},
                                                {2.0,1.0,2.0,1.0}};

    // Пример обычной матрицы
    static float instance_matrix[4][5] = {      {5.0,4.0,3.0,2.0,30.0},
                                                {1.0,3.0,3.0,6.0,40.0},
                                                {1.0,1.0,1.0,1.0,10.0},
                                                {0.0,4.0,5.0,-2.0,15.0}};

    // Пример треугольной матрицы, полученной после Гауссова исключения
    static float instanse_tr_matrix[4][5] = {0};

    printf("\n");

    // Определитель квадратной матрицы
    printf("Determinant of matrix 4 x 4 = %f\n", Define_determinant(&instance_sq_matrix[0][0], 4U,1.E-6));
    // Ранг матрицы и преобразование матрицы к треугольному виду
    printf("Rank of matrix 4 x 5 = %d\n", Gaussian_elimination(&instance_matrix[0][0], 4U, 5U, 1.E-6, &instanse_tr_matrix[0][0]));

    printf("\nTriangle matrix 4 x 5\n");
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            printf("%f\t", instanse_tr_matrix[i][j]);
        }
        printf("\n");
    }
    
    // Умножение матрицы на число
    Matrix_mul_by_num(&instance_sq_matrix[0][0], 4U, 4U, 2.0);
    printf("\nSquare matrix * 2\n");
    printf("Determinant = %f\n", Define_determinant(&instance_sq_matrix[0][0], 4U,1.E-6));
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            printf("%f\t", instance_sq_matrix[i][j]);
        }
        printf("\n");
    }

    // Умножение матриц
    printf("\n( matrix 4 x 4 ) X ( matrix 4 x 5 )\n");
    Matrices_mul(&instance_sq_matrix[0][0], &instance_matrix[0][0], 4U, 4U, 5U, &instanse_tr_matrix[0][0]);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            printf("%f\t", instanse_tr_matrix[i][j]);
        }
        printf("\n");
    }

    printf("\n");
}

// Примеры работы с функциями модуля "complex"
static void Show_complex_numbers(void)
{
    // Комплексные числа
    static complex32_t a = {1.0, 0.0};
    static complex32_t b = {0.0, 1.0};
    static complex32_t c = {0};

    printf("\n");

    // Аргумент и модуль комплексных чисел
    printf("arg(a) = %f ; abs(a) = %f\n", Complex_arg(a), Complex_abs(a));
    printf("arg(b) = %f ; abs(b) = %f\n", Complex_arg(b), Complex_abs(b));

    // Сложение комплексных чисел
    c = Complex_add(a,b);
    printf("a + b = %f%+fi\n", c.real, c.imag);

    // Вычитание комплексных чисел
    c = Complex_sub(a,b);
    printf("a - b = %f%+fi\n", c.real, c.imag);

    // Умножение комплексных чисел
    c = Complex_mul(a,b);
    printf("a * b = %f%+fi\n", c.real, c.imag);

    // Деление комплексных чисел
    c = Complex_div(a,b);
    printf("a / b = %f%+fi\n", c.real, c.imag);

    // Возведение в степень комплексного числа
    c = Complex_pow(Complex_add(a,b) , 2);
    printf("(a + b) ^ 2 = %f%+fi\n", c.real, c.imag);

    // Извлечение корня комплексного числа
    printf("(a - b) ^ (1/4) = ");
    for (int i = 0; i < 4; i++)
    {
        c = Complex_root(Complex_sub(a,b), 4, i);
        printf("%f%+fi ", c.real, c.imag);
    }
    printf("\n");

    // Получение комплексно-сопряженного числа
    c = Complex_conjugate(b);
    printf("conjugate(b)= %f%+fi\n", c.real, c.imag);

    // Сравнение комплексных чисел
    c = Complex_mul(a,b);
    a = Complex_div(a,b);
    printf("a * b == a / b = ");
    if ( Is_complex_equal(a, c, 1.E-6) )
        printf("True\n");
    else
        printf("False\n");

    a = Complex_conjugate(a);
    printf("a * b == conjugate(a / b) = ");
    if ( Is_complex_equal(a, c, 1.E-6) )
        printf("True\n");
    else
        printf("False\n");

    printf("\n");
}

// Примеры работы с функциями модуля "equation"
static void Show_equations(void)
{
    // Массив решений уравнений
    static complex32_t solutions[3] = {0};

    printf("\n");

    // Решение квадратных уравнений
    // Дискриминанты > 0
    printf("Equation x^2 + 3x - 4 = 0\n");
    printf("Discriminant = %f\n", Solve_quadratic_equation(1.0, 3.0, -4.0, &solutions[0]));

    for (int i = 0; i < 2; i++)
    {
        printf("x%d =  %f%+fi\n", i+1, solutions[i].real,solutions[i].imag);
    }   

    // Дискриминанты = 0
    printf("\nEquation 2x^2 - 12x + 18 = 0\n");
    printf("Discriminant = %f\n", Solve_quadratic_equation(2.0, -12.0, 18.0, &solutions[0]));

    for (int i = 0; i < 2; i++)
    {
        printf("x%d =  %f%+fi\n", i+1, solutions[i].real,solutions[i].imag);
    }   

    // Дискриминанты < 0
    printf("\nEquation x^2 + 4 = 0\n");
    printf("Discriminant = %f\n", Solve_quadratic_equation(1.0, 0.0, 4.0, &solutions[0]));

    for (int i = 0; i < 2; i++)
    {
        printf("x%d =  %f%+fi\n", i+1, solutions[i].real,solutions[i].imag);
    }   

    // Решение кубических уравнений
    // Дискриминанты > 0
    printf("\nEquation x^3 - 6x^2 + 11x - 6 = 0\n");
    printf("Discriminant = %f\n", Solve_cubic_equation(1.0, -6.0, 11.0, -6.0, &solutions[0]));

    for (int i = 0; i < 3; i++)
    {
        printf("x%d =  %f%+fi\n", i+1, solutions[i].real,solutions[i].imag);
    }   

    // Дискриминанты = 0
    printf("\nEquation -2x^3 + 18x^2 - 48x + 32 = 0\n");
    printf("Discriminant = %f\n", Solve_cubic_equation(-2.0, 18.0, -48.0, 32.0, &solutions[0]));

    for (int i = 0; i < 3; i++)
    {
        printf("x%d =  %f%+fi\n", i+1, solutions[i].real,solutions[i].imag);
    }   

    // Дискриминанты < 0
    printf("\nEquation -x^3 + 3x^2 - 9x + 27 = 0\n");
    printf("Discriminant = %f\n", Solve_cubic_equation(-1.0, 3.0, -9.0, 27.0, &solutions[0]));

    for (int i = 0; i < 3; i++)
    {
        printf("x%d =  %f%+fi\n", i+1, solutions[i].real,solutions[i].imag);
    }   

    printf("\n");
}

// Примеры работы с функциями модуля "gis"
static void Show_gis(void)
{
    // Координаты Москвы, [градусы]
    float lat_Mos = 55.755797;
    float long_Mos = 37.617729;
    // Координаты Санкт-Петербурга, [градусы]
    float lat_SPb = 59.939037;
    float long_SPb = 30.315828;
    // Координаты Сочи, [градусы]
    float lat_Soc = 43.585326;
    float long_Soc = 39.722453;
    // Координаты Владивосток, [градусы]
    float lat_Vlk = 43.115295;
    float long_Vlk = 131.885453;

    printf("\n");
    
    // Вычисление расстояний
    printf("Distance: Moscow - Saint-Peterburg = %f [m]\n", Get_distance_from_geopoints(lat_Mos, long_Mos, lat_SPb, long_SPb, EARTH_RADIUS));
    printf("Distance: Sochi - Vladivostok = %f [m]\n", Get_distance_from_geopoints(lat_Soc, long_Soc, lat_Vlk, long_Vlk, EARTH_RADIUS));

    // Вычисление направлений
    printf("Course: Saint-Peterburg - Sochi = %f [grad]\n", Get_direction_from_geopoints(lat_SPb, long_SPb, lat_Soc, long_Soc));
    printf("Course: Vladivostok - Moscow = %f [grad]\n", Get_direction_from_geopoints(lat_Vlk, long_Vlk, lat_Mos, long_Mos));

    // Вычисление координат по направлению и расстоянию
    Get_geopoint_to_direction(lat_Soc, long_Soc, -5.591860, 1360820.75, EARTH_RADIUS, &lat_Mos, &long_Mos);
    Get_geopoint_to_direction(lat_SPb, long_SPb, 56.738583, 6533325.5, EARTH_RADIUS, &lat_Vlk, &long_Vlk);
    printf("Coordinates Moscow: N%f E%f\n", lat_Mos, long_Mos);
    printf("Coordinates Vladivostok: N%f E%f\n", lat_Vlk, long_Vlk);

    printf("\n");
}

// Примеры работы с функциями модуля "polygon"
static void Show_polygon(void)
{
    // Тестовый многоугольник
    static float test_polygon[6][2] = { {40.0   , 10.0  },
                                        {25.0   , 35.98 },
                                        {-5     , 35.98 },
                                        {-20    , 10.0  },
                                        {-5     , -15.98},
                                        {25.0   , -15.98}};

    // Координаты тестовой точки
    float x_of_point = 10.0;
    float y_of_point = 10.0;

    printf("\n");

    // Определение положения точки относительно многоугольника
    for (int i = 0; i < 10; i++)
    {
        if (Is_point_in_polygon(x_of_point, y_of_point, &test_polygon[0][0], 6U, 1.E-6))
        {
            printf("P%d\t(%f, %f)\tInside poligon.\n", i+1, x_of_point, y_of_point);
        }
        else
        {
            printf("P%d\t(%f, %f)\tOutside poligon.\n", i+1, x_of_point, y_of_point);
        }
        x_of_point += 4.0;
        y_of_point += 4.0;
    }

    // Вычисление площади многоугольника
    printf("\nArea of polygon: S = %f\n", Calculate_polygon_area(&test_polygon[0][0], 6U));

    // Разворот многоугольника
    for (int i = 0; i < 3; i++)
    {
        x_of_point = test_polygon[i][0];
        y_of_point = test_polygon[i][1];
        test_polygon[i][0] = test_polygon[5 - i][0];
        test_polygon[i][1] = test_polygon[5 - i][1];
        test_polygon[5 - i][0] = x_of_point;
        test_polygon[5 - i][1] = y_of_point;
    }
    // Вычисление площади развернутого многоугольника
    printf("\nArea of reversed polygon: S = %f\n", Calculate_polygon_area(&test_polygon[0][0], 6U));

    printf("\n");
}

// Примеры работы с функциями модуля "plane_curve_approximation"
static void Show_plane_curve_approximation(void)
{
    // Набор точек (вершины правильного многоугольника)
    static float set_of_points[6][2] = {{40.0   , 10.0  },
                                        {25.0   , 35.98 },
                                        {-5     , 35.98 },
                                        {-20    , 10.0  },
                                        {-5     , -15.98},
                                        {25.0   , -15.98}};
    
    // Переменные для хранения результатов
    // Координаты точки
    float x_of_point = 0.0;
    float y_of_point = 0.0;
    // Координаты вектора
    float x_of_vector = 0.0;
    float y_of_vector = 0.0;

    printf("\n");

    // Аппроксимация точек окружностью (первый набор)
    printf("Circle 1: Radius = %f, ", Circle_approximation(&set_of_points[0][0], 6U, &x_of_point, &y_of_point));
    printf("Center = (%f, %f)\n", x_of_point, y_of_point);

    // Аппроксимация точек прямой (первый набор)
    printf("Line 1: Vector length = %f, ", Line_approximation(&set_of_points[0][0], 6U, &x_of_point, &y_of_point, &x_of_vector, &y_of_vector));
    printf("Starting point = (%f, %f), Vector = (%f, %f)\n", x_of_point, y_of_point, x_of_vector, y_of_vector);

    // Сжатие точек к оси 0Х
    set_of_points[1][1] = 25.0;
    set_of_points[2][1] = 25.0;
    set_of_points[4][1] = -5.0;
    set_of_points[5][1] = -5.0;

    // Аппроксимация точек окружностью (второй набор)
    printf("\nCircle 2: Radius = %f, ", Circle_approximation(&set_of_points[0][0], 6U, &x_of_point, &y_of_point));
    printf("Center = (%f, %f)\n", x_of_point, y_of_point);

    // Аппроксимация точек прямой (второй набор)
    printf("Line 2: Vector length = %f, ", Line_approximation(&set_of_points[0][0], 6U, &x_of_point, &y_of_point, &x_of_vector, &y_of_vector));
    printf("Starting point = (%f, %f), Vector = (%f, %f)\n", x_of_point, y_of_point, x_of_vector, y_of_vector);

    printf("\n");
}

// Примеры работы с функциями модуля "space_curve_approximation"
static void Show_space_curve_approximation(void)
{
    // Набор точек (вершины правильного многоугольника)
    static float set_of_points[36][3] = {{ 10.0  , 10.0  , 10.0},
                                        { 20.0  , 10.0  , 10.0},
                                        { 20.0  , 20.0  , 10.0},
                                        { 10.0  , 20.0  , 10.0},
                                        { 10.0  , 10.0  , 20.0},
                                        { 20.0  , 10.0  , 20.0},
                                        { 20.0  , 20.0  , 20.0},
                                        { 10.0  , 20.0  , 20.0}};

    // Переменные для хранения результатов
    static float curve_parameters[2][3] = {0};

    printf("\n");

    // Аппроксимация точек прямой (первый набор)
    printf("Line_3D 1: Vector length = %f, ", Line_approximation_3D(&set_of_points[0][0], 8U, &curve_parameters[0][0]));
    printf("Starting point = (%f, %f, %f), ", curve_parameters[0][0], curve_parameters[0][1], curve_parameters[0][2]);
    printf("Vector = (%f, %f, %f)\n", curve_parameters[1][0], curve_parameters[1][1], curve_parameters[1][2]);

    // Изменение координат точек
    for (int i = 0; i < 36; i++)
    {
        set_of_points[i][0] = 2.0 * sin((double)(M_PI * i )/1.0) + 5.0 + (float)(1.0 * i); 
        set_of_points[i][1] = 2.0 * cos((double)(M_PI * i )/1.0) + 5.0 + (float)(1.0 * i); 
        set_of_points[i][2] = 10.0 + (float)(1.0 * i); 
    }

    // Аппроксимация точек прямой (второй набор)
    printf("Line_3D 2: Vector length = %f, ", Line_approximation_3D(&set_of_points[0][0], 36U, &curve_parameters[0][0]));
    printf("Starting point = (%f, %f, %f), ", curve_parameters[0][0], curve_parameters[0][1], curve_parameters[0][2]);
    printf("Vector = (%f, %f, %f)\n", curve_parameters[1][0], curve_parameters[1][1], curve_parameters[1][2]);

    printf("\n");
}

// Примеры работы с функциями модуля "vector"
static void Show_vector(void)
{
    // Векторы
    vector32_3D_t a = {1.0,1.0,0.0};
    vector32_3D_t b = {-1.0,1.0,0.0};
    vector32_3D_t c = {0.0,0.0,1.0};
    vector32_3D_t d = {0.0,0.0,0.0};

    // Массив для матрицы поворота
    static float matrix[3][3] = {0};

    printf("\n");

    
    // Длина векторов
    printf("a{%f, %f, %f} |a| = %f\n", a.i, a.j, a.k, Vector_abs(a));
    printf("b{%f, %f, %f} |b| = %f\n", b.i, b.j, b.k, Vector_abs(b));
    printf("c{%f, %f, %f} |c| = %f\n", c.i, c.j, c.k, Vector_abs(c));

    // Сложение векторов
    d = Vector_add(a,b);
    printf("a + b{%f, %f, %f} |a + b| = %f\n", d.i, d.j, d.k, Vector_abs(d));
    // Вычитание векторов
    d = Vector_sub(c,a);
    printf("c - a{%f, %f, %f} |c - a| = %f\n", d.i, d.j, d.k, Vector_abs(d));
    // Умножение вектора на число
    d = Vector_mul(c, 15.0);
    printf("15 * c = {%f, %f, %f} |15 * c| = %f\n", d.i, d.j, d.k, Vector_abs(d));
    // Векторное произведение векторов
    d = Vector_cross_mul(c,b);
    printf("c X b{%f, %f, %f} |c X b| = %f\n", d.i, d.j, d.k, Vector_abs(d));
    // Скалярное произведение векторов
    printf("d{%f, %f, %f} |d| = %f a * d = %f\n", d.i, d.j, d.k, Vector_abs(d), Vector_scalar_mul(a,d));
    // Смешанное произведение векторов
    printf("a * (b X c) = %f\n", Vector_mixed_mul(a,b,c));
    // Вычисление нормализованного вектора
    d = Vector_normalize(d);
    printf("d / |d| {%f, %f, %f} |d / |d|| = %f\n", d.i, d.j, d.k, Vector_abs(d));
    // Вычисление угла между векторами
    printf("b^c = %f [grad]\t", RAD_TO_GRAD(Define_angle_of_vectors(b,c)));
    printf("a^d = %f [grad]\n", RAD_TO_GRAD(Define_angle_of_vectors(a,d)));

    // Матрицы поворота
    Get_rotation_matrix((vector32_3D_t) {1.0, 0.0, 0.0}, GRAD_TO_RAD(-45.0), &matrix[0][0]);
    printf("\nX-axis rotation matrix by -45 [grad]\n");
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }

    Get_rotation_matrix((vector32_3D_t){0.0, 1.0, 0.0}, GRAD_TO_RAD(90.0), &matrix[0][0]);
    printf("\nY-axis rotation matrix by 90 [grad]\n");
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }

    Get_rotation_matrix((vector32_3D_t){0.0, 0.0, 1.0}, GRAD_TO_RAD(180.0), &matrix[0][0]);
    printf("\nZ-axis rotation matrix by 180 [grad]\n");
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }

    d = Vector_cross_mul(a,b);
    Get_rotation_matrix(d, Define_angle_of_vectors(a,b), &matrix[0][0]);
    printf("\nFree-axis {%f, %f, %f} rotation matrix by %f [grad]\n", d.i, d.j, d.k, RAD_TO_GRAD(Define_angle_of_vectors(a,b)));
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }

    // Сравнение векторов
    c = Vector_sub(a,b);
    d = Vector_sub(b,a);
    printf("\nc{%f, %f, %f} == d{%f, %f, %f} = ", c.i, c.j, c.k, d.i, d.j, d.k);
    if ( Is_vector_equal(c, d, 1.E-6) )
        printf("True\n");
    else
        printf("False\n");
    
    d = Vector_mul(d, -1.0);
    printf("c{%f, %f, %f} == d{%f, %f, %f} = ", c.i, c.j, c.k, d.i, d.j, d.k);
    if ( Is_vector_equal(c, d, 1.E-6) )
        printf("True\n");
    else
        printf("False\n");

    printf("\n");
}

// Примеры работы с функциями модуля "quaternon"
static void Show_quaternion(void)
{
    // Кватернионы
    quaternion32_t a = {1.0, {2.0, 3.0, 4.0}};
    quaternion32_t b = {1.0, {0}};
    quaternion32_t c = {0};
    quaternion32_t d = {0};

    // Вектор для вращения
    vector32_3D_t axis = {0.0, 0.0, 1.0};

    printf("\n");

    printf("a{%f, %fi, %fj, %fk} |a| = %f\n", a.scalar, a.vector.i, a.vector.j, a.vector.k, Quaternion_abs(a));
    printf("b{%f, %fi, %fj, %fk} |b| = %f\n", b.scalar, b.vector.i, b.vector.j, b.vector.k, Quaternion_abs(b));

    // Сложение кватернионов
    c = Quaternion_add(a,b);
    printf("a + b {%f, %fi, %fj, %fk} |a + b| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    // Вычитание кватернионов
    c = Quaternion_sub(a,b);
    printf("a - b {%f, %fi, %fj, %fk} |a - b| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    // Умножение кватернионов
    c = Quaternion_mul(a,a);
    printf("a * a {%f, %fi, %fj, %fk} |a * a| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    // Деление кватернионов
    c = Quaternion_div(b,a);
    printf("b / a {%f, %fi, %fj, %fk} |b / a| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    // Обратный кватернион
    d = Quaternion_reciprocal(a);
    printf("a ^ (-1) {%f, %fi, %fj, %fk} |a ^ (-1)| = %f\n", d.scalar, d.vector.i, d.vector.j, d.vector.k, Quaternion_abs(d));
    // Сравнивание кватернионов
    printf("b / a == a ^ (-1) = ");
    if ( Is_quaternion_equal(c, d, 1.E-6) )
        printf("True\n");
    else
        printf("False\n");
    // Сопряженный кватернион
    c = Quaternion_conjugate(c);
    printf("conjugate(b / a) {%f, %fi, %fj, %fk} |conjugate(b / a)| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    printf("conjugate(b / a) == a ^ (-1) = ");
    if ( Is_quaternion_equal(c, d, 1.E-6) )
        printf("True\n");
    else
        printf("False\n");
    // Нормализация кватерниона
    c = Quaternion_normalize(d);
    printf("a ^ (-1) / |a ^ (-1)| {%f, %fi, %fj, %fk} |a ^ (-1) / |a ^ (-1)|| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    // Кватернион поворота
    c = Quaternion_rotor(axis, GRAD_TO_RAD(90.0));
    printf("rotor{%f, %fi, %fj, %fk} |rotor| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));
    // Поворот кватерниона
    c = Rotate_quaternion(a, axis, GRAD_TO_RAD(90.0));
    printf("rot(a){%f, %fi, %fj, %fk} |rot(a)| = %f\n", c.scalar, c.vector.i, c.vector.j, c.vector.k, Quaternion_abs(c));


    printf("\n");
}
