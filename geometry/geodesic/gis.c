#include "gis.h"
#include "common.h"

/***********************************************************
* Функция вычисления расстояния между точками, заданными геодезическими координатами
* Входные данные:
* first_latitude и first_longitude - Геодезические координаты широты и долготы первой точки (градусы)
* second_latitude и second_longitude - Геодезические координаты широты и долготы второй точки (градусы)
* radius - радиус сфероида (требуемые динейные единицы измерения: м, км и т.д.)
* Выходные данные:
* Расстояние между точками (требуемые динейные единицы измерения: м, км и т.д.)
************************************************************/
float Get_distance_from_geopoints(float first_latitude, float first_longitude, float second_latitude, float second_longitude, float radius)
{
    double argument = 0.0;

    argument = sin( GRAD_TO_RAD(second_latitude) ) * sin( GRAD_TO_RAD(first_latitude) ) 
    + cos( GRAD_TO_RAD(second_latitude) ) * cos( GRAD_TO_RAD(first_latitude) ) * cos( GRAD_TO_RAD(first_longitude - second_longitude) );
    
    return (radius * acos(argument));
}

/***********************************************************
* Функция вычисления азимута из первой точки на вторую, заданными геодезическими координатами
* Входные данные:
* first_latitude и first_longitude - Геодезические координаты широты и долготы первой точки (градусы)
* second_latitude и second_longitude - Геодезические координаты широты и долготы второй точки (градусы)
* Выходные данные:
* Азимут (градусы)
************************************************************/
float Get_direction_from_geopoints(float first_latitude, float first_longitude, float second_latitude, float second_longitude)
{
    double numenator = 0.0;
    double denominator = 0.0;

    numenator = sin(GRAD_TO_RAD(second_longitude - first_longitude)) * cos( GRAD_TO_RAD(second_latitude));

    denominator = cos( GRAD_TO_RAD(first_latitude) ) * sin( GRAD_TO_RAD(second_latitude) );
    denominator -= sin( GRAD_TO_RAD(first_latitude) ) * cos( GRAD_TO_RAD(second_latitude) ) * cos( GRAD_TO_RAD(second_longitude - first_longitude) );

    return ((float) RAD_TO_GRAD(atan2(numenator, denominator)));
}

/***********************************************************
* Функция вычисления координат второй точки относительно первой по направлению и расстоянию
* Входные данные:
* latitude и longitude - Геодезические координаты широты и долготы первой точки (градусы)
* direction - Азимут (градусы)
* distance - Расстояние до второй точки (требуемые динейные единицы измерения: м, км и т.д.)
* radius - радиус сфероида (требуемые динейные единицы измерения: м, км и т.д.)
* Выходные данные:
* out_latitude и out_longitude - Указатели для записи вычисленных значений широты и долготы (градусы)
************************************************************/
void Get_geopoint_to_direction(float latitude, float longitude, float direction, float distance, float radius, float *out_latitude, float *out_longitude)
{
    double argument = 0.0;

    argument = ( sin( GRAD_TO_RAD(latitude) ) * cos( distance / radius ) 
    + cos( GRAD_TO_RAD(latitude) ) * sin( distance / radius ) * cos( GRAD_TO_RAD(direction) ) );
    
    *out_latitude = RAD_TO_GRAD( asin(argument) );

    argument = ( cos( GRAD_TO_RAD(latitude) ) * cos( distance / radius ) 
    - sin( GRAD_TO_RAD(latitude) ) * sin( distance / radius ) * cos( GRAD_TO_RAD(direction) ) );

    *out_longitude = longitude + RAD_TO_GRAD( atan2( (sin( distance / radius ) * sin( GRAD_TO_RAD(direction) )), argument ) ) ;
}
