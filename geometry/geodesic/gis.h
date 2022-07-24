#ifndef __GIS_H
#define __GIS_H

#define EARTH_RADIUS 6367449.1458234153093  // Радиус Земли в системе WGS-84 (м)

// Функция вычисления расстояния между точками, заданными геодезическими координатами
extern float Get_distance_from_geopoints(float first_latitude, float first_longitude, float second_latitude, float second_longitude, float radius);

// Функция вычисления азимута из первой точки на вторую, заданными геодезическими координатами
extern float Get_direction_from_geopoints(float first_latitude, float first_longitude, float second_latitude, float second_longitude);

//Функция вычисления координат второй точки относительно первой по направлению и расстоянию
extern void Get_geopoint_to_direction(float latitude, float longitude, float direction, float distance, float radius, float *out_latitude, float *out_longitude);

#endif //__GIS_H
