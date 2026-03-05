#ifndef _OPTICALFLOW_H_
#define _OPTICALFLOW_H_

#ifdef __cplusplus
extern "C" {
#endif

// Структура входных данных Обычного Лукаса-Канаде[in] 
typedef struct {
    unsigned char* frame_prev;  //Предыдущий кадр
    unsigned char* frame_curr;  //Текущий кадр
    int frame_width;            // Ширина кадра
    int frame_height;           //Высота кадра
    int points_num;             //Количество точек для обработки
    int* keypoints_in;          //Ключевые точки в формате {x0, y0, x1, y1, ...}
    int num_cores;              //Количество ядер процессора(если 0 -> распаралеливание не выполняется)
} OpticalFlowLKIn;

// Структура входных данных Пирамидного Лукаса-Канаде[in]
typedef struct{
    unsigned char* frame_prev;  //Предыдущий кадр
    unsigned char* frame_curr;  //Текущий кадр
    int frame_width;            // Ширина кадра
    int frame_height;           //Высота кадра
    int points_num;             //Количество точек для обработки
    int* keypoints_in;          //Ключевые точки в формате {x0, y0, x1, y1, ...}
    int pyramid_level;          //Количество уровней пирамид Гаусса (не больше 4 иначе работать не будет)
    int img_resize_method;      //1 - метод блочного усреднения, 2 - метод ближайшего соседа
    int num_cores;              //Количество ядер процессора(если 0 -> распаралеливание не выполняется)
} OpticalFlowLKPiramidalIn;

//Выходная структура данных[out]
typedef struct {
    float* keypoints_shift; //Смещение точек в формате {dx0, dy0, dx1, dy1, ...}. Max кол-во = 2 * points_num_in
} OpticalFlowOut;

// Структура для хранения результатов смещения одной точки
typedef struct {
    float dx; //Смещение точки по х
    float dy; //Смещение точки по у
} FlowResult;

/// @brief Функция для вычисления выборочного(заданное количество точек) оптического потока алгоритмом Лукаса-Канаде 
/// @param Input_data Структура входных данных [in]
/// @param Output_data Структура выходных данных[out]
void Lukas_Kanade(const OpticalFlowLKIn* Input_data, OpticalFlowOut* Output_data);

/// @brief Функция для вычисления оптического потока в точке  методом Лукаса-Канаде
/// @param frame1   Предыдущий кадр
/// @param frame2   Текущий кадр
/// @param width    Высота кадра
/// @param height   Ширина кадра
/// @param x        x - координата точки (ширина)
/// @param y        y - координата точки(высота)
/// @return         структура FlowResult[out]
FlowResult Lukas_Kanade_point(unsigned char* frame1, unsigned char* frame2, int width, int height, int x, int y);

/// @brief Функция вычисления оптического потока алгоритмом Пирамидального Лукаса-Канаде
/// @param Input_data   Структура входных данных[in]
/// @param Output_data  Структура выходных данных[out]
/// @return             STATUS - OK или ERROR
STATUS Lukas_Kanade_piramidal(const OpticalFlowLKPiramidalIn* Input_data, OpticalFlowOut* Output_data);

//Структура входных данных алгоритма Farneback
typedef struct {
    unsigned char* frame_prev;  //Предыдущий кадр
    unsigned char* frame_curr;  //Текущий кадр
    int frame_width;            // Ширина кадра
    int frame_height;           //Высота кадра
    int num_cores;              //Количество ядер процессора(если 0 -> распаралеливание не выполняется)
} FarnebackInput;

// It isn`t Farneback method. It is Lukas-Kanade dence method
/// @brief Функция для вычисления плотного оптического потока алгоритмом Фарнебека
/// @param Input_data   Структура входных данных[in]
/// @param flow_x       Указатель на массив смещений по х[out]
/// @param flow_y       Указатель на массив смещений по y[out]
/// @return             STATUS - OK или ERROR
STATUS Farneback(const FarnebackInput* Input_data, float* flow_x, float* flow_y);

#ifdef __cplusplus
}
#endif

#endif
