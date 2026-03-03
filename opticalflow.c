#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "opticalflow.h"
#include <multicore.h> //Оптимизация для распаралеливания(будующая)
#include "img_size_reducion.h" //Построение уменьшеных изображений для пирамид

#define WINDOW_SIZE 1 //радиус окружности для расчета(1 -> 3х3)
#define Koef_K 2 //Коэффициент сжатия изображения
#define Max_Pyramid_Level 4 //Максимальноек колличество уровней пирамиды Гаусса
#define PYRAMID_LEVELS 3 //Количество уровней пирамид для плотноного оптического потока методом SAD

// Веса окна по Гаусианне(нормализованные, ср.кв.откл=0.8)
static const float GaussianWeights[3][3] = {
    {0.061f, 0.135f, 0.061f},
    {0.135f, 0.434f, 0.135f},
    {0.061f, 0.135f, 0.061f}
};

// Структура для трехмерного массива градиента
typedef struct {
    float data[3][3][3];
} Gradient3x3x3;

// Структура для пирамид Гаусса(SAD)
typedef struct {
    unsigned char* data;  // Выделяется нами
    int width;
    int height;
} ImagePyramidLevel;

// Подсчет значений частных производных(градиента)
static Gradient3x3x3 calculate_gradient(unsigned char* frame1, unsigned char* frame2, 
                                        int width, int height, int x, int y) {
    Gradient3x3x3 grad;
    float Ix[3][3], Iy[3][3], It[3][3];
    
    // Проверка, что окно не выходит за границы
    if (x < WINDOW_SIZE || x >= width - WINDOW_SIZE ||
        y < WINDOW_SIZE || y >= height - WINDOW_SIZE) {
        // Заполняем нулями
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    grad.data[k][i][j] = 0.0f;
        return grad;
    }
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int di = i - 1;
            int dj = j - 1;
            int px = x + dj;
            int py = y + di;
            
            // Индексы с проверкой границ
            int px_left = (px > 0) ? px - 1 : 0;
            int px_right = (px < width - 1) ? px + 1 : width - 1;
            int py_up = (py > 0) ? py - 1 : 0;
            int py_down = (py < height - 1) ? py + 1 : height - 1;
            
            Ix[i][j] = (frame1[py * width + px_right] - frame1[py * width + px_left]) / 2.0f;
            Iy[i][j] = (frame1[py_down * width + px] - frame1[py_up * width + px]) / 2.0f;
            It[i][j] = (float)(frame2[py * width + px] - frame1[py * width + px]) / 255.0f;
        }
    }
    
    // Копирование в структуру
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            grad.data[0][i][j] = Ix[i][j];
            grad.data[1][i][j] = Iy[i][j];
            grad.data[2][i][j] = It[i][j];
        }
    }
    return grad;
}

// Алгоритм лукаса-канаде для одной точки
FlowResult Lukas_Kanade_point(unsigned char* frame1, unsigned char* frame2, int width, int height, int x, int y) {
    FlowResult result = {0, 0}; // Инициализация структуры результата
    Gradient3x3x3 grad = calculate_gradient(frame1, frame2, width, height, x, y); //Вычисление градиента
    float (*Ix)[3] = grad.data[0];
    float (*Iy)[3] = grad.data[1]; 
    float (*It)[3] = grad.data[2];
    float det = 0, Ix2 = 0, Iy2 = 0, IxIy = 0, IxIt = 0, IyIt = 0;

    for (int i = 0; i < WINDOW_SIZE*2+1; i++) {
        for (int j = 0; j < WINDOW_SIZE*2+1; j++) {
            const float w = GaussianWeights[i][j];//Инициализация веса
            Ix2 += Ix[i][j] * Ix[i][j] * w;
            Iy2 += Iy[i][j] * Iy[i][j] * w;
            IxIy += Ix[i][j] * Iy[i][j] * w;
            IxIt += Ix[i][j] * It[i][j] * w;
            IyIt += Iy[i][j] * It[i][j] * w;
            }
        }

    det = Ix2 * Iy2 - IxIy * IxIy;
    if (fabs(det) < 1e-6) {
        return result; //Система не имеет решений либо беск.много решений
    }

    //Вычисление результатом методом крамера
    result.dx = - (Iy2*IxIt -IxIy * IyIt) / det;
    result.dy = (IxIy*IxIt - Ix2*IyIt) / det;

    return result;
}

//Заглавная функция для обработки массы контрольных точек
void Lukas_Kanade(const OpticalFlowLKIn* Input_data, OpticalFlowOut* Output_data){
    //Инициализация входных данных
    const int widht = Input_data->frame_width;
    const int height = Input_data->frame_height;
    const int point_num = Input_data->points_num;

    int* points_coords = Input_data->keypoints_in;
    unsigned char* frame1 = Input_data->frame_prev;
    unsigned char* frame2 = Input_data->frame_curr;

    //Инициализируем циклом алгоритм Лукаса-Канаде для каждой ключевой точки
    for (int i = 0; i < point_num; i++){
        const int x = points_coords[i*2];
        const int y = points_coords[i*2+1];
        FlowResult LKresult = Lukas_Kanade_point(frame1, frame2, widht, height, x, y);
        //Возвращаем полученные смещения в структуру вывода
        Output_data->keypoints_shift[i*2] = LKresult.dx;
        Output_data->keypoints_shift[i*2+1] = LKresult.dy;
    }
}

// Функция для очистки памяти пирамид
static void pyramidMemoryDelete(unsigned char* pyramid1[Max_Pyramid_Level], unsigned char* pyramid2[Max_Pyramid_Level], int levels) {
    for (int level = 1; level < levels; level++) {
        if (pyramid1[level] != NULL) {
            free(pyramid1[level]);
            pyramid1[level] = NULL;
        }
        if (pyramid2[level] != NULL) {
            free(pyramid2[level]);
            pyramid2[level] = NULL;
        }
    }
}

STATUS Lukas_Kanade_piramidal(const OpticalFlowLKPiramidalIn* Input_data, OpticalFlowOut* Output_data){
    //Инициализируем входные данные
    const int width = Input_data->frame_width;
    const int height = Input_data->frame_height;
    const int point_num = Input_data->points_num;
    const int pyramid_lvl = Input_data->pyramid_level;
    const int resize_method = Input_data->img_resize_method;
    
    int* points_coords = Input_data->keypoints_in;
    unsigned char* frame1 = Input_data->frame_prev;
    unsigned char* frame2 = Input_data->frame_curr;
    
    // Массивы для хранения пирамид изображений
    unsigned char* pyramid1[Max_Pyramid_Level] = {NULL, NULL, NULL, NULL};
    unsigned char* pyramid2[Max_Pyramid_Level] = {NULL, NULL, NULL, NULL};
    int pyramid_width[Max_Pyramid_Level], pyramid_height[Max_Pyramid_Level];
    
    // Уровень 0 - оригинальные изображения
    pyramid1[0] = frame1;
    pyramid2[0] = frame2;
    pyramid_width[0] = width;
    pyramid_height[0] = height;
    
    // Создаем пирамиду для frame1(предыдущего кадра)
    for (int level = 1; level < pyramid_lvl; level++) {
        InputImage input;
        OutputImage output;
        
        input.data = pyramid1[level - 1];
        input.width = pyramid_width[level - 1];
        input.height = pyramid_height[level - 1];
        
        output.width = input.width / Koef_K;
        output.height = input.height / Koef_K;
        
        pyramid_width[level] = output.width;
        pyramid_height[level] = output.height;
        
        output.data = (unsigned char*)malloc(output.width * output.height * sizeof(unsigned char));
        if (output.data == NULL){
            pyramidMemoryDelete(pyramid1, pyramid2, pyramid_lvl);
            return ERROR;
        }
        pyramid1[level] = output.data;
        
        //Разделение на методы обработки изображения(в зависимости от заданных параметров)
        switch (resize_method) {
        case 1:
            reduce_image_by_averaging(&input, Koef_K, &output);
            break;
        case 2:
            reduce_image_by_nearest_neighbor(&input, Koef_K, &output);
            break;
        default:
            return ERROR;
        }
    }
    
    // Создаем пирамиду для frame2(текущего кадра)
    for (int level = 1; level < pyramid_lvl; level++) {
        InputImage input;
        OutputImage output;
        
        input.data = pyramid2[level - 1];
        input.width = pyramid_width[level - 1];
        input.height = pyramid_height[level - 1];
        
        output.width = input.width / Koef_K;
        output.height = input.height / Koef_K;
        
        output.data = (unsigned char*)malloc(output.width * output.height * sizeof(unsigned char));
        if (output.data == NULL){
            pyramidMemoryDelete(pyramid1, pyramid2, pyramid_lvl);
            return ERROR;
        }
        pyramid2[level] = output.data;
        
        //Разделение на методы обработки изображения(в зависимости от заданных параметров)
        switch (resize_method) {
        case 1:
            reduce_image_by_averaging(&input, Koef_K, &output);
            break;
        case 2:
            reduce_image_by_nearest_neighbor(&input, Koef_K, &output);
            break;
        default:
            return ERROR;
        }
    }
    
    // Обрабатываем каждую ключевую точку
    for (int i = 0; i < point_num; i++) {
        float dx = 0.0f, dy = 0.0f;
        int x = points_coords[i * 2];
        int y = points_coords[i * 2 + 1];
        
        // Проходим по уровням пирамиды сверху вниз
        for (int level = pyramid_lvl - 1; level >= 0; level--) {
            // Масштабируем исходные координаты для текущего уровня
            float scale_down = 1.0f / (1 << level);
            int level_x = (int)(x * scale_down);
            int level_y = (int)(y * scale_down);
    
            // Добавляем предполагаемое смещение с предыдущего уровня
            if (level < pyramid_lvl - 1) {
                level_x += (int)(dx * scale_down);
                level_y += (int)(dy * scale_down);
            }
    
            // Проверка границ
            if (level_x < 1 || level_x >= pyramid_width[level] - 1 ||
                level_y < 1 || level_y >= pyramid_height[level] - 1) {
                continue;
            }
    
            // Вычисляем смещение на текущем уровне
            FlowResult LKshift = Lukas_Kanade_point(
                pyramid1[level], pyramid2[level],
                pyramid_width[level], pyramid_height[level],
                level_x, level_y
            );
    
            // Накопление с масштабированием
            if (level == 0) {
                dx += LKshift.dx;
                dy += LKshift.dy;
            }
            else {
                dx = (dx + LKshift.dx) * Koef_K;
                dy = (dy + LKshift.dy) * Koef_K;
            }
        }
        
        // Записываем результат
        Output_data->keypoints_shift[i * 2] = dx;
        Output_data->keypoints_shift[i * 2 + 1] = dy;
    }
    
    // Очищаем память
    pyramidMemoryDelete(pyramid1, pyramid2, pyramid_lvl);
    return OK;
}

// Решение системы 2x2 по правилу Крамера(для решения уравнения в алгоритме Фернебека)
void solve_2x2(float A[2][2], float b[2], float* x, float* y) {
    float det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (det*det < 1e-8) {
        *x = 0;
        *y = 0;
        return;
    }
    float inv_det = 1.0f / det;
    *x = (b[0]*A[1][1] - b[1]*A[0][1]) * inv_det;
    *y = (A[0][0]*b[1] - A[1][0]*b[0]) * inv_det;
}

// вычисляем градиенты для всего изображения(для определения плотного оптического потока)
void compute_gradients(const unsigned char* img, int width, int height, float* grad_x, float* grad_y) {
    for (int y=1; y<height-1; y++) {
        for (int x=1; x<width-1; x++) {
            int idx = y*width + x;
            float gx = (float)(img[idx+1] - img[idx-1]) * 0.5f;
            float gy = (float)(img[(y+1)*width + x] - img[(y-1)*width + x]) * 0.5f;
            grad_x[idx] = gx;
            grad_y[idx] = gy;
        }
    }
    // границы
    for (int x=0; x<width; x++) {
        grad_x[x] = 0;
        grad_x[(height-1)*width + x] = 0;
        grad_y[x] = 0;
        grad_y[(height-1)*width + x] = 0;
    }
    for (int y=0; y<height; y++) {
        grad_x[y*width] = 0;
        grad_x[y*width + (width-1)] = 0;
        grad_y[y*width] = 0;
        grad_y[y*width + (width-1)] = 0;
    }
}

//Основная функция нахождения плотного оптичесокго потока алгоритмом Farneback
STATUS Farneback(const unsigned char* img1, const unsigned char* img2, const int width, const int height,
                float* flow_x, float* flow_y) {
    
    int size = width * height;

    // Заполнение смещения нулями
    memset(flow_x, 0, size * sizeof(float));
    memset(flow_y, 0, size * sizeof(float));

    // Выделяем память под градиентов всего изображения
    float* Ix = (float*)malloc(size * sizeof(float));
    if (Ix == NULL){
        return ERROR;
    }
    float* Iy = (float*)malloc(size * sizeof(float));
    if (Iy == NULL){
        return ERROR;
    }

    // Вычислить градиенты для первого(предыдущего) кадра
    compute_gradients(img1, width, height, Ix, Iy);

    // Вычисляем оптический потко для каждой точки кадра
    for (int y=1; y<height-1; y++) {
        for (int x=1; x<width-1; x++) {
            int idx = y*width + x;

            float gx = Ix[idx];
            float gy = Iy[idx];
            float It = (float)(img2[idx] - img1[idx]);

            // Предвычисляем квадраты для повторного использования
            const float gx2 = gx * gx;
            const float gy2 = gy * gy;
            const float gxgy = gx * gy;

            // Заполнение матрицы А и b
            float A[2][2] = {
                {gx2, gxgy},
                {gxgy, gy2}
            };
            
            // Вектор b
            float b[2] = {
                -gx * It,
                -gy * It
            };

            //Инициализируем функцию решения матричного уравнения
            solve_2x2(A, b, &flow_x[idx], &flow_y[idx]);
        }
    }
    // Освобождаем память
    free(Ix);
    free(Iy);

    return OK;
}
