#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <utils/mvt/core_sched.h> //Автомотическое определение распределения данных по ядрам
#include <cache.h>

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

// Структура для передачи параметров (минимальная)
typedef struct {
    int start_y;
    int end_y;
    int width;
    int height;
    unsigned char* img1;
    unsigned char* img2;
    bool is_add_core;
    float* flow_x;
    float* flow_y;
} ThreadData;

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

int process_slice(int param) {
    ThreadData* data = (ThreadData*)param;
    // Сброc кеша для доп. ядер
    if (data->is_add_core){
        inval_cache((unsigned long)data, (unsigned long)sizeof(*data));
        inval_cache((unsigned long)(data->img1+data->start_y*data->width), (unsigned long)((data->end_y-data->start_y)*data->width));
        inval_cache((unsigned long)(data->img2+data->start_y*data->width), (unsigned long)((data->end_y-data->start_y)*data->width));
    }
    
    // Предварительное вычисление указателей
    int width = data->width;
    unsigned char* img1 = data->img1;
    unsigned char* img2 = data->img2;
    float* flow_x = data->flow_x;
    float* flow_y = data->flow_y;
    
    for (int y = data->start_y; y < data->end_y; y++) {
        // Индексы строк для быстрого доступа
        int row_offset = y * width;
        int row_above = (y-1) * width;
        int row_below = (y+1) * width;
        
        for (int x = 1; x < width - 1; x++) {
            int idx = row_offset + x;

            // Градиенты
            float Ix = (float)(img1[idx+1] - img1[idx-1]) * 0.5f;
            float Iy = (float)(img1[row_below + x] - img1[row_above + x]) * 0.5f;
            float It = (float)(img2[idx] - img1[idx]);

            // Вычисление значений матрицы
            float gx2 = Ix * Ix;
            float gy2 = Iy * Iy;
            float gxgy = Ix * Iy;

            float A00 = gx2;
            float A01 = gxgy;
            float A11 = gy2;
            
            float b0 = -Ix * It;
            float b1 = -Iy * It;

            // Решение методом Крамера (оптимизированное)
            float det = A00 * A11 - A01 * A01;
            if (det*det > 1e-8) {
                float inv_det = 1.0f / det;
                flow_x[idx] = (b0 * A11 - b1 * A01) * inv_det;
                flow_y[idx] = (A00 * b1 - A01 * b0) * inv_det;
            } else {
                flow_x[idx] = 0;
                flow_y[idx] = 0;
            }
        }
    }
    // Сбрасываем в оперативку инфу c доп. ядер(вычисленные значения)
    if (data->is_add_core){
        flush_cache((unsigned long)(flow_x+data->start_y*data->width), (unsigned long)(((data->end_y-data->start_y)*data->width)*sizeof(float)));
        flush_cache((unsigned long)(flow_y+data->start_y*data->width), (unsigned long)(((data->end_y-data->start_y)*data->width)*sizeof(float)));
    }
    
    return OK;
}

STATUS Lukas_Kanade_dense(const LKDenseInput* Input_data, float* flow_x, float* flow_y) {
    const int width = Input_data->frame_width;
    const int height = Input_data->frame_height;
    unsigned char* img1 = Input_data->frame_prev;
    unsigned char* img2 = Input_data->frame_curr;
    const int cores_to_use = Input_data->num_cores;

    int size = width * height;
    memset(flow_x, 0, size * sizeof(float));
    memset(flow_y, 0, size * sizeof(float));

    // Последовательное выполнение(1 ядро)
    if (cores_to_use <= 1) {
        ThreadData main_data = {
            .start_y = 1,
            .end_y = height - 1,
            .width = width,
            .height = height,
            .img1 = img1,
            .img2 = img2,
            .flow_x = flow_x,
            .flow_y = flow_y,
            .is_add_core = FALSE,
        };
        return process_slice((int)&main_data) == OK ? OK : ERROR;
    }

    // Параллельное выполнение для 2+ ядер
    ThreadData core_data[4];
    int results[4] = {0};
    int data_per_core[4];
    
    // Оптимальное разделение
    int dist_error = distribute_data_cores(data_per_core, height-1, cores_to_use);
    if (dist_error != 0){
        return ERROR;
    }
    int current_y = 1;
    int additional_data = 0;
    
    // Запускаем задачи(заполняем структуры)
    for (int i = 0; i < cores_to_use; i++) {
        core_data[i] = (ThreadData){
            .start_y = current_y,
            .end_y = current_y + data_per_core[i],
            .width = width,
            .height = height,
            .img1 = img1,
            .img2 = img2,
            .flow_x = flow_x,
            .flow_y = flow_y
        };
        current_y = core_data[i].end_y;

        if (i == 0){
            core_data[i].is_add_core = FALSE;
        }
        else {
            core_data[i].is_add_core = TRUE;
            additional_data += data_per_core[i];
        }
    }

    //Передаем информацию для доп ядер в оперативку
    flush_cache((unsigned long)&core_data[1], (unsigned long)((cores_to_use-1)*sizeof(ThreadData)));
    flush_cache((unsigned long)(img1+data_per_core[0]*width), (unsigned long)(additional_data*width));
    flush_cache((unsigned long)(img2+data_per_core[0]*width), (unsigned long)(additional_data*width));
    
    //Загружаем доп. ядра
    for (int i = 1; i < cores_to_use; i++){
        eCoreNumber core_num;
        if (i == 1){
            core_num = core_1;
        }
        else if (i == 2){
            core_num = core_2;
        }
        else{
            core_num = core_3;
        }
        if (coreExecute(core_num, process_slice, (int)&core_data[i]) != OK) {
            return ERROR;
        }
    }

    //Запускаем основное ядро
    process_slice((int)&core_data[0]);
    
    // Ждем результат с дополнительных ядер
    for (int i = 1; i < cores_to_use; i++) {
        eCoreNumber core_num = (i == 1) ? core_1 : 
                               (i == 2) ? core_2 : core_3;
        
        if (coreWait(core_num, 250, &results[i]) != OK || results[i] != OK) {
            return ERROR;
        }
    }
    //Очистка кеша основного ядра
    inval_cache((unsigned long)(flow_x+data_per_core[0]*width), (unsigned long)(additional_data*width*sizeof(float)));
    inval_cache((unsigned long)(flow_y+data_per_core[0]*width), (unsigned long)(additional_data*width*sizeof(float)));

    return OK;
}
