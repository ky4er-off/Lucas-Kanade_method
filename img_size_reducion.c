#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "img_size_reducion.h"

// Оптимизированная версия для k=2 
static ResizeErrorCode reduce_by_2_fast(const InputImage* input, OutputImage* output) {
    const int input_width = input->width;
    const unsigned char* src = input->data;
    unsigned char* dst = output->data;
    
    const int output_width = output->width;
    const int output_height = output->height;
    
    // Используем uint32_t для загрузки 4 пикселей сразу
    for (int y = 0; y < output_height; y++) {
        const int src_y = y * 2;
        
        for (int x = 0; x < output_width; x++) {
            const int src_x = x * 2;
            
            // Загружаем 4 пикселя как 32-битное слово
            uint32_t block = *(uint32_t*)&src[src_y * input_width + src_x];
            
            // Извлекаем отдельные байты и суммируем
            unsigned int sum = ((block >> 0) & 0xFF) +
                               ((block >> 8) & 0xFF) +
                               ((block >> 16) & 0xFF) +
                               ((block >> 24) & 0xFF);
            
            dst[y * output_width + x] = (unsigned char)(sum >> 2); // деление на 4 сдвигом
        }
    }
    
    return RESIZE_SUCCESS;
}

//Оптимизированная версия для любого k
static ResizeErrorCode reduce_by_k_fast(const InputImage* input, int k, OutputImage* output) {
    const int input_width = input->width;
    const unsigned char* src = input->data;
    unsigned char* dst = output->data;
    
    const int output_width = output->width;
    const int output_height = output->height;
    const int block_size = k * k;
    
    // Предвычисляем множитель для быстрого деления
    const unsigned int inv_block_size = (1 << 16) / block_size; // для фиксированной точки
    
    for (int y = 0; y < output_height; y++) {
        const int src_y_start = y * k;
        
        for (int x = 0; x < output_width; x++) {
            const int src_x_start = x * k;
            unsigned int sum = 0;
            
            // Оптимизированный внутренний цикл
            for (int dy = 0; dy < k; dy++) {
                const unsigned char* src_row = &src[(src_y_start + dy) * input_width + src_x_start];
                
                // Частичная развертка цикла для улучшения производительности
                int dx = 0;
                // Обрабатываем по 4 пикселя за итерацию если k >= 4
                for (; dx <= k - 4; dx += 4) {
                    sum += src_row[dx] + src_row[dx + 1] + src_row[dx + 2] + src_row[dx + 3];
                }
                // Обрабатываем остаток
                for (; dx < k; dx++) {
                    sum += src_row[dx];
                }
            }
            
            // Быстрое деление через умножение и сдвиг для фиксированной точки
            // (sum * inv_block_size) >> 16 эквивалентно sum / block_size
            dst[y * output_width + x] = (unsigned char)((sum * inv_block_size) >> 16);
        }
    }
    
    return RESIZE_SUCCESS;
}

// Главная функция с выбором оптимального метода
ResizeErrorCode reduce_image_by_averaging(const InputImage* input, int k, OutputImage* output) {
    // Проверка входных параметров
    if (input == NULL || output == NULL) {
        return RESIZE_ERROR_NULL_POINTER;
    }
    
    if (input->data == NULL || k <= 0) {
        return RESIZE_ERROR_INVALID_PARAMETER;
    }
    
    // Проверка, что размеры кратны k
    if (input->width % k != 0 || input->height % k != 0) {
        return RESIZE_ERROR_SIZE_MISMATCH;
    }
    
    // Вычисляем размеры выходного изображения
    output->width = input->width / k;
    output->height = input->height / k;
    
    // Проверка, что после уменьшения есть хотя бы один пиксель
    if (output->width <= 0 || output->height <= 0) {
        return RESIZE_ERROR_INVALID_PARAMETER;
    }
    
    // Проверка, что выходной буфер выделен
    if (output->data == NULL) {
        return RESIZE_ERROR_NULL_POINTER;
    }
    
    // Выбираем самый быстрый метод в зависимости от k
    switch (k) {
        case 2:
            return reduce_by_2_fast(input, output);
        default:
            return reduce_by_k_fast(input, k, output);
    }
}

ResizeErrorCode reduce_image_by_nearest_neighbor(const InputImage* input, int k, OutputImage* output) {
    // Проверка входных параметров
    if (input == NULL || output == NULL) {
        return RESIZE_ERROR_NULL_POINTER;
    }
    
    if (input->data == NULL || k <= 0) {
        return RESIZE_ERROR_INVALID_PARAMETER;
    }
    
    // Проверка, что размеры кратны k
    if (input->width % k != 0 || input->height % k != 0) {
        return RESIZE_ERROR_SIZE_MISMATCH;
    }
    
    // Вычисляем размеры выходного изображения
    output->width = input->width / k;
    output->height = input->height / k;
    
    // Проверка, что после уменьшения есть хотя бы один пиксель
    if (output->width <= 0 || output->height <= 0) {
        return RESIZE_ERROR_INVALID_PARAMETER;
    }
    
    // Проверка, что выходной буфер выделен
    if (output->data == NULL) {
        return RESIZE_ERROR_NULL_POINTER;
    }
    
    const unsigned char* src = input->data;
    unsigned char* dst = output->data;
    const int src_width = input->width;
    const int dst_width = output->width;
    const int dst_height = output->height;
    
    //реализация алгоритма "ближайшего соседа" (Оптимизация: используем указатели для обхода массивов)
    unsigned char* dst_ptr = dst;
    
    for (int y = 0; y < dst_height; y++) {
        // Вычисляем индекс начала строки в исходном изображении
        int src_row_start = (y * k) * src_width;
        
        for (int x = 0; x < dst_width; x++) {
            // Копируем пиксель и сразу перемещаем указатель
            *dst_ptr++ = src[src_row_start + (x * k)];
        }
    }
    
    return RESIZE_SUCCESS;
}