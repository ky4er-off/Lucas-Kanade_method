#ifndef IMG_SIZE_REDUCTION_H
#define IMG_SIZE_REDUCTION_H

#ifdef __cplusplus
extern "C" {
#endif

// Структура для входного изображения
typedef struct {
    const unsigned char* data;   // Массив пикселей - входное изображение
    int width;                   // Ширина изображения
    int height;                  // Высота изображения
} InputImage;

// Структура для выходного изображения
typedef struct {
    unsigned char* data;  // Массив пикселей - выходное изображение
    int width;            // Ширина изображения
    int height;           // Высота изображения
} OutputImage;

// Коды ошибок
typedef enum {
    RESIZE_SUCCESS = 0,
    RESIZE_ERROR_NULL_POINTER = -1,
    RESIZE_ERROR_INVALID_PARAMETER = -2,
    RESIZE_ERROR_SIZE_MISMATCH = -3,
    RESIZE_ERROR_MEMORY_ALLOCATION = -4
} ResizeErrorCode;

/// @brief Функция уменьшения изображения в k(целое число) раз методом усреднения по блокам
/// @param input    структура входных данных[in] (память выделять самотсоятельно)
/// @param k        коэффициент уменьшения(целое, положительное)
/// @param output   структура выходных данных[out] (память выделять самостоятельно)
/// @return         код ошибки (0 при успехе)
ResizeErrorCode reduce_image_by_averaging(const InputImage* input, int k, OutputImage* output);

/// @brief Функция уменьшения изображения в k(целое число) раз методом "ближайшего соседа"(выкидываем лишние пиксели)
/// @param input    входное изображение(структура)[in] 
/// @param k        коэффициент уменьшения (целое, положительное)
/// @param output   выходное изображение(структура)[out]
/// @return 
ResizeErrorCode reduce_image_by_nearest_neighbor(const InputImage* input, int k, OutputImage* output);

#ifdef __cplusplus
}
#endif

#endif