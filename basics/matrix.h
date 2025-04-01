#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Matrix operations
void matrix_print(double *matrix, int rows, int cols);
void matrix_vandermonde(double *matrix, size_t rows, size_t cols, double *x_values);
void matrix_swap_rows(double *matrix, size_t row1, size_t row2, size_t cols);
void matrix_enhance(double *matrix, size_t rows, size_t cols, double *identity);
void matrix_gauss_inversion(double *matrix, size_t rows, size_t cols, double *inverse);
void matrix_transpose(double *matrix, size_t rows, size_t cols, double *result);
int matrix_multiplication(double *m1, size_t rows1, size_t cols1, double *m2, size_t cols2, double *result);

// Array operations
void array_product(double *array, size_t length, double value);
void array_clone(double *src, double *dst, size_t length);
void array_subtract(double *arr1, double *arr2, size_t length);

#endif /* MATRIX_H */

