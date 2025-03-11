#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

//[VANDER^T * VANDER]^-1 * VANDER^T * ARRAY_Y = [poly_coefficients]

/*
|1 x1 x1^2 x1^3 ... x1^n|   |a0|   |y1|
|1 x2 x2^2 x2^3 ... x2^n|   |a1|   |y2|
|1 x3 x3^2 x3^3 ... x3^n| * |a2| = |y3|
|1 x4 x4^2 x4^3 ... x4^n|   |a3|   |y4|
|1 x5 x5^2 x5^3 ... x5^n|   |a4|   |y5|
...
|1 xn xn^2 xn^3 ... xn^n|   |an|   |yn|
*/

void poly_fitting(double *arrayx, double *arrayy, size_t lenghtx, double *poly_coefficients, size_t degree)
{
    //alloco la matrice varmondiana con dimensioni (campioni in x * grado del polinomio)
    double *vandermonde_matrix = (double *) malloc(lenghtx * degree * sizeof(double));
    //alloco la matrice varmondiana trasposta, dimensioni opposte a vandermonde_matrix
    double *matrice_trasposta = (double *) malloc(degree * lenghtx * sizeof(double));
    //alloco la matrice risultante
    double *matrice_prodotto = (double *) malloc(degree * degree * sizeof(double));

    // Controllo se la memoria è stata allocata correttamente
    if (!vandermonde_matrix || !matrice_trasposta || !matrice_prodotto) {
        // Handle memory allocation failure
        fprintf(stderr, "Memory allocation failed\n");
        // Free any successfully allocated memory
        // Return or exit
    }

    //creo la matrice varmondiana con dimensioni (campioni in x * grado del polinomio)
    matrix_vandermonde(vandermonde_matrix, lenghtx, degree, arrayx);
    //creo la matrice vermondiana trasposta
    matrix_transpose(vandermonde_matrix, lenghtx, degree, matrice_trasposta);
    //alloco la matrice risultante
    //VANDER^T * VANDER
    matrix_multiplication(matrice_trasposta, degree, lenghtx, vandermonde_matrix, degree, matrice_prodotto);
    //[VANDER^T * VANDER]^-1
    matrix_gauss_inversion(matrice_prodotto, degree, degree, matrice_prodotto);
    // ... * ARRAY_Y 
    matrix_multiplication(matrice_prodotto, degree, degree, matrice_trasposta, lenghtx, matrice_trasposta);
    //... * ARRAY_Y = [poly_coefficients]
    matrix_multiplication(matrice_trasposta, degree, lenghtx, arrayy, 1, poly_coefficients);
    //printo il risultato
    matrix_print(poly_coefficients, degree, 1);
    //libero la memoria
    free(vandermonde_matrix);
    free(matrice_trasposta);
    free(matrice_prodotto);
}