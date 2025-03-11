#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

void matrix_vandermonde(double *mat, size_t righe, size_t colonne, double *valorix)
//funzione che crea la matrice di Vandermonde, ovvero la matrice con 1, x, x^2, poi la seconda riga 1, x, x^2 etc

{
    for (size_t i = 0; i < righe; i++)
    {
        mat[i * colonne] = 1; //aggiungo l'1 come primo elemento
        for (size_t j = 0; j < colonne; j++)
        {
            mat[i * colonne + j] = pow(valorix[i], j); //prendo i valori di x e li elevo al numero di colonna ovvero j
            //è però computing intensive!!
        }
    }
}

int matrix_multiplication(double *m1, size_t rows1, size_t cols1, double *m2, size_t cols2, double *result) {
    if (m1 == NULL || m2 == NULL || result == NULL) {
        return -1;
    }

    for (size_t i = 0; i < rows1; i++) {
        for (size_t j = 0; j < cols2; j++) {
            result[i * cols2 + j] = 0;
            for (size_t k = 0; k < cols1; k++) {
                result[i * cols2 + j] += m1[i * cols1 + k] * m2[k * cols2 + j];
            }
        }
    }
    return 0;
}

void matrix_transpose(double *matrice1, size_t righe1, size_t colonne_1, double *matrice2)
{
    for (size_t i = 0; i < colonne_1; i++)
    {
        for (size_t j = 0; j < righe1; j++)
        {
            matrice2[i * righe1 + j] = matrice1[j * colonne_1 + i];
        }
    }
}

void matrix_print(double *matrice, int righe, int colonne)
{
    for (int i = 0; i < righe; i++)
    {
        for (int j = 0; j < colonne; j++)
        {
            printf("%f ", matrice[i * colonne + j]);
        }
        printf("\n");
    }
}

void matrix_swap_rows(double *mat, size_t riga_1, size_t riga_2, size_t colonne)
{
    double holder;
    for (size_t j = 0; j < colonne; j++)
    {
        holder = mat[riga_1 * colonne + j];
        mat[riga_1 * colonne + j] = mat[riga_2 * colonne + j];
        mat[riga_2 * colonne + j] = holder;
    }
}

void array_product(double *arr, size_t lenght, double value)
{
    for (size_t i = 0; i < lenght; i++)
    {
        arr[i] *= value;
    }
}

void array_clone(double *arr1, double *arr2, size_t lenght)
{
    for (size_t i = 0; i < lenght; i++)
    {
        arr1[i] = arr2[i];
    }
}

void matrix_enhance(double *mat, size_t righe, size_t colonne, double *mat_id)
//clono la matrice
{
    for (size_t i = 0; i < righe; i++)
    {
        for (size_t j = 0; j < colonne; j++)
        {
            mat_id[i * colonne * 2 + j] = mat[i * colonne + j];//i*colonne+j perche la matrice in realtà è un arraymat_id[i * colonne * 2 + colonne + j] = (double)(i == j);
        }
    }
//creo la matrice identità (FATTO SOPRA)
/*
if (i == j)
{
else
{
    mat_id[i * colonne*2 + colonne + j] = 0;
}
*/
//posso fare una cosa più furba

}

void array_subtract(double *arr1, double *arr2, size_t lenght)
{
    for (size_t i = 0; i < lenght; i++)
    {
        arr1[i] -= arr2[i];
    }
}

void matrix_gauss_inversion(double *mat, size_t righe, size_t colonne, double *mat_inv)
{
   double *mat_id = (double *) malloc(righe * colonne * 2 * sizeof(double));
   double *row_copy = (double *) malloc(colonne * 2 * sizeof(double));//allocato ora ma serve dopo
    
    if (!mat_id || !row_copy) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mat_id);
        free(row_copy);
        return;
    }

    matrix_enhance(mat, righe, colonne, mat_id);//riempio con la matrice identità
    //cerco il massimo nella colonna e lo metto in cima
    for (size_t i= 0; i<righe; i++)
    {
        double best_val= mat_id[i * colonne * 2 + i];//prendo il valore della diagonale, che è il valore massimo
        size_t best_row = i;
        for (size_t j = i+1; j<righe; j++)
        {
            if (fabs(mat_id[j * colonne * 2 + i]) > fabs(best_val))
            {
                best_val = mat_id[j * colonne * 2 + i];
                best_row = j;
            }
        }
        if (best_row != i)
        {
            matrix_swap_rows(mat_id, i, best_row, colonne * 2);
        }
        //normalizzo la riga
        //multiply_array(mat_id + i * colonne * 2, colonne * 2, 1.0 / mat_id[i * colonne * 2 + i]);
        array_product(&mat_id[i * colonne * 2], colonne * 2, 1.0 / mat_id[i * colonne * 2 + i]);

        //ora devo fare la sottrazione delle righe
        for (size_t k = 0; k < righe; k++)
        {
            if (k==i)
         {
              continue;
            }
    
        array_clone(row_copy, mat_id + i * colonne * 2, colonne * 2);
        array_product(row_copy, colonne * 2, mat_id[k * colonne * 2 + i]);
        array_subtract(mat_id + k * colonne * 2, row_copy, colonne * 2);
        }

    }
    //ora devo copiare la matrice inversa
    for (size_t i = 0; i < righe; i++)
    {
        for (size_t j = 0; j < colonne; j++)
        {
            mat_inv[i * colonne + j] = mat_id[i * colonne * 2 + colonne + j];
        }
    }
    free(mat_id);
    free(row_copy);
}