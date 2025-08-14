#ifndef _GADGETS_H
#define _GADGETS_H

#include <stdint.h>
#include <stdio.h>
#include "masked_vector_arith.h"
#include "matrix_arith.h"
#include "random.h"
#include "params.h"
#include "masked_gfarith.h"
#include "masking_interface.h"


/**
 * @brief Performs a t-NI  refresh operation on a masked matrix.
 *
 * This function applies a linear arithmetic refresh to the given masked matrix `matA`.
 * The refresh operation ensures that the masked representation remains secure by 
 * re-randomizing the shares while preserving the underlying secret.
 *
 * @param matA Pointer to the masked matrix to be refreshed.
 * @param n_A_vec_byte The number of bytes in each vector of the matrix.
 * @param n_A_width The width (number of columns) of the matrix.
 */
void mat_linear_arithmetic_refresh(Masked_matrix *matA, unsigned n_A_vec_byte, unsigned n_A_width);

/**
 * @brief Performs a t-SNI refresh operation on a masked matrix.
 *
 * This function applies a arithmetic refresh to the given masked matrix `matA`.
 * The refresh operation ensures that the masked representation remains secure by 
 * re-randomizing the shares while preserving the underlying secret.
 *
 * @param matC Pointer to the masked matrix whose masks will be refreshed.
 * @param matA Pointer to the masked matrix used as a reference for refreshing.
 * @param n_A_vec_byte The number of bytes in each vector of matrix A.
 * @param n_A_width The width (number of columns) of matrix A.
 */
void mat_refresh_masks(Masked_matrix *matC, Masked_matrix *matA, unsigned n_A_vec_byte, unsigned n_A_width);

void gf256_linear_arithmetic_refresh(Masked *a);
void gf256_linear_arithmetic_refresh_u32(Masked_u32 *a);

/// @brief Refresh a masked value using arithmetic refreshing
/// @param[out] r - the refreshed output
/// @param[in] a - the input masked value to refresh
void gf256_arithmetic_refresh(Masked *z, Masked *a);

/// @brief Refresh a masked vector using arithmetic refreshing
/// @param[out] r - the refreshed output vector
/// @param[in] a - the input masked vector to refresh
/// @param[in] n_vec_byte - the number of elements in the vector
void vec_arithmetic_refresh(Masked *r, Masked *a, unsigned n_vec_byte);


/**
 * @brief Computes the product of a GF(256) matrix and a masked vector.
 *
 * This function performs the multiplication of a GF(256) matrix (`matA`)
 * with a masked vector (`masked_b`) and stores the result in the masked vector
 * `masked_c`. The matrix is represented in col-major order.
 *
 * @param[out] masked_c       Pointer to the output masked vector to store the result.
 * @param[in]  matA           Pointer to the GF(256) matrix in col-major order.
 * @param[in]  n_A_vec_byte   Number of bytes in each row of the matrix.
 * @param[in]  n_A_width      Number of columns in the matrix.
 * @param[in]  masked_b       Pointer to the input masked vector.
 */
void half_masked_gf256mat_prod(Masked *masked_c,
                               const uint8_t *matA,
                               unsigned n_A_vec_byte,
                               unsigned n_A_width,
                               const Masked *masked_b);

/**
 * @brief Performs a masked multiplication of two square matrices over GF(256).
 *
 * This function computes the product of two masked matrices `matA` and `matB`,
 * and stores the result in the masked matrix `matC`. All matrices are assumed
 * to be square and represented in the masked format.
 *
 * @param[out] matC Pointer to the masked matrix where the result will be stored.
 * @param[in]  matA The first masked input matrix.
 * @param[in]  matB The second masked input matrix.
 */
void masked_gf256mat_mul_square(Masked_matrix *matC, const Masked_matrix matA, const Masked_matrix matB);

void masked_gf256mat_mul_square_lowlevel(Masked_matrix *matC, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked_matrix matB);

void masked_gf256mat_mul_square_ISW(Masked_matrix *matC, const Masked_matrix matA, const Masked_matrix matB);

/**
 * @brief Solves the linear equation Ax = b, where A is a masked matrix, x is a masked vector, 
 *        and b is a masked vector.
 *
 * @param[out] masked_vec_r  Pointer to the masked vector that will store the solution x.
 * @param[in]  masked_matA   The masked matrix A.
 * @param[in]  masked_vec_b  Pointer to the masked vector b.
 *
 * @return An unsigned integer indicating the success or failure of the operation.
 *         0 indicates success, while non-zero values indicate on solution.
 */
// solving Ax = b, where A is a masked matrix, x is a masked vector, and b is a masked vector

unsigned masked_linear_equation_solver(Masked *masked_vec_r, const Masked_matrix masked_matA, 
                                   const Masked *masked_vec_b);
                                   
void masked_linear_equation_solver_rankcheck(Masked *masked_vec_r, const Masked_matrix masked_matA, 
                                   const Masked *masked_vec_b);



void masked_linear_equation_solver_alternative(Masked *masked_vec_r, const Masked_matrix masked_matA, 
                                    const Masked *masked_vec_b);
#endif
