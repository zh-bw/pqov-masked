#ifndef _MATRIX_ARITH_H
#define _MATRIX_ARITH_H

#include <stdint.h>
#include "src/blas_matrix.h"
#include "masked_vector_arith.h"
#include "params.h"
#include "masking_interface.h"




///////////////// Section: multiplications  ////////////////////////////////


// column-major matrix multiplication
/// @brief matrix-matrix multiplication:  C = A * B , in GF(256)
void gf256mat_mul_square(uint8_t *matC, const uint8_t *matA, const uint8_t *matB, unsigned width);

// column-major matrix addition
void gf256mat_add(uint8_t *r, const uint8_t *mat_A, const uint8_t *mat_B, 
                  unsigned n_A_vec_byte, unsigned n_A_width);



//////////////////////////  Gaussian for solving lienar equations ///////////////////////////


/// @brief Inversion of a square matrix, in GF(256)
unsigned gf256mat_inv(uint8_t *A, uint8_t *A_inv, unsigned n);

// generate a random invertible square matrix
void gf256mat_gen_upper(uint8_t *mat, unsigned n);
// generate a random square matrix
void gf256mat_gen(uint8_t *mat, unsigned n);

/* masked matrix-vector product */
// c = A * b
void masked_gf256mat_prod(Masked *c, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked *b);

void masked_gf256mat_prod_mtimes(Masked *c, const Masked* matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked *b);

void masked_gf256mat_prod_unbatched(Masked *c, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked *b);

void half_masked_gf256matvec_prod(Masked *c, const Masked* matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

unsigned check_fullrank(const uint8_t *mat, unsigned height, unsigned width);

// r = vT * P1 * v
void masked_quad_trimat_eval(Masked *r, const uint8_t *mat_P1, Masked *v,unsigned dim, unsigned size_batch);

/* ------------- used in keygen ----------------- */
void masked_calculate_F2_P3(Masked *masked_S, unsigned char *P3, const unsigned char *P1, const Masked *P2, Masked *masked_O);

// bC += btriA * B, bC and B are masked, btriA is unmasked 
void masked_batch_trimat_madd(Masked *bC, const unsigned char *btriA, 
                            const Masked *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch );

// bC += btriA^Tr * B, bC and B are masked, btriA is unmasked
void masked_batch_trimatTr_madd(Masked *bC, const unsigned char *btriA, 
                            const Masked *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch );

// bC = A^Tr * bB, all masked
void masked_batch_upper_matTr_x_mat( Masked *bC,
                                    const Masked *A_to_tr, unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
                                    const Masked *bB, unsigned Bwidth, unsigned size_batch );


#endif
