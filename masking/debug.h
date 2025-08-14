#ifndef DEBUG_H
#define DEBUG_H
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "gadgets.h"
// #include "gfarith.h" 
#include "masked_vector_arith.h"
#include "matrix_arith.h"
#include "masking_interface.h"

void print_masked_bool(Masked* y);
void print_bitstring(uint8_t * bs, int len);
void unmask_bitstring(uint8_t * bs, int len);

void combine_arith_shares(uint8_t *unmasked_vec, uint8_t *masked_vec, unsigned n_A_vec_byte);
void vec_combine_arith_shares(uint8_t *unmasked_vec, const Masked *masked_vec, unsigned n_A_vec_byte);
void print_masked_arith(Masked* x);
void combine_boolean_shares(uint8_t *unmasked_buf, uint8_t *bs, int len);
void print_matrix(const uint8_t *matrix, unsigned n_rows, unsigned n_cols);
void print_matrix_rowmajor(const uint8_t *matrix, unsigned n_rows, unsigned n_cols);
void mat_combine_arith_shares(uint8_t *unmasked_mat, const Masked_matrix masked_mat, unsigned n_A_vec_byte, unsigned n_A_width);
#endif