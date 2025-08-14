#ifndef MASK_VECTOR_ARITH_H
#define MASK_VECTOR_ARITH_H

#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "gadgets.h"
#include "masked_gfarith.h"
#include "params.h"
#include "masking_interface.h"
#include "src/blas.h"






// masked vector arithmetic
void masked_gf256v_madd(Masked *accu_c, const Masked *a, const Masked gf256_b, unsigned _num_byte);

void masked_gf256_madd(Masked *accu_c, const Masked *a, const Masked gf256_b, unsigned _num_byte);

// accu_c = accu_c + a
void masked_gf256v_add(Masked *accu_c, const Masked *a, unsigned _num_byte);
#endif
