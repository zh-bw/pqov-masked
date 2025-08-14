#ifndef MASK_GFARITH_H
#define MASK_GFARITH_H

#include <stddef.h>
#include <stdint.h>
#include "random.h"
#include "gadgets.h"
#include "src/gf16.h"
#include "masking_interface.h"
#include "params.h"



#ifdef COUNT
extern uint64_t count_op;
#endif
//masked GF(256) arithmetics
void masked_gf256_mul(Masked *z, const Masked a, const Masked b);
void masked_gf256_sqr(Masked *z, const Masked a);
void masked_gf256_inv(Masked *z, const Masked a);


void masked_gf256v_mul(Masked_u32 *z, const Masked_u32 a, const Masked b);
void masked_gf256v_sqr(Masked_u32 *z, const Masked_u32 a);
void masked_gf256v_mul_u32_u32(Masked_u32 *z, const Masked_u32 a, const Masked_u32 b);
#endif
