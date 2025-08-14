#ifndef MASK_INTER_H
#define MASK_INTER_H

#include <stdint.h>
#include "src/params.h"

#define MAX_O_BYTE 96

typedef struct
{
    uint8_t shares[N_SHARES][MAX_O_BYTE * MAX_O_BYTE]; // 96 is MAX_O_BYTE
} Masked_matrix;

typedef struct Masked{
  uint8_t shares[N_SHARES];
} Masked;

typedef struct Masked_u32{
  uint32_t shares[N_SHARES];
} Masked_u32;

#endif
