#include "masked_vector_arith.h"
#include "debug.h"
#include "src/ov_blas.h"

// we assume _num_byte is multiple of 4
// Todo: fix when _num_byte is not multiple of 4
void masked_gf256v_madd(Masked *accu_c, const Masked *a, const Masked gf256_b, unsigned _num_byte){
    int i = 0,k;
    Masked_u32 u32;
    Masked u8[4];

    while (_num_byte >=4)
    {
        // batch processing, we process 4 blocks at a time
        for (k = 0; k < N_SHARES; k++) {
            u32.shares[k] = ((uint32_t)a[i + 3].shares[k] << 24) | 
                            ((uint32_t)a[i + 2].shares[k] << 16) |
                            ((uint32_t)a[i + 1].shares[k] << 8)  |
                            ((uint32_t)a[i + 0].shares[k]);
        }
        
        masked_gf256v_mul(&u32, u32, gf256_b);



        // load back
        for (k = 0; k < N_SHARES; k++)
        {
            
            accu_c[i + 0].shares[k] ^= (uint8_t)(u32.shares[k] & 0xFF);
            accu_c[i + 1].shares[k] ^= (uint8_t)((u32.shares[k] >> 8) & 0xFF);
            accu_c[i + 2].shares[k] ^= (uint8_t)((u32.shares[k] >> 16) & 0xFF);
            accu_c[i + 3].shares[k] ^= (uint8_t)((u32.shares[k] >> 24) & 0xFF);
            #ifdef COUNT
            count_op += 4;
            #endif
        }

        i += 4;
        _num_byte -= 4;
    }

    if (0 == _num_byte){
        return;
    }
}


void masked_gf256v_add(Masked *accu_c, const Masked *a, unsigned _num_byte){
    int k;
    for (k = 0; k < N_SHARES; k++){
        uint8_t temp_accu[_num_byte], temp_a[_num_byte];
        
        // Extract the k-th share from both vectors
        for (int i = 0; i < _num_byte; i++){
            temp_accu[i] = accu_c[i].shares[k];
            temp_a[i] = a[i].shares[k];
        }
        
        // Apply gf256v_add to the k-th share
        gf256v_add(temp_accu, temp_a, _num_byte);
        
        // Store the result back to the k-th share
        for (int i = 0; i < _num_byte; i++){
            accu_c[i].shares[k] = temp_accu[i];
        }
    }
}

void masked_gf256_madd(Masked *accu_c, const Masked *a, const Masked gf256_b, unsigned _num_byte){
    int i;
    Masked temp;

    for (i = 0; i < _num_byte; i++){
        masked_gf256_mul(&temp, a[i], gf256_b);
        for (int k = 0; k < N_SHARES; k++){
            accu_c[i].shares[k] ^= temp.shares[k];
            #ifdef COUNT
            count_op ++;
            #endif
        }
    }
}