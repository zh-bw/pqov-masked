#include <math.h>
#include <stdio.h>
#include <string.h>
#include "masking/random.h"
#include "masked_gfarith.h"
#include "debug.h"
#include "masked_vector_arith.h"
#include "matrix_arith.h"
#include "gadgets.h"
#include "cpucycles.h"
#include "masking_interface.h"
#include "src/ov_blas.h"
#include "src/ov.h"
#include "masking/masked_ov.h"
#include "src/api.h"
#include "utils/fips202.h"
#include "utils/utils_randombytes.h"
#include "src/blas_matrix.h"

#define ITER 500
uint64_t start, stop;

void timing_gadgets(){
    Masked_matrix masked_matA, masked_matD;
    uint8_t unmasked_matA[_O_BYTE*_O_BYTE], matC[_O_BYTE*_O_BYTE], unmasked_matD[_O_BYTE*_O_BYTE];
    uint8_t unmasked_matA_inv[_O_BYTE*_O_BYTE], unmasked_matD_inv[_O_BYTE*_O_BYTE];
    uint8_t r[_O_BYTE*N_SHARES];
    Masked random_matrix[_O_BYTE * _O_BYTE];
    unsigned ok = 0;
    Masked masked_vec_b[_O_BYTE];
    for (int i = 0; i < _O_BYTE; i++){
        for (int k = 0; k < N_SHARES; k++){
            randombytes(&masked_vec_b[i].shares[k], 1);
            // masked_vec_b[i].shares[k] = rand8();
        }
    }

    Masked a[_O_BYTE], b[_O_BYTE];
    for (int i = 0; i < _O_BYTE; i++) {
        for (int j = 0; j < N_SHARES; j++) {
            randombytes(&a[i].shares[j], 1);
            randombytes(&b[i].shares[j], 1);
            // a[i].shares[j] = rand8();
            // b[i].shares[j] = rand8();
        }
    }

    // Try repeatedly until the combined (unmasked) matrix is invertible.
    do {
        for (int i = 0; i < N_SHARES; i++){
            randombytes(masked_matA.shares[i], _O_BYTE * _O_BYTE);
        }

        memcpy(unmasked_matA, masked_matA.shares[0], _O_BYTE * _O_BYTE);
        for (int i = 1; i < N_SHARES; i++) {
            gf256mat_add(unmasked_matA, unmasked_matA, masked_matA.shares[i], _O_BYTE, _O_BYTE);
            // gf256mat_add(unmasked_matD, unmasked_matD, masked_matD.shares[i], _O_BYTE, _O_BYTE);
        }
        ok = gf256mat_inv(unmasked_matA, unmasked_matA_inv, _O_BYTE); 
        // & gf256mat_inv(unmasked_matD, unmasked_matD_inv, _O_BYTE);
    } while (!ok);
    puts("\n**************** Timing masked linear equations solver *******************\n");
    // printf("\n* ------------------- Timing masked linear equations solver -------------------\n");
    Masked masked_vec_r[_O_BYTE];
    printf("Masking order: %d\n", MASKING_ORDER);
    

    start = cpucycles();
    for (int i = 0; i < ITER; i++) {
        masked_linear_equation_solver_rankcheck(masked_vec_r, masked_matA, masked_vec_b);
    }
    stop = cpucycles();
    printf("\n* Avg speed masked linear equations solver via multiplicative shares %.1f cycles.\n", (double)(stop-start)/(ITER));

    start = cpucycles();
    for (int i = 0; i < ITER; i++) {
        masked_linear_equation_solver_alternative(masked_vec_r, masked_matA, masked_vec_b);
    }
    stop = cpucycles();
    // printf("\n* ------------------- Timing masked linear equations solver alternative -------------------\n");
    printf("\n* Avg speed masked linear equations solver via additive shares: %.1f cycles.\n", (double)(stop-start)/(ITER ));
    puts("\n*********************************************************\n");
}

static inline
void random_boolean_mask(unsigned char *masked_m, unsigned char *m, int length)
{
  int i, k;
  unsigned char share;

  // generate random 8bits mask
  for (i = 0; i < length; i++)
  {
    share = m[i];
    for (k = 0; k < MASKING_ORDER; k++)
    {
      masked_m[i + length*k] = rand8();
      share ^= masked_m[i + length*k];
    }
    masked_m[i + length*MASKING_ORDER] = share;
  }

}
static inline
void convert_to_Masked(Masked *masked_r, uint8_t *r, int length)
{
    size_t bufsize = (size_t)length * (size_t)N_SHARES;
    uint8_t *temp = malloc(bufsize);
    if (!temp) {
        perror("malloc for temp in convert_to_Masked");
        exit(1);
    }

    random_boolean_mask(temp, r, length);

    for (int i = 0; i < length; i++) {
        for (int k = 0; k < N_SHARES; k++) {
            masked_r[i].shares[k] = temp[i + k * length];
        }
    }

    free(temp);
}


void timing_masked_ov(){
        // generate keypair
    cpk_t pk;
    sk_t sk;

    unsigned char sk_seed[LEN_SKSEED];
    randombytes(sk_seed, LEN_SKSEED);

    generate_keypair_pkc(&pk, &sk, sk_seed);

    // avoid stack overflow
    static masked_sk_t masked_sk;

    uint8_t masked_sk_seed[LEN_SKSEED * N_SHARES];
    random_boolean_mask(masked_sk_seed, sk_seed, LEN_SKSEED);

    printf("\n* ------------------- Timing unmasked ov-pkc keygen -------------------\n");
    start = cpucycles();
    for (int i = 0; i < 50; i++) {
        generate_keypair_pkc(&pk, &sk, sk_seed);    
    }
    stop = cpucycles();
    printf("\n* Avg speed unmasked keygen: %.1f cycles.\n", (double)(stop-start)/(50));

    puts("\n****************** Timing masked ov-pkc ******************\n");
    printf("Masking order: %d\n", MASKING_ORDER);
    start = cpucycles();
    for (int i = 0; i < ITER; i++) {
        masked_generate_keypair_pkc(&pk, &masked_sk, masked_sk_seed);
    }
    stop = cpucycles();
    printf("\n* Avg speed masked keygen: %.1f cycles.\n", (double)(stop-start)/(ITER));
    // generate message
    unsigned char m[53];
    unsigned long long mlen = 53;
    for (unsigned i = 0; i < 53; i++) m[i] = i;

    // sign
    uint8_t signature[_PUB_N_BYTE + _SALT_BYTE];
    
    printf("\n* ------------------- Timing unmasked ov-pkc sign -------------------\n");
    
    start = cpucycles();
    for (int i = 0; i < 50; i++) {
        ov_sign(signature, &sk, m, mlen);
    }
    stop = cpucycles();
    printf("\n* Avg speed unmasked signing: %.1f cycles.\n", (double)(stop-start)/(50));

    // printf("\n* ------------------- Timing masked ov-pkc -------------------\n");
    // printf("Masking order: %d\n", MASKING_ORDER);
    Masked masked_signature[_PUB_N_BYTE + _SALT_BYTE];
    
    start = cpucycles();
    for (int i = 0; i < ITER; i++) {
        masked_ov_sign(signature, &masked_sk, m, mlen);
    }
    stop = cpucycles();
    printf("\n* Avg speed masked signing: %.1f cycles.\n", (double)(stop-start)/(ITER));

    puts("\n*********************************************************\n");
}

void timing_randomness(){
    uint8_t r;
    r = rand8();

    start = cpucycles();
    for (int i = 0; i < 20000; i++){
        r = rand8();
    }
    stop = cpucycles();
    printf("\n* ------------------- Timing randomness -------------------\n");
    printf("\n* Avg speed randomness: %.1f cycles.\n", (double)(stop-start)/(20000));
}


int main() {
    // timing_gadgets();
    // timing_masked_ov();
    timing_randomness();
    return 0;
}