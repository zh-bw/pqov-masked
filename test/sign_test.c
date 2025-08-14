#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "src/ov.h"
#include "utils/utils_randombytes.h"
#include "src/params.h"
#include "masking/random.h"
#include "masking/debug.h"
#include "masking/masked_ov.h"
#include "masking/cpucycles.h"

#define ITER 500
uint64_t start, stop;

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


// pkc mode
void test_maskedov() {
    // generate keypair
    cpk_t pk, pk_test;
    sk_t sk;
    unsigned char sk_seed[LEN_SKSEED];
    randombytes(sk_seed, LEN_SKSEED);

    generate_keypair_pkc(&pk, &sk, sk_seed);

    uint8_t masked_sk_seed[LEN_SKSEED * N_SHARES];

    random_boolean_mask(masked_sk_seed, sk_seed, LEN_SKSEED);
    // avoid stack overflow
    static masked_sk_t masked_sk_test;

    printf("start masked keypair generation...\n");
    masked_generate_keypair_pkc(&pk_test, &masked_sk_test, masked_sk_seed);
    // check if pk_test matches pk
    if (memcmp(pk.pk_seed, pk_test.pk_seed, LEN_PKSEED) != 0) {
        printf("Error: pk_seed does not match\n");
        return;
    }
    if (memcmp(pk.P3, pk_test.P3, LEN_PKSEED) != 0) {
        printf("Error: pk.P3 does not match\n");
        return;
    }

    // generate message
    unsigned char m[53];
    unsigned long long mlen = 53;
    for (unsigned i = 0; i < 53; i++) {
        m[i] = i;
    }

    // sign
    uint8_t signature[_PUB_N_BYTE + _SALT_BYTE];
    ov_sign(signature, &sk, m, mlen);

    printf("======= Masked OV Signature Test ========\n");

    // Masked masked_signature[_PUB_N_BYTE + _SALT_BYTE];
    uint8_t rec_sig[_PUB_N_BYTE + _SALT_BYTE];
    masked_ov_sign(rec_sig, &masked_sk_test, m, mlen);

    for (int i = 0; i < _PUB_N_BYTE + _SALT_BYTE; i++) {
        if (rec_sig[i] != signature[i]) {
            printf("Error at index %d: expected %d, got %d\n", i, signature[i], rec_sig[i]);
            return;
        }
    }
    printf("Masked keypair generated successfully.\n");
    printf("Masked signature matches unmasked signature.\n");
}



void test_UOV() {
    // generate keypair
    cpk_t pk;

    uint8_t masked_sk_seed[LEN_SKSEED * N_SHARES];
    printf("Generating masked sk_seed...\n");
    randombytes(masked_sk_seed, LEN_SKSEED * N_SHARES);
    // avoid stack overflow
    static masked_sk_t masked_sk;


    printf("start masked keypair generation...\n");
    masked_generate_keypair_pkc(&pk, &masked_sk, masked_sk_seed);

    // generate message
    unsigned char m[53];
    unsigned long long mlen = 53;
    for (unsigned i = 0; i < 53; i++) {
        m[i] = i;
    }

    // sign
    uint8_t signature[_PUB_N_BYTE + _SALT_BYTE];

    printf("======= Masked OV Signature Test ========\n");

    masked_ov_sign(signature, &masked_sk, m, mlen);
    printf("Masked signature generated successfully.\n");
    // verify in pkc mode
    int rc = ov_expand_and_verify(m, mlen, signature, &pk);
    if (rc != 0) {
        printf("Signature verification failed.\n");
    } else {
        printf("Signature verification succeeded.\n");
    }
}

int main(){
    test_UOV();
    return 0;
}