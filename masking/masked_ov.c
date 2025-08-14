#include "masked_ov.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "utils/fips202.h"
#include "utils/utils_randombytes.h"
#include "utils/utils_prng.h"
#include "debug.h"

#define MAX_ATTEMPT_VINEGAR  256


int masked_ov_sign( uint8_t *signature, const masked_sk_t *masked_sk, const uint8_t *message, size_t mlen ){
    int i, j, k;
    uint8_t salt[_SALT_BYTE];
    // fix salt for testing
    // for (i = 0; i < _SALT_BYTE; i++)
    //     salt[i] = i+1;
    randombytes( salt, _SALT_BYTE );
    
    uint8_t m_salt[_SALT_BYTE + mlen];
    uint8_t masked_m_salt_secret[(mlen + _SALT_BYTE + LEN_SKSEED) * N_SHARES];
    uint8_t masked_h_vinegar[(mlen + _SALT_BYTE + LEN_SKSEED + 1) * N_SHARES]; // m||salt||sk_seed||ctr
    uint8_t vinegar[_V_BYTE * N_SHARES], x_o1[_O_BYTE];
    Masked mat_l1[_O * _O_BYTE], masked_r_l1_F1[_O_BYTE], x[_O_BYTE];
    Masked_matrix masked_mat_l1;

    memcpy(m_salt, message, mlen);
    memcpy(m_salt + mlen, salt, _SALT_BYTE);

    uint8_t y[_PUB_M_BYTE];
    shake256(y, _PUB_M_BYTE, m_salt, mlen + _SALT_BYTE);  // H(M||salt)

    Masked masked_y[_PUB_M_BYTE], masked_vinegar[_V_BYTE];

    unsigned n_attempt = 0;
    while (MAX_ATTEMPT_VINEGAR > n_attempt) {
        uint8_t ctr = n_attempt & 0xff;                    // counter for generating vinegar
        n_attempt++;
        // construct masked_h_vinegar
        memset(masked_h_vinegar, 0, (mlen + _SALT_BYTE + LEN_SKSEED + 1) * N_SHARES);

        memcpy(masked_h_vinegar, message, mlen);
        memcpy(masked_h_vinegar + mlen, salt, _SALT_BYTE);
        for (k = 0; k < N_SHARES; k++)
            memcpy(masked_h_vinegar + mlen + _SALT_BYTE + k * (mlen + _SALT_BYTE + LEN_SKSEED +1), 
                    masked_sk->masked_sk_seed + k * LEN_SKSEED, LEN_SKSEED);

        masked_h_vinegar[mlen + _SALT_BYTE + LEN_SKSEED] = ctr; 
        shake256_masked(vinegar, _V_BYTE, masked_h_vinegar, mlen + _SALT_BYTE + LEN_SKSEED + 1);


        // extract masked_vinegar
        for (k = 0; k < N_SHARES; k++) {
            for (i = 0; i < _V_BYTE; i++) {
                masked_vinegar[i].shares[k] = vinegar[i + k * _V_BYTE];
            }
        }

        // generate linear system
        // matrix
        masked_gf256mat_prod_mtimes(mat_l1, masked_sk->masked_S, _O * _O_BYTE, _V, masked_vinegar);
        // SecMatVec(mat_l1, masked_sk->masked_S, _O * _O_BYTE, _V, masked_vinegar);
        for (k = 0; k < N_SHARES; k++){
            for (i = 0; i < _O * _O_BYTE; i++){
                masked_mat_l1.shares[k][i] = mat_l1[i].shares[k];
            }
        }
        
        // constant
        // Given vinegars, evaluate P1 with the vinegars
        masked_quad_trimat_eval(masked_r_l1_F1, masked_sk->P1, masked_vinegar, _V, _O_BYTE);

        // substract the contribution from vinegar variables
        for (i = 0; i < _O_BYTE; i++) {
            masked_r_l1_F1[i].shares[0] ^= y[i];
        }

        // solve linear system
        unsigned l1_succ = masked_linear_equation_solver(x, masked_mat_l1, masked_r_l1_F1);
        if (!l1_succ) {
            continue;
        }
        break; // successfully solved the linear system
    }

    // failed with limitted attempts
    if ( MAX_ATTEMPT_VINEGAR <= n_attempt ) {
        return -1;
    }

    // secure recombine shares of x_o1
    Masked x_refreshed[_O_BYTE];
    vec_arithmetic_refresh(x_refreshed, x, _O_BYTE);
    vec_combine_arith_shares(x_o1, x_refreshed, _O_BYTE);

    //  w = T^-1 * x
    Masked masked_signature[_PUB_N_BYTE + _SALT_BYTE]; // signature = w || masked_salt
    Masked *w = masked_signature;    // [_PUB_N_BYTE];
    
    // initialize w

    for (i = 0; i < _V_BYTE; i++){
        for (k = 0; k < N_SHARES; k++){
            w[i].shares[k] = masked_vinegar[i].shares[k];
        }
    }


    for (i = 0; i < _O_BYTE; i++){
        for (k = 0; k < N_SHARES; k++){
            w[i + _V_BYTE].shares[k] = x[i].shares[k];
        }
    }
    

    // Computing the O part of T.
    half_masked_gf256matvec_prod(masked_y, masked_sk->masked_O, _V_BYTE, _O, x_o1);
    // masked_gf256mat_prod_mtimes(masked_y, masked_sk->masked_O, _V_BYTE, _O, x);

    masked_gf256v_add(w, masked_y, _V_BYTE);

    // return: masked signature <- w || masked_salt, signature length =  _PUB_N_BYTE + _SALT_BYTE
    for (i = 0; i < _SALT_BYTE; i++) {
            masked_signature[_PUB_N_BYTE + i].shares[0] = salt[i];
    }
    for (k = 1; k < N_SHARES; k++) {
        for (i = 0; i < _SALT_BYTE; i++) {
            masked_signature[_PUB_N_BYTE + i].shares[k] = 0; // other shares are 0
        }
    }

    Masked masked_signature_refreshed[_PUB_N_BYTE + _SALT_BYTE];
    // refresh the masked signature
    vec_arithmetic_refresh(masked_signature_refreshed, masked_signature, _PUB_N_BYTE + _SALT_BYTE);
    // combine the shares of masked signature
    vec_combine_arith_shares(signature, masked_signature_refreshed, _PUB_N_BYTE + _SALT_BYTE);

    return 0;
}


int masked_generate_keypair_pkc( cpk_t *pk, masked_sk_t *masked_sk, const unsigned char *masked_sk_seed ){
    int i, j, k;
    // copy sk_seed
    memcpy(masked_sk->masked_sk_seed, masked_sk_seed, LEN_SKSEED * N_SHARES);

    // pk_seed || O
    size_t bufsize = (size_t) (LEN_PKSEED + _V_BYTE * _O) * (size_t) N_SHARES;
    // unsigned char buf[(LEN_PKSEED + _V_BYTE * _O) * N_SHARES];
    unsigned char *buf = malloc(bufsize);
    // prng for sk
    shake256_masked(buf, LEN_PKSEED + _V_BYTE * _O, masked_sk_seed, LEN_SKSEED);
    uint8_t masked_pk_seed[LEN_PKSEED*N_SHARES]; 
    size_t masked_O_size = (size_t) _V_BYTE * _O * N_SHARES;
    // uint8_t masked_O[_V_BYTE * _O * N_SHARES];
    uint8_t *masked_O = malloc(masked_O_size);
    // extract masked_pk_seed and masked_O from buf
    for (int i = 0; i < LEN_PKSEED; i++){
        for (int k = 0; k < N_SHARES; k++){
            masked_pk_seed[i + k * LEN_PKSEED] = buf[i + k * (LEN_PKSEED + _V_BYTE * _O)];
        }
    }

    for (int i = 0; i < _V_BYTE * _O; i++){
        for (int k = 0; k < N_SHARES; k++){
            masked_O[i + k * _V_BYTE * _O] = buf[LEN_PKSEED + i + k * (LEN_PKSEED + _V_BYTE * _O)];
        }
    }

    // unmasked pk_seed
    uint8_t pk_seed[LEN_PKSEED] = {0};
    for (i = 0; i < LEN_PKSEED; i++){
        for (k = 0; k < N_SHARES; k++){
            pk_seed[i] ^= masked_pk_seed[i + k * LEN_PKSEED];
        }
    }
    // load to pk
    memcpy(pk->pk_seed, pk_seed, LEN_PKSEED);

    // load to masked_sk
    for (i = 0; i < _V_BYTE * _O; i++){
        for (k = 0; k < N_SHARES; k++){
            masked_sk->masked_O[i].shares[k] = masked_O[i + k * _V_BYTE * _O];
        }
    }

    // uint8_t P2[_PK_P2_BYTE];
    size_t P2_size = (size_t) _PK_P2_BYTE;
    uint8_t *P2 = malloc(P2_size);
    // printf("P2 successfully allocated.\n");
    // prng for pk
    prng_publicinputs_t prng1;
    prng_set_publicinputs( &prng1, pk_seed );
    // P1
    prng_gen_publicinputs(&prng1, masked_sk->P1, _PK_P1_BYTE );
    // P2
    prng_gen_publicinputs(&prng1, P2, _PK_P2_BYTE );
    prng_release_publicinputs(&prng1);
    // initialize S
    for (i = 0; i < _PK_P2_BYTE; i++)
            masked_sk->masked_S[i].shares[0] = P2[i];
    
    for (i = 0; i < _PK_P2_BYTE; i++){
        for (k = 1; k < N_SHARES; k++){
            masked_sk->masked_S[i].shares[k] = 0; // other shares are 0
        }
    }
    masked_calculate_F2_P3(masked_sk->masked_S, pk->P3, masked_sk->P1, masked_sk->masked_S, masked_sk->masked_O);
    free(buf);
    free(masked_O);
    free(P2);
    
    return 0;
}
