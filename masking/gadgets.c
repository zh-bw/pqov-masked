#include <stdint.h>
#include <string.h>
#include "gadgets.h"
#include "debug.h"

/* ================================================== REFRESHING =================================================== */
// GF(256)
// O(n* (3m^2))
void mat_linear_arithmetic_refresh(Masked_matrix *matA, unsigned n_A_vec_byte, unsigned n_A_width){
    uint8_t matR[n_A_vec_byte * n_A_width];
    int i;
    for (i = 0; i < N_SHARES - 1; i++){
        gf256mat_gen(matR, n_A_width); // sample random matrix R_i
        gf256mat_add(matA->shares[i], matA->shares[i], matR, n_A_vec_byte, n_A_width); // A_i = A_i + R_i
        gf256mat_add(matA->shares[N_SHARES - 1], matA->shares[N_SHARES - 1], matR, n_A_vec_byte, n_A_width); // A_n = A_n - R_i
    }
}

void mat_refresh_masks(Masked_matrix *matC, Masked_matrix *matA, unsigned n_A_vec_byte, unsigned n_A_width){
    uint8_t matR[n_A_vec_byte * n_A_width];
    int i, j;
    for (i = 0; i < N_SHARES; i++) memcpy(matC->shares[i], matA->shares[i], n_A_vec_byte * n_A_width);
    
    for (i = 0; i < N_SHARES ; i++){
        for (j = i + 1; j < N_SHARES; j++)
        {
            gf256mat_gen(matR, n_A_width); // sample random matrix R_i
            gf256mat_add(matC->shares[i], matC->shares[i], matR, n_A_vec_byte, n_A_width); // C_i = C_i + R_i
            gf256mat_add(matC->shares[j], matC->shares[j], matR, n_A_vec_byte, n_A_width); // C_j = C_j - R_i
        }
    }
}


// O(n)
void gf256_linear_arithmetic_refresh(Masked *a){
    uint8_t r;
    int i;
    for (i = 0; i < N_SHARES - 1; i++){
        r = rand8(); // sample random r_i
        a->shares[i] ^= r; // a_i = a_i + r_i
        a->shares[N_SHARES - 1] ^= r; // a_n = a_n - r_i
    }
}

void gf256_linear_arithmetic_refresh_u32(Masked_u32 *a){
    uint32_t r;
    int i;
    for (i = 0; i < N_SHARES - 1; i++){
        r = rand32(); // sample random r_i
        a->shares[i] ^= r; // a_i = a_i + r_i
        a->shares[N_SHARES - 1] ^= r; // a_n = a_n - r_i
    }
}

// O(n^2)
void gf256_arithmetic_refresh(Masked *z, Masked *a){
    uint8_t r;
    int i, j;
    // copy a to z
    for (i = 0; i < N_SHARES; i++)
        z->shares[i] = a->shares[i];
    
    for (i = 0; i < N_SHARES; i++){
        for (j = i + 1; j < N_SHARES; j++)
        {
            r = rand8();
            z->shares[i] ^= r;
            z->shares[j] ^= r;
        }
    }
}

// O(n^2)
void vec_arithmetic_refresh(Masked *r, Masked *a, unsigned n_vec_byte){
    int i;
    
    // Apply arithmetic refresh to each element of the vector
    for (i = 0; i < n_vec_byte; i++) {
        gf256_arithmetic_refresh(&r[i], &a[i]);
    }
}

void vec_linear_arithmetic_refresh(Masked *r, unsigned n_vec_byte){
    // Apply linear arithmetic refresh to each element of the vector
    for (int i = 0; i < n_vec_byte; i++) {
        gf256_linear_arithmetic_refresh(&r[i]);
    }
}


void half_masked_gf256mat_prod(Masked *masked_c,
                               const uint8_t *matA,
                               unsigned n_A_vec_byte,
                               unsigned n_A_width,
                               const Masked *masked_b)
{
    // Temporary buffers for one share at a time
    uint8_t b_share[n_A_width];
    uint8_t c_share[n_A_vec_byte];

    // For each share k
    for (unsigned k = 0; k < N_SHARES; k++) {
        // extract the k‑th share of b
        for (unsigned j = 0; j < n_A_width; j++) {
            b_share[j] = masked_b[j].shares[k];
        }

        // compute unmasked product: c_share = A * b_share
        gf256mat_prod(c_share, matA, n_A_vec_byte, n_A_width, b_share);

        // scatter c_share back into the k‑th share of masked_c
        for (unsigned i = 0; i < n_A_vec_byte; i++) {
            masked_c[i].shares[k] = c_share[i];
        }
    }
}

// =================================================== main gadgets ===================================================

static inline 
void A2M_rankcheck(uint8_t *matC, Masked_matrix* matR, Masked_matrix* matU, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width){
    Masked_matrix matB;
    uint8_t temp[n_A_vec_byte * n_A_width], temp2[n_A_vec_byte * n_A_width];
    uint8_t R[n_A_vec_byte * n_A_width];
    uint8_t U[n_A_vec_byte * n_A_width];
    unsigned R_is_fullrank = 0, U_is_fullrank = 0;
    int i, k;
    // (B_1, ..., B_n) <- (A_1, ..., A_n)
    for (k = 0; k < N_SHARES; k++)
        memcpy(matB.shares[k], matA.shares[k], _O_BYTE * _O_BYTE);

    for (i = 0; i < N_SHARES; i++){
        while (!R_is_fullrank)
        {
            gf256mat_gen(R, n_A_width); // sample invertible matrix R_i
            R_is_fullrank = check_fullrank(R, n_A_width, n_A_width );
        }

        while (!U_is_fullrank)
        {
            gf256mat_gen(U, n_A_width); // sample invertible matrix U_i
            U_is_fullrank = check_fullrank(U, n_A_width, n_A_width );
        }


        memcpy(matR->shares[i], R, n_A_vec_byte * n_A_width);
        memcpy(matU->shares[i], U, n_A_vec_byte * n_A_width);
        // (B_1, ..., B_n) <- U_i * (B_1, ..., B_n) * R_i
        for (k = 0; k < N_SHARES; k++)
        {
            // gf256mat_mul_square(temp, matB.shares[k], R, n_A_width);
            gf256mat_mul_square(temp, matB.shares[k], R, n_A_width);
            gf256mat_mul_square(temp2, U, temp, n_A_width);
            // gf256mat_mul_square_tri(temp, matB.shares[k], R, n_A_width);
            memcpy(matB.shares[k], temp2, n_A_vec_byte * n_A_width);
        }
        // linear_arithmetic_refresh of matB
        mat_linear_arithmetic_refresh(&matB, n_A_vec_byte, n_A_width);
    }
    // C <- (B_1 + ...+ B_n)
    for (k = 0; k < N_SHARES; k++)
        gf256mat_add(matC, matC, matB.shares[k], n_A_vec_byte, n_A_width);
}

// for benchmark test
void masked_linear_equation_solver_rankcheck(Masked *masked_vec_r, const Masked_matrix masked_matA, 
                                   const Masked *masked_vec_b)
{
    Masked_matrix matR, matU;
    uint8_t matC[_O_BYTE * _O_BYTE];
    uint8_t matC_inv[_O_BYTE * _O_BYTE];
    memset(matC, 0, _O_BYTE * _O_BYTE);
    unsigned ok;

    
    // A2M_rankcheck(matC, &matR, masked_matA, _O_BYTE, _O_BYTE);
    A2M_rankcheck(matC, &matR, &matU, masked_matA, _O_BYTE, _O_BYTE);
    ok = gf256mat_inv_nonconst(matC, matC_inv, _O_BYTE);

    // uint8_t masked_vec_temp[_O_BYTE * N_SHARES];
    
    half_masked_gf256mat_prod(masked_vec_r, matU.shares[0], _O_BYTE, _O_BYTE, masked_vec_b);
    vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);

    for (int i = 1; i < N_SHARES; i++){
        half_masked_gf256mat_prod(masked_vec_r, matU.shares[i], _O_BYTE, _O_BYTE, masked_vec_r);
        vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);
    }

    
    // O(nm^2)
    half_masked_gf256mat_prod(masked_vec_r, matC_inv, _O_BYTE, _O_BYTE, masked_vec_r);
    vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);

    // O(n^2 * m^2)
    for (int i = N_SHARES - 1; i >= 0; i--){
        half_masked_gf256mat_prod(masked_vec_r, matR.shares[i], _O_BYTE, _O_BYTE, masked_vec_r);
        vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);
    }

}



// implementated with rankcheck, used in masked ov sign
unsigned masked_linear_equation_solver(Masked *masked_vec_r, const Masked_matrix masked_matA, 
                                   const Masked *masked_vec_b)
{
    Masked_matrix matR, matU;
    uint8_t matC[_O_BYTE * _O_BYTE];
    uint8_t matC_inv[_O_BYTE * _O_BYTE];
    memset(matC, 0, _O_BYTE * _O_BYTE);
    unsigned ok;

    
    A2M_rankcheck(matC, &matR, &matU, masked_matA, _O_BYTE, _O_BYTE);
    ok = gf256mat_inv_nonconst(matC, matC_inv, _O_BYTE);

    if (!ok)
    {
        // printf("Matrix is not invertible, retrying...\n");
        return ok;
    }

    half_masked_gf256mat_prod(masked_vec_r, matU.shares[0], _O_BYTE, _O_BYTE, masked_vec_b);
    vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);

    for (int i = 1; i < N_SHARES; i++){
        half_masked_gf256mat_prod(masked_vec_r, matU.shares[i], _O_BYTE, _O_BYTE, masked_vec_r);
        vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);
    }

    
    // O(nm^2)
    half_masked_gf256mat_prod(masked_vec_r, matC_inv, _O_BYTE, _O_BYTE, masked_vec_r);
    vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);

    // O(n^2 * m^2)
    for (int i = N_SHARES - 1; i >= 0; i--){
        half_masked_gf256mat_prod(masked_vec_r, matR.shares[i], _O_BYTE, _O_BYTE, masked_vec_r);
        vec_linear_arithmetic_refresh(masked_vec_r, _O_BYTE);
    }


    return ok;
}


// Randomness: O(n*m^2 + n(n-1)/2 * m^2)
// Complexity: O(m^2 * (3.5n^2 + 0.5n)+ m^3 * (n^2 + 2) + m * (2.5n^2 - 2.5n))
// for benchmark test
void masked_linear_equation_solver_alternative(Masked *masked_vec_r, const Masked_matrix masked_matA, 
                                    const Masked *masked_vec_b)
{
    Masked_matrix matB, masked_matC, masked_matU, temp;
    unsigned ok = 0;
    uint8_t unmasked_matB[_O_BYTE * _O_BYTE] = {0}, unmasked_matC [_O_BYTE * _O_BYTE], unmasked_matB_inv[_O_BYTE * _O_BYTE];
    Masked masked_vec_temp[_O_BYTE];

    uint8_t repeat = 0;
    while (!ok){
        repeat++;
        // generate random matrix C and U
        for (int k = 0; k < N_SHARES; k++){
            gf256mat_gen(masked_matC.shares[k], _O_BYTE);
            gf256mat_gen(masked_matU.shares[k], _O_BYTE);
        }
        masked_gf256mat_mul_square(&matB, masked_matA, masked_matC);
        masked_gf256mat_mul_square(&matB, masked_matU, matB);
        mat_refresh_masks(&temp, &matB, _O_BYTE, _O_BYTE);


        for (int k = 0; k < N_SHARES; k++)
            gf256mat_add(unmasked_matB, unmasked_matB, temp.shares[k], _O_BYTE, _O_BYTE);
        
        
        ok = gf256mat_inv(unmasked_matB, unmasked_matB_inv, _O_BYTE);
        if (repeat > 16){
            printf("Matrix is not invertible after 16 retries, exiting...\n");
            return;
        }

    }
    masked_gf256mat_prod(masked_vec_temp, masked_matU, _O_BYTE, _O_BYTE, masked_vec_b);

    half_masked_gf256mat_prod(masked_vec_temp, unmasked_matB_inv, _O_BYTE, _O_BYTE, masked_vec_temp);

    masked_gf256mat_prod(masked_vec_r, masked_matC, _O_BYTE, _O_BYTE, masked_vec_temp);
}



// improved ISW from Cor14, with space complexity O(m^2 * n)
// Complexity: O(n^2 * (2.5m^2 + m^3) - 2.5n*m^2) 
void masked_gf256mat_mul_square(Masked_matrix *matC, const Masked_matrix matA, const Masked_matrix matB){
    int i,j;
    uint8_t matR[2][_O_BYTE * _O_BYTE];
    uint8_t temp[_O_BYTE * _O_BYTE];
    // printf("square matrix multiplication with ISW\n");
    for (i = 0; i < N_SHARES; i++)
        gf256mat_mul_square(matC->shares[i], matA.shares[i], matB.shares[i], _O_BYTE); // C_i = A_i * B_i
    for (i = 0; i < N_SHARES - 1; i++){
        for (j = i+1; j < N_SHARES; j++){
            gf256mat_mul_square(temp, matA.shares[i], matB.shares[j], _O_BYTE);
            gf256mat_gen(matR[0], _O_BYTE); // sample random matrix R_i
            gf256mat_add(matR[1], matR[0], temp, _O_BYTE, _O_BYTE);
            gf256mat_mul_square(temp, matA.shares[j], matB.shares[i], _O_BYTE);
            gf256mat_add(matR[1], matR[1], temp, _O_BYTE, _O_BYTE);

            gf256mat_add(matC->shares[i], matC->shares[i], matR[0], _O_BYTE, _O_BYTE); // C_i = C_i + R0_i
            gf256mat_add(matC->shares[j], matC->shares[j], matR[1], _O_BYTE, _O_BYTE); // C_j = C_j - R1_i
        }
    }
}

// original ISW, with space complexity O(n^2)
void masked_gf256mat_mul_square_ISW(Masked_matrix *matC, const Masked_matrix matA, const Masked_matrix matB){
    int i, j;
    uint8_t matR[N_SHARES][N_SHARES][_O_BYTE * _O_BYTE];
    uint8_t temp[_O_BYTE * _O_BYTE];

    for (i = 0; i < N_SHARES; i++){
        for (j = i + 1; j < N_SHARES; j++){
            gf256mat_gen(matR[i][j], _O_BYTE); // sample random matrix R_ij
            gf256mat_mul_square(temp, matA.shares[i], matB.shares[j], _O_BYTE);
            gf256mat_add(matR[j][i], matR[i][j], temp, _O_BYTE, _O_BYTE);
            gf256mat_mul_square(temp, matA.shares[j], matB.shares[i], _O_BYTE);
            gf256mat_add(matR[j][i], matR[j][i], temp, _O_BYTE, _O_BYTE);
        }
    }
    for (i = 0; i < N_SHARES; i++){
        gf256mat_mul_square(matC->shares[i], matA.shares[i], matB.shares[i], _O_BYTE); // C_i = A_i * B_i
        for (j = 0; j < N_SHARES; j++){
            if (i != j){
                gf256mat_add(matC->shares[i], matC->shares[i], matR[i][j], _O_BYTE, _O_BYTE); // C_i = A_i * B_i + \sum{R_ij}
            }
        }
    }
}


// implemented from masked mat-vec mult
void masked_gf256mat_mul_square_lowlevel(Masked_matrix *matC, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked_matrix matB){
    int i, k;
    Masked r[n_A_vec_byte], col[n_A_vec_byte];
    for (int j = 0; j < n_A_width; j++){
        // extract the j-th column of B
        for (i = 0; i < n_A_vec_byte; i++){
            for (k = 0; k < N_SHARES; k++)
                col[i].shares[k] = matB.shares[k][i + n_A_vec_byte*j];
            }
        masked_gf256mat_prod(&r, matA, n_A_vec_byte, n_A_width, col);
        // store the result in the j-th column of C
        for (i = 0; i < n_A_vec_byte; i++){
            for (k = 0; k < N_SHARES; k++)
                matC->shares[k][i + n_A_vec_byte*j] = r[i].shares[k];
        }
    }
}




