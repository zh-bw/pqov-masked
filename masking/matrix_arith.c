#include <stdint.h>
#include <string.h>
#include "matrix_arith.h"
#include "random.h"
#include "debug.h"

///////////  matrix-matrix  multiplications  ////////////////////////////////

void gf256mat_mul_square(uint8_t *matC, const uint8_t *matA, const uint8_t *matB, unsigned width) {
    for (int i = 0; i < width; i++)
        gf256mat_prod(matC + i * width, matA, width, width, matB + i * width);
}

void gf256mat_add(uint8_t *r, const uint8_t *mat_A, const uint8_t *mat_B, 
                  unsigned n_A_vec_byte, unsigned n_A_width) {
    for (unsigned col = 0; col < n_A_width; col++) {
        unsigned offset = col * n_A_vec_byte;
        // Copy column from mat_A into r
        memcpy(r + offset, mat_A + offset, n_A_vec_byte);
        // Add the corresponding column from mat_B (GF(256) addition is XOR)
        gf256v_add(r + offset, mat_B + offset, n_A_vec_byte);
    }
}



//////////////////    Gaussian elimination + Back substitution for matrix inversion  //////////////////

static inline
unsigned gf256mat_gauss_elim_row_echolen( uint8_t *mat, unsigned h, unsigned w ) {
    unsigned r8 = 1;

    for (unsigned i = 0; i < h; i++) {
        uint8_t *ai = mat + w * i;
        //unsigned i_start = i-(i&(_BLAS_UNIT_LEN_-1));
        unsigned i_start = i;
        for (unsigned j = i + 1; j < h; j++) {
            uint8_t *aj = mat + w * j;
            gf256v_conditional_add( ai + i_start, !gf256_is_nonzero(ai[i]), aj + i_start, w - i_start );
        }
        r8 &= gf256_is_nonzero(ai[i]);
        uint8_t pivot = ai[i];
        pivot = gf256_inv( pivot );
        gf256v_mul_scalar( ai + i_start, pivot, w - i_start );
        for (unsigned j = i + 1; j < h; j++) {
            uint8_t *aj = mat + w * j;
            gf256v_madd( aj + i_start, ai + i_start, aj[i], w - i_start );
        }
    }
    return r8;
}




// Back Substitution on an Augmented Matrix (row-major)
// Assumes matrix is of size n x (2*n) and in row-echelon form.
static void gf256mat_back_substitute_augmented(uint8_t *aug, unsigned n) {
    unsigned width = 2 * n;
    for (int i = n - 1; i > 0; i--) {
        for (int j = 0; j < i; j++) {
            uint8_t factor = aug[j * width + i];
            for (unsigned k = i; k < width; k++) {
                aug[j * width + k] ^= gf256_mul(factor, aug[i * width + k]);
            }
        }
    }
}

// O(2m^3)
unsigned gf256mat_inv(uint8_t *A, uint8_t *A_inv, unsigned n) {
    const unsigned MAX_H = 96;
    uint8_t aug[MAX_H * (2 * MAX_H)];  // Augmented matrix buffer; ensure n <= MAX_H.
    unsigned height = n;
    unsigned width = 2 * n;
    
    // Build the augmented matrix [A | I] in row-major order.
    // Transpose A from column-major to row-major and build I.
    for (unsigned i = 0; i < height; i++) {
        uint8_t *row = aug + i * width;
        // Left half: transpose of A.
        for (unsigned j = 0; j < n; j++) {
            row[j] = A[j * n + i];  // A(i,j) = A[j*n + i] since A is in column-major.
        }
        // Right half: identity matrix.
        for (unsigned j = 0; j < n; j++) {
            row[n + j] = (i == j) ? 1 : 0;
        }
    }
    
    // Forward elimination: convert augmented matrix to row-echelon form.
    unsigned ok = gf256mat_gauss_elim_row_echolen(aug, height, width);
    if (!ok) {
        // A is singular.
        return 0;
    }
    
    // Back substitution: clear above-diagonal entries.
    gf256mat_back_substitute_augmented(aug, n);
    
    // Extract A_inv from aug and convert it back to column-major order.
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            A_inv[j * n + i] = aug[i * width + n + j];
        }
    }
    return 1;
}

void gf256mat_gen_upper(uint8_t *mat, unsigned n) {
    for (unsigned j = 0; j < n; j++) {
        for (unsigned i = 0; i < n; i++) {
            if (i > j) {
                // Below the diagonal: set to 0.
                mat[i + j * n] = 0;
            } 
            else if (i == j) {
                // Diagonal: choose a random nonzero value.
                uint8_t val = 0;
                while (val == 0) {
                    val = rand8();
                }
                mat[i + j * n] = val;
            } 
            else { // i < j, above the diagonal.
                // Any random value (could be zero).
                mat[i + j * n] = rand8();
            }
        }
    }
}

void gf256mat_gen(uint8_t *mat, unsigned n) {
    for (unsigned j = 0; j < n; j++) {
        for (unsigned i = 0; i < n; i++) {
            mat[i + j * n] = rand8();
            // randombytes(&mat[i + j * n], 1);
        }
    }
}

// check if a col-major matrix is full rank
unsigned check_fullrank(const uint8_t *mat, unsigned height, unsigned width) {
    uint8_t rowmajor[96*96]; // max size
    unsigned r8;
    int i,j;
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            rowmajor[i * width + j] = mat[j * height + i];
        }
    }

    r8 = gf256mat_gauss_elim_row_echolen(rowmajor, height, width);
    return r8;
}