#include "matrix_arith.h"
#include "debug.h"
#include "gadgets.h"
#include "src/ov_blas.h"

// x = A * b with ISW style
void masked_gf256mat_prod(Masked *c, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked *b){
    int s, i, j, k;

    // 1) flat share‐planes for c and b
    uint8_t c_plane[N_SHARES][n_A_vec_byte];
    uint8_t b_plane[N_SHARES][n_A_width];

    // 2) unpack b → b_plane
    for (s = 0; s < N_SHARES; s++) {
        for (k = 0; k < n_A_width; k++) {
            b_plane[s][k] = b[k].shares[s];
        }
    }

    // 3) sharewise matrix–vector product
    for (s = 0; s < N_SHARES; s++) {
        // matA.shares[s] is a flat matrix (row‐major, row count = n_A_vec_byte,
        // col count = n_A_width), b_plane[s] is flat vector length n_A_width
        // result goes into c_plane[s], length n_A_vec_byte
        gfmat_prod(
          c_plane[s],
          matA.shares[s],
          n_A_vec_byte,
          n_A_width,
          b_plane[s]
        );
    }

    // 4) ISW cross‐terms
    //    for each pair of distinct shares (i<j), add a masked cross‐term
    uint8_t tmp[n_A_vec_byte];
    uint8_t r0 [n_A_vec_byte];
    uint8_t r1 [n_A_vec_byte];

    for (i = 0; i < N_SHARES - 1; i++) {
        for (j = i + 1; j < N_SHARES; j++) {
            // sample fresh mask r0
            for (k = 0; k < n_A_vec_byte; k++) {
                r0[k] = rand8();
            }

            // tmp ← A_i * b_j
            gfmat_prod(tmp, matA.shares[i], n_A_vec_byte, n_A_width, b_plane[j]);
            // r1 = r0 ⊕ tmp
            for (k = 0; k < n_A_vec_byte; k++) {
                r1[k] = r0[k] ^ tmp[k];
            }

            // tmp ← A_j * b_i
            gfmat_prod(tmp, matA.shares[j], n_A_vec_byte, n_A_width, b_plane[i]);
            // r1 ⊕= tmp
            for (k = 0; k < n_A_vec_byte; k++) {
                r1[k] ^= tmp[k];
            }

            // fold into c_plane
            for (k = 0; k < n_A_vec_byte; k++) {
                c_plane[i][k] ^= r0[k];  // share i += r0
                c_plane[j][k] ^= r1[k];  // share j += (tmp⊕r0⊕tmp')
            }
        }
    }

    // 5) repack c_plane → interleaved Masked *c
    for (k = 0; k < n_A_vec_byte; k++) {
        for (s = 0; s < N_SHARES; s++) {
            c[k].shares[s] = c_plane[s][k];
        }
    }
}


// used in generating the linear system in Sign.
void masked_gf256mat_prod_mtimes(Masked *c, const Masked* matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked *b){
    int s, i, j, k;
    size_t planeA_size = (size_t)N_SHARES * n_A_vec_byte * n_A_width;
    size_t planeb_size = (size_t)N_SHARES * n_A_width;
    size_t planec_size = (size_t)N_SHARES * n_A_vec_byte;
    
    // 1) heap‐allocate flat planes:
    uint8_t *A_plane = malloc(planeA_size);
    uint8_t *b_plane = malloc(planeb_size);
    uint8_t *c_plane = malloc(planec_size);
    uint8_t *tmp     = malloc(n_A_vec_byte);
    uint8_t *r0      = malloc(n_A_vec_byte);
    uint8_t *r1      = malloc(n_A_vec_byte);
    if (!A_plane || !b_plane || !c_plane || !tmp || !r0 || !r1) {
        perror("malloc");
        exit(1);
    }
    
    // Helpers to index A_plane[s][k] etc:
    #define A_PLANE(s,k)  A_plane[((s)*(n_A_vec_byte*n_A_width)) + (k)]
    #define b_PLANE(s,k)  b_plane[((s)*n_A_width) + (k)]
    #define c_PLANE(s,k)  c_plane[((s)*n_A_vec_byte) + (k)]

    // 2) unpack matA to A_plane
    for (s = 0; s < N_SHARES; s++) {
        for (k = 0; k < n_A_vec_byte * n_A_width; k++) {
            A_PLANE(s,k) = matA[k].shares[s];
        }
    }
    // 3) unpack b to b_plane
    for (s = 0; s < N_SHARES; s++) {
        for (k = 0; k < n_A_width; k++) {
            b_PLANE(s,k) = b[k].shares[s];
        }
    }

    // 4) sharewise mat‐vec
    for (s = 0; s < N_SHARES; s++) {
        gfmat_prod(
          &c_PLANE(s,0),          // pointer into share-s result plane
          &A_PLANE(s,0),          // pointer into share-s matrix plane
          n_A_vec_byte,
          n_A_width,
          &b_PLANE(s,0)           // pointer into share-s input vector
        );
    }

    // 5) ISW cross‐terms
    for (i = 0; i < N_SHARES - 1; i++) {
        for (j = i + 1; j < N_SHARES; j++) {
            for (k = 0; k < n_A_vec_byte; k++) {
                r0[k] = rand8();        // random mask
            }
            gfmat_prod(tmp, &A_PLANE(i,0), n_A_vec_byte, n_A_width, &b_PLANE(j,0));
            for (k = 0; k < n_A_vec_byte; k++) r1[k] = r0[k] ^ tmp[k];
            gfmat_prod(tmp, &A_PLANE(j,0), n_A_vec_byte, n_A_width, &b_PLANE(i,0));
            for (k = 0; k < n_A_vec_byte; k++) r1[k] ^= tmp[k];
            for (k = 0; k < n_A_vec_byte; k++) {
                c_PLANE(i,k) ^= r0[k];
                c_PLANE(j,k) ^= r1[k];
            }
        }
    }

    // 6) repack c_plane to interleaved c
    for (k = 0; k < n_A_vec_byte; k++) {
        for (s = 0; s < N_SHARES; s++) {
            c[k].shares[s] = c_PLANE(s,k);
        }
    }

    // 7) clean up
    free(A_plane); free(b_plane);
    free(c_plane); free(tmp);
    free(     r0); free(r1);

    #undef A_PLANE
    #undef b_PLANE
    #undef c_PLANE

}

// share-wise matrix-vector product
// c = A * b, where A is a masked matrix and b is a vector
void half_masked_gf256matvec_prod(Masked *c, const Masked* matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    int i, k, j;

    // Initialize c to zero
    for (k = 0; k < N_SHARES; k++) {
        for (i = 0; i < n_A_vec_byte; i++) {
            c[i].shares[k] = 0;
        }
    }

    for (k = 0; k < N_SHARES; k++) { // Iterate over each share
        for (i = 0; i < n_A_width; i++) { // Iterate over matrix columns
            const Masked *current_column = matA + i * n_A_vec_byte;

            // Extract the k-th share of the entire column
            uint8_t column_share[n_A_vec_byte];
            for (j = 0; j < n_A_vec_byte; j++) {
                column_share[j] = current_column[j].shares[k];
            }

            // Use gf256v_madd to multiply column by scalar and add to c
            uint8_t c_temp[n_A_vec_byte];
            for (j = 0; j < n_A_vec_byte; j++) {
                c_temp[j] = c[j].shares[k];
            }
            
            gf256v_madd(c_temp, column_share, b[i], n_A_vec_byte);
            
            for (j = 0; j < n_A_vec_byte; j++) {
                c[j].shares[k] = c_temp[j];
            }
        }
    }
}


void masked_gf256mat_prod_unbatched(Masked *c, const Masked_matrix matA, unsigned n_A_vec_byte, unsigned n_A_width, const Masked *b){
    int i, k, j;
    // initialize c to zero
    for (k = 0; k < N_SHARES; k++){
        for (i = 0; i < n_A_width; i++){
            c[i].shares[k] = 0;
        }
    }
    for (i = 0; i < n_A_width; i++){
        Masked temp_col[n_A_vec_byte];
        for (j = 0; j < n_A_vec_byte; j++){
            for (k = 0; k < N_SHARES; k++){
                temp_col[j].shares[k] = matA.shares[k][i * n_A_vec_byte + j];
            }
        }
        masked_gf256_madd(c, temp_col, b[i], n_A_vec_byte);
    }

}

// y = vT * P1 * v
void masked_quad_trimat_eval(Masked *y, const uint8_t *trimat, Masked *x, unsigned dim, unsigned size_batch){
    int i, j, k;
    // initialize y to zero
    for (k = 0; k < N_SHARES; k++){
        for (i = 0; i < size_batch; i++)
            y[i].shares[k] = 0;
    }
    /// assert (dim <= 256) ;
    Masked tmp[size_batch];
    
    // refresh the input vector
    Masked s[dim];
    vec_arithmetic_refresh(s, x, dim); // t-SNI refresh

    for (i = 0; i < dim - 15; i++){
        half_masked_gf256mat_prod(tmp, trimat, size_batch, dim - i, x + i);
        masked_gf256v_madd(y, tmp, s[i], size_batch);
        trimat += (dim - i) * size_batch;
    }

    Masked quad_terms[128];
    // Initialize quad_terms to zero
    for (i = 0; i < 128; i++) {
        for (k = 0; k < N_SHARES; k++) {
            quad_terms[i].shares[k] = 0;
        }
    }
    
    int idx = 0;
    for (i = dim - 15; i < dim; i++) {
        masked_gf256_madd(quad_terms + idx, x + i, s[i], dim - i);
        idx += dim - i;
    }

    half_masked_gf256mat_prod(tmp, trimat, size_batch, 120, quad_terms);
    masked_gf256v_add(y, tmp, size_batch);
}

// bC += btriA * B, bC and B are masked, btriA is unmasked
void masked_batch_trimat_madd(Masked *bC, const unsigned char *btriA, 
                            const Masked *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch ){
    int i, j;
    Masked tmp_c[96]; // MAX_O_BYTE

    // access fixed positions of destination matrix C
    unsigned Aheight = Bheight;
    for (i = 0; i < Aheight; i++){
        for (j = 0; j < Bwidth; j++){
            half_masked_gf256mat_prod( tmp_c, btriA, size_batch, Aheight - i, B + j * size_Bcolvec + i );
            masked_gf256v_add(bC, tmp_c, size_batch);
            bC += size_batch;
        }
        btriA += size_batch * (Aheight - i);
    }
}

// bC += btriA^Tr * B, bC and B are masked, btriA is unmasked
void masked_batch_trimatTr_madd(Masked *bC, const unsigned char *btriA, 
                            const Masked *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch ){
    int i, j;
    Masked tmp_c[96]; // MAX_O_BYTE
    uint8_t tmp_Arow[148 * 96]; // MAX_V_BYTE * MAX_O_BYTE

    // access fixed positions of destination matrix C
    unsigned Aheight = Bheight;
    for (i = 0; i < Aheight; i++){
        const uint8_t *ptr = btriA + i*size_batch;
        for (j = 0; j < i; j++){
            memcpy( tmp_Arow + j * size_batch, ptr, size_batch );
            ptr += (Aheight - j - 1) * size_batch;
        }
        memcpy( tmp_Arow + i * size_batch, ptr, size_batch );

        for (j = 0; j < Bwidth; j++){
            half_masked_gf256mat_prod(tmp_c, tmp_Arow, size_batch, i + 1, B + j * size_Bcolvec);
            masked_gf256v_add(bC, tmp_c, size_batch);
            bC += size_batch;
        }
    }
}

// bC = A^Tr * bB, all masked
void masked_batch_upper_matTr_x_mat( Masked *bC,
                               const Masked *A_to_tr, unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
                               const Masked *bB, unsigned Bwidth, unsigned size_batch ){
    int i, j, k;
    Masked row[96 * 96]; // MAX_O_BYTE * MAX_O_BYTE
    unsigned Atr_height = Awidth;
    unsigned Atr_width = Aheight;
    
    for (i = 0; i < Atr_height; i++){
        masked_gf256mat_prod_mtimes( row, bB, Bwidth * size_batch, Atr_width, A_to_tr + size_Acolvec * i );
        Masked *ptr = bC + i * size_batch;
        for (j = 0; j < i; j++){
            masked_gf256v_add(ptr, row + size_batch * j, size_batch);
            ptr += (Bwidth - j - 1) * size_batch;
        }
        
        for (k = 0; k < N_SHARES; k++){
            for (j = 0; j < size_batch * (Bwidth - i); j++){
                ptr[j].shares[k] = row[size_batch * i + j].shares[k];
            }
        } 
    }
}


void masked_calculate_F2_P3(Masked *masked_S, unsigned char *P3, const unsigned char *P1, const Masked *P2, Masked *masked_O){
    size_t P3_size = (size_t) _PK_P3_BYTE;
    // Masked masked_P3[_PK_P3_BYTE];
    Masked *masked_P3 = malloc(P3_size * sizeof(Masked));
    // S = P1 * O + P2
    masked_batch_trimat_madd(masked_S, P1, masked_O, _V, _V_BYTE, _O, _O_BYTE);
    // refresh O
    Masked masked_O_refreshed[_V_BYTE * _O];
    vec_arithmetic_refresh(masked_O_refreshed, masked_O, _V_BYTE * _O);
    // P3 = Upper(O^Tr * (P1*O + P2))
    masked_batch_upper_matTr_x_mat(masked_P3, masked_O_refreshed, _V, _V_BYTE, _O, masked_S, _O, _O_BYTE);
    // S += P1 * O + P2 + P1^Tr * O
    masked_batch_trimatTr_madd(masked_S, P1, masked_O, _V, _V_BYTE, _O, _O_BYTE);

    // recombine shares of P3
    // Masked masked_P3_refreshed[_PK_P3_BYTE];
    Masked *masked_P3_refreshed = malloc(P3_size * sizeof(Masked));
    vec_arithmetic_refresh(masked_P3_refreshed, masked_P3, _PK_P3_BYTE);
    // combine shares of P3
    vec_combine_arith_shares(P3, masked_P3_refreshed, _PK_P3_BYTE);
    
    free(masked_P3);
    free(masked_P3_refreshed);
}