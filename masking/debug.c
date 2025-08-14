#include "debug.h"
#include "gadgets.h"

void print_masked_bool(Masked* y){
  int t=0;
  printf(" (");
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("%i, ",y->shares[i]);
    t ^= y->shares[i];
  }
  t ^= y->shares[MASKING_ORDER];
  printf("%i) = (", y->shares[MASKING_ORDER]);
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("0x%X, ", y->shares[i]);
  }
  printf("0x%x) = %i\n", y->shares[MASKING_ORDER], t);
}

void print_bitstring(uint8_t * bs, int len){
  for(int i=0; i < len; ++i) printf("%02X ",bs[i]);
  printf("\n");
}

void unmask_bitstring(uint8_t * bs, int len){
  
  unsigned char unmasked_buf[len];

  for(int i=0; i < len; ++i) unmasked_buf[i] = 0;
  for(int k=0; k < MASKING_ORDER+1; ++k){
    for(int i=0; i < len; ++i) unmasked_buf[i] ^= bs[i+k*len];
  }
  for(int i=0; i < len; ++i) printf("%02X ",unmasked_buf[i]);
  printf("\n");
}

void combine_boolean_shares(uint8_t *unmasked_buf, uint8_t *bs, int len)
{
  // unsigned char unmasked_buf[len];

  for(int i=0; i < len; ++i) unmasked_buf[i] = 0;
  for(int k=0; k < MASKING_ORDER+1; ++k){
    for(int i=0; i < len; ++i) unmasked_buf[i] ^= bs[i+k*len];
  }
}

// print column-major matrix
void print_matrix(const uint8_t *matrix, unsigned n_rows, unsigned n_cols) {
    // Iterate over rows.
    for (unsigned i = 0; i < n_rows; i++) {
        // Iterate over columns.
        for (unsigned j = 0; j < n_cols; j++) {
            unsigned index = j * n_rows + i;
            printf("%02X ", matrix[index]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_rowmajor(const uint8_t *matrix, unsigned n_rows, unsigned n_cols) {
    // Iterate over rows.
    for (unsigned i = 0; i < n_rows; i++) {
        // Iterate over columns.
        for (unsigned j = 0; j < n_cols; j++) {
            unsigned index = i * n_cols + j;
            printf("%02X ", matrix[index]);
        }
        printf("\n");
    }
    printf("\n");
}


void print_masked_arith(Masked* x){
  printf(" (");
  int t=0;
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("0x%X, ",x->shares[i]);
    t ^= x->shares[i];
  }
  t ^= x->shares[MASKING_ORDER];
  printf("0x%X) = 0x%X \n", x->shares[MASKING_ORDER], t);

}

void combine_arith_shares(uint8_t *unmasked_vec, uint8_t *masked_vec, unsigned n_A_vec_byte){
    int i, k;
    for (i = 0; i < n_A_vec_byte; i++){
        int t = 0;
        for (k = 0; k < MASKING_ORDER+1; k++){
            t ^= masked_vec[k*n_A_vec_byte + i]; // GF(256) addition
        }
        unmasked_vec[i] = t;
    }
}

void mat_combine_arith_shares(uint8_t *unmasked_mat, const Masked_matrix masked_mat, unsigned n_A_vec_byte, unsigned n_A_width){
    int i, k;
    uint8_t matR[n_A_vec_byte * n_A_width];

    memcpy(matR, masked_mat.shares[0], n_A_vec_byte * n_A_width);
    // for (i = 0; i < n_A_vec_byte * n_A_width; i++)
    
    for (k = 1; k < N_SHARES; k++){
        // uint8_t temp_mat[n_A_vec_byte * n_A_width];
        // for (i = 0; i < n_A_vec_byte * n_A_width; i++)
        //   temp_mat[i] = masked_mat[i].shares[k];
        gf256mat_add(matR, matR, masked_mat.shares[k], n_A_vec_byte, n_A_width);
    }
    memcpy(unmasked_mat, matR, n_A_vec_byte * n_A_width);
}

void vec_combine_arith_shares(uint8_t *unmasked_vec, const Masked *masked_vec, unsigned n_A_vec_byte){
    int i, k;
    for (i = 0; i < n_A_vec_byte; i++){
        int t = 0;
        for (k = 0; k < MASKING_ORDER+1; k++){
            t ^= masked_vec[i].shares[k]; // GF(256) addition
        }
        unmasked_vec[i] = t;
    }
}