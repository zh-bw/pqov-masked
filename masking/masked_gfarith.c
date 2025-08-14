#include "masked_gfarith.h"
#ifdef COUNT
uint64_t count_op = 0;
#endif
static inline uint32_t _gf256v_mul_u32_u32(
    uint32_t a0, uint32_t a1, uint32_t a2, uint32_t a3,
    uint32_t a4, uint32_t a5, uint32_t a6, uint32_t a7,
    uint32_t b0, uint32_t b1, uint32_t b2, uint32_t b3,
    uint32_t b4, uint32_t b5, uint32_t b6, uint32_t b7)
{
    uint32_t c0  = a0 & b0;
    uint32_t c1  = (a0 & b1) ^ (a1 & b0);
    uint32_t c2  = (a0 & b2) ^ (a1 & b1) ^ (a2 & b0);
    uint32_t c3  = (a0 & b3) ^ (a1 & b2) ^ (a2 & b1) ^ (a3 & b0);
    uint32_t c4  = (a0 & b4) ^ (a1 & b3) ^ (a2 & b2) ^ (a3 & b1) ^ (a4 & b0);
    uint32_t c5  = (a0 & b5) ^ (a1 & b4) ^ (a2 & b3) ^ (a3 & b2) ^ (a4 & b1) ^ (a5 & b0);
    uint32_t c6  = (a0 & b6) ^ (a1 & b5) ^ (a2 & b4) ^ (a3 & b3) ^ (a4 & b2) ^ (a5 & b1) ^ (a6 & b0);
    uint32_t c7  = (a0 & b7) ^ (a1 & b6) ^ (a2 & b5) ^ (a3 & b4) ^ (a4 & b3) ^ (a5 & b2) ^ (a6 & b1) ^ (a7 & b0);
    uint32_t c8  = (a1 & b7) ^ (a2 & b6) ^ (a3 & b5) ^ (a4 & b4) ^ (a5 & b3) ^ (a6 & b2) ^ (a7 & b1);
    uint32_t c9  = (a2 & b7) ^ (a3 & b6) ^ (a4 & b5) ^ (a5 & b4) ^ (a6 & b3) ^ (a7 & b2);
    uint32_t c10 = (a3 & b7) ^ (a4 & b6) ^ (a5 & b5) ^ (a6 & b4) ^ (a7 & b3);
    uint32_t c11 = (a4 & b7) ^ (a5 & b6) ^ (a6 & b5) ^ (a7 & b4);
    uint32_t c12 = (a5 & b7) ^ (a6 & b6) ^ (a7 & b5);
    uint32_t c13 = (a6 & b7) ^ (a7 & b6);
    uint32_t c14 = a7 & b7;

    return c0 ^ (c1 << 1) ^ (c2 << 2) ^ (c3 << 3) ^ (c4 << 4) ^ (c5 << 5) ^ (c6 << 6) ^ (c7 << 7) ^
           (c8 * 0x1b) ^ (c9 *0x36) ^ (c10 * 0x6c) ^ (c11 * 0xd8) ^ (c12 * 0xab) ^ (c13 * 0x4d) ^ (c14 * 0x9a);
}


static inline uint32_t gf256v_mul_u32_u32(uint32_t a, uint32_t b) {
    uint32_t a0 =  a         & 0x01010101;
    uint32_t a1 = (a >> 1)   & 0x01010101;
    uint32_t a2 = (a >> 2)   & 0x01010101;
    uint32_t a3 = (a >> 3)   & 0x01010101;
    uint32_t a4 = (a >> 4)   & 0x01010101;
    uint32_t a5 = (a >> 5)   & 0x01010101;
    uint32_t a6 = (a >> 6)   & 0x01010101;
    uint32_t a7 = (a >> 7)   & 0x01010101;
    
    uint32_t b0 =  b         & 0x01010101;
    uint32_t b1 = (b >> 1)   & 0x01010101;
    uint32_t b2 = (b >> 2)   & 0x01010101;
    uint32_t b3 = (b >> 3)   & 0x01010101;
    uint32_t b4 = (b >> 4)   & 0x01010101;
    uint32_t b5 = (b >> 5)   & 0x01010101;
    uint32_t b6 = (b >> 6)   & 0x01010101;
    uint32_t b7 = (b >> 7)   & 0x01010101;
    
    return _gf256v_mul_u32_u32(a0, a1, a2, a3, a4, a5, a6, a7,
                                b0, b1, b2, b3, b4, b5, b6, b7);
}


// gf256 := gf2[X]/ (x^8+x^4+x^3+x+1)   // 0x11b , AES field

// z = a * b in arithmetic shares
// O(n+n^2) ~ O(n^2)
void masked_gf256_mul(Masked *z, const Masked a, const Masked b){
    int i, j;
    uint8_t r[2];
    for (i = 0; i < N_SHARES; i++){
        z->shares[i] = gf256_mul(a.shares[i], b.shares[i]); // z_i = a_i * b_i
        #ifdef COUNT
        count_op++;
        #endif
    }
    for (i = 0; i < N_SHARES; i++){
        for (j = i + 1; j < N_SHARES; j++){
            r[0] = rand8(); // sample random r_ij
            r[1] = r[0] ^ (gf256_mul(a.shares[i], b.shares[j])) ^ (gf256_mul(a.shares[j], b.shares[i]));
            z->shares[i] ^= r[0];
            z->shares[j] ^= r[1];
            #ifdef COUNT
            count_op += 7;
            #endif
        }
    }
}

// z = a * b, batch version
// O(n + n(n-1)/2 * ())
// space complexity O(n)
void masked_gf256v_mul(Masked_u32 *z, const Masked_u32 a, const Masked b){
    int i, j;
    uint32_t r[2];
    for (i = 0; i < N_SHARES; i++){
        z->shares[i] = gf256v_mul_u32(a.shares[i], b.shares[i]); // z_i = a_i * b_i
    }
    for (i = 0; i < N_SHARES; i++){
        for (j = i + 1; j < N_SHARES; j++){
            r[0] = rand32(); // sample random r_ij
            r[1] = r[0] ^ (gf256v_mul_u32(a.shares[i], b.shares[j])) ^ (gf256v_mul_u32(a.shares[j], b.shares[i]));
            z->shares[i] ^= r[0];
            z->shares[j] ^= r[1];
        }
    }
}

void masked_gf256v_mul_u32_u32(Masked_u32 *z, const Masked_u32 a, const Masked_u32 b){
    int i, j;
    uint32_t r[2];
    for (i = 0; i < N_SHARES; i++)
        z->shares[i] = gf256v_mul_u32_u32(a.shares[i], b.shares[i]); // z_i = a_i * b_i
    
    for (i = 0; i < N_SHARES - 1; i++){
        for (j = i + 1; j < N_SHARES; j++){
            r[0] = rand32(); // sample random r_ij
            r[1] = r[0] ^ (gf256v_mul_u32_u32(a.shares[i], b.shares[j])) ^ (gf256v_mul_u32_u32(a.shares[j], b.shares[i]));
            z->shares[i] ^= r[0];
            z->shares[j] ^= r[1];
        }
    }
}

// z = a^2 in arithmetic shares
void masked_gf256_sqr(Masked *z, const Masked a){
    Masked fresh_a;
    for (int i = 0; i < N_SHARES; i++)
        fresh_a.shares[i] = a.shares[i];
    gf256_linear_arithmetic_refresh(&fresh_a);
    masked_gf256_mul(z, a, fresh_a);
}

// z = a^2 in arithmetic shares, batch version
void masked_gf256v_sqr(Masked_u32 *z, const Masked_u32 a){
    Masked_u32 fresh_a;
    for (int i = 0; i < N_SHARES; i++)
        fresh_a.shares[i] = a.shares[i];
    gf256_linear_arithmetic_refresh_u32(&fresh_a);
    masked_gf256v_mul_u32_u32(z, a, fresh_a);
}

// z = a^-1 in arithmetic shares, using fermat's little theorem a^-1 = a^254
void masked_gf256_inv(Masked *z, const Masked a){
    // Masked a2 = masked_gf256_sqr(a);
    Masked a2, a4, a4_2, a8_4_2, a8, a64, a64_2, a128;
    Masked t;

    masked_gf256_sqr(&a2, a);                   // a^2
    masked_gf256_sqr(&a4, a2);                 // a^4
    masked_gf256_sqr(&a8, a4);                 // a^8
    masked_gf256_mul(&a4_2, a4, a2);          // a^4 * a^2 = a^6
    masked_gf256_mul(&a8_4_2, a4_2, a8);      // a^6 * a^8 = a^14
    masked_gf256_sqr(&a64, a8_4_2);            // a^14 * a^14 = a^28
    masked_gf256_sqr(&t, a64);                 // a^28 * a^28 = a^56
    masked_gf256_sqr(&a64, t);                 // a^56 * a^56 = a^112
    masked_gf256_mul(&a64_2, a64, a8_4_2);    // a^112 * a^14 = a^126
    masked_gf256_sqr(&a128, a64_2);            // a^126 * a^126 = a^252
    masked_gf256_mul(z, a2, a128);            // z = a^2 * a^252 = a^254 = a^-1
}